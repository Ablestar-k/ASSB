# Set environment variable for PyTorch memory management
import os
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
import sys
import torch

import numpy as np
from ase import Atoms
from ase import units
from ase.io import read, write, Trajectory
from ase.geometry import wrap_positions
from ase.md.langevin import Langevin
#from ase.md.npt import NPT 
from ase.md.verlet import VelocityVerlet # For NVE
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.optimize import FIRE
from mattersim.forcefield.potential import Potential
from mattersim.forcefield import MatterSimCalculator

# --- Simulation Parameters ---

ENSEMBLE_INDEX = sys.argv[1]

# General
TIMESTEP_FS = 2.0

# Energy minimazation parameter
MIN_STEPS = 50000
fmax_criteria = 0.005 # ev/A

# Heating parameter (NVT)
TARGET_PRESSURE = 1.0 # atm
ATM_TO_GPA = 101325 * 1e-9    # GPa/atm
HEATING_START_TEMP = 300.0 # Kelvin
HEATING_TARGET_TEMP = 1200.0
HEATING_STEPS = 200000 # Heating for 200 ps

# Equilibration parameter (NVT)
HEATING_EQ_NVT_STEPS = 250000 # Equilibration at NVT for 0.5 ns to stabilize pressure

# Quenching parameter(NVT)
QUENCHING_TARGET_TEMP = 300.0
QUENCHING_RATE = 5 # K/ps

# Quenching Equilibration (NVT) - Changed from NPT to NVT
QUENCHING_EQ_NVT_STEPS = 250000 # Production for 0.5 ns

# Production (NVE) - Changed from NVT to NVE
PROD_NVE_RUN_STEPS = 500000 # Production for 1.0 ns

# I/O Parameters
THERMO_FREQ = 1000  # Frequency for thermodynamic data (Every 2 ps)
DUMP_FREQ = 1000    # Frequency for trajectory data (Every 2 ps)

# --- File and Path Definitions ---
SIMULATION_NAME = f'NTOC_ver1_{ENSEMBLE_INDEX}'
INITIAL_PATH = f'../initial'
DUMP_PATH = f'../dump/dump_{SIMULATION_NAME}'
THERMO_PATH = f'../thermo/thermo_{SIMULATION_NAME}'
DATA_PATH = f'../data/data_{SIMULATION_NAME}'

INITIAL_STRUCTURE_FILE = f'{INITIAL_PATH}/NTOC_{ENSEMBLE_INDEX}.xyz'

# --- ML Potential ---
CHECKPOINT_PATH = "./mattersim-v1.0.0-5M.pth"

# Output files
MINIMIZED_DATA_FILE = f'{DATA_PATH}/{SIMULATION_NAME}_minimized.data'
HEATING_DATA_FILE = f'{DATA_PATH}/{SIMULATION_NAME}_heating.data'
HEATING_EQ_NVT_DATA_FILE = f'{DATA_PATH}/{SIMULATION_NAME}_heating_eq_nvt.data'
QUENCHING_DATA_FILE = f'{DATA_PATH}/{SIMULATION_NAME}_quench.data'
QUENCHING_EQ_DATA_FILE = f'{DATA_PATH}/{SIMULATION_NAME}_quench_eq_nvt.data' 
PRODUCT_DATA_FILE = f'{DATA_PATH}/{SIMULATION_NAME}_product_nve.data' 

HEATING_LOG_FILE = f'{THERMO_PATH}/{SIMULATION_NAME}_heating.thermo'
HEATING_EQ_NVT_LOG_FILE = f'{THERMO_PATH}/{SIMULATION_NAME}_heating_eq_nvt.thermo'
QUENCHING_LOG_FILE = f'{THERMO_PATH}/{SIMULATION_NAME}_quench.thermo'
QUENCHING_EQ_LOG_FILE = f'{THERMO_PATH}/{SIMULATION_NAME}_quench_eq_nvt.thermo' 
PRODUCT_LOG_FILE = f'{THERMO_PATH}/{SIMULATION_NAME}_product_nve.thermo'

HEATING_TRAJECTORY_FILE = f'{DUMP_PATH}/{SIMULATION_NAME}_heating.traj'
HEATING_EQ_NVT_TRAJECTORY_FILE = f'{DUMP_PATH}/{SIMULATION_NAME}_heating_eq_nvt.traj'
QUENCHING_TRAJECTORY_FILE = f'{DUMP_PATH}/{SIMULATION_NAME}_quench.traj'
QUENCHING_EQ_TRAJECTORY_FILE = f'{DUMP_PATH}/{SIMULATION_NAME}_quench_eq_nvt.traj'
PRODUCT_TRAJECTORY_FILE = f'{DUMP_PATH}/{SIMULATION_NAME}_product_nve.traj'

# --- Hardware Setup ---
if torch.cuda.is_available():
    device = 'cuda'
    print(f"CUDA found. Using GPU: {torch.cuda.get_device_name(0)}")
else:
    device = 'cpu'
    print("CUDA not found. Using CPU.")


def setup_thermo_logger(atoms_obj, filename, thermo_freq, start_step=0):
    """Sets up a logger to print thermodynamic data during a simulation."""
    if not os.path.exists(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename), exist_ok=True)
    if start_step == 0 and os.path.exists(filename):
        os.remove(filename)

    with open(filename, 'a') as f:
        if f.tell() == 0:
            f.write(f"{'Step':>12s} {'Time[ps]':>12s} {'Temp[K]':>10s} "
                    f"{'Press[atm]':>12s} {'PotEng[eV]':>15s} "
                    f"{'KinEng[eV]':>15s} {'TotEng[eV]':>15s} "
                    f"{'Volume[A^3]':>15s} {'Density[g/cm^3]':>18s}\n")

    step_counter = [start_step]

    def log_data():
        step = step_counter[0]
        if step % int(thermo_freq) == 0:
            pe = atoms_obj.get_potential_energy()
            ke = atoms_obj.get_kinetic_energy()
            temp = ke / (1.5 * units.kB * len(atoms_obj))
            time_ps = step * TIMESTEP_FS / 1000.0
            vol = atoms_obj.get_volume()
            density = sum(atoms_obj.get_masses()) / vol * 1.66054

            try:
                stress = atoms_obj.get_stress(voigt=False)
                pressure_gpa = -np.trace(stress) / 3.0 / units.GPa
                pressure_atm = pressure_gpa * 9869.23
            except Exception:
                pressure_atm = 0.0

            with open(filename, 'a') as f:
                f.write(f"{step:>12d} {time_ps:12.4f} {temp:10.2f} "
                        f"{pressure_atm:12.4f} {pe:15.6f} {ke:15.6f} "
                        f"{(pe+ke):15.6f} {vol:15.3f} {density:18.6f}\n")
        
        step_counter[0] += 1

    return log_data

# ==============================================================================
# ==============================================================================
# Main Script
# ==============================================================================
# ==============================================================================

# =======================================================
# === 1. Setup Directories and Load Initial Structure ===
# =======================================================

print("--- Setting up simulation environment ---")
for path in [DUMP_PATH, THERMO_PATH, DATA_PATH]:
    os.makedirs(path, exist_ok=True)

print(f"Reading initial structure from: '{INITIAL_STRUCTURE_FILE}'")
if not os.path.exists(INITIAL_STRUCTURE_FILE):
    raise FileNotFoundError(f"Error: The initial structure file '{INITIAL_STRUCTURE_FILE}' was not found.")

unit_cell = read(INITIAL_STRUCTURE_FILE, format='xyz')
print(f"Successfully loaded {len(unit_cell)} atoms from the unit cell file.")

box_size = [48.78, 48.78, 48.78]
unit_cell.set_cell(box_size)
unit_cell.set_pbc(True)
unit_cell.set_positions(
    wrap_positions(unit_cell.get_positions(), unit_cell.get_cell(), pbc=[True]*3)
)
atoms = unit_cell

print("\n--- Supercell Information ---")
print(f"Total atoms: {len(atoms)}")
print("Supercell vectors:\n", atoms.get_cell())

# ============================
# === 2. Load ML Potential ===
# ============================

print(f"\n--- Loading ML potential ---")
print(f"Checkpoint path: {CHECKPOINT_PATH}")
if not os.path.exists(CHECKPOINT_PATH):
    raise FileNotFoundError(f"Error: Checkpoint file not found at '{CHECKPOINT_PATH}'.")
torch.cuda.empty_cache()
potential = Potential.from_checkpoint(CHECKPOINT_PATH, device=device)
atoms.calc = MatterSimCalculator(potential=potential)
print(f"Successfully attached potential to device: '{device}'")


# ===============================
# === 3. Energy Minimization ====
# ===============================

print(f"\n--- Starting Energy Minimization (Steps={MIN_STEPS}, f_max={fmax_criteria:.2e} eV/A) ---")
dyn_opt = FIRE(atoms, maxstep=0.2, logfile=f'{THERMO_PATH}/{SIMULATION_NAME}_minimization.log')
dyn_opt.run(fmax=fmax_criteria, steps=MIN_STEPS)
write(MINIMIZED_DATA_FILE, atoms, format='lammps-data')
print("--- Energy Minimization Finished ---")
print(f"Final box dimensions after minimization: {atoms.get_cell().lengths()}")

# ==========================
# === 4. Heating Process ===
# ==========================

# --- 4.1. Gradual Heating (NVT)---
total_steps_so_far = 0

print(f"\n--- Starting Gradual Heating (NVT) ---")
MaxwellBoltzmannDistribution(atoms, temperature_K=HEATING_START_TEMP)
print(f"Initial velocities set to {HEATING_START_TEMP} K.")

print(f"\n--- Heating gradually: {HEATING_START_TEMP} K -> {HEATING_TARGET_TEMP} K for {HEATING_STEPS * TIMESTEP_FS / 1000:.1f} ps ---")
temp_increment = 50.0
total_temp_diff = HEATING_TARGET_TEMP - HEATING_START_TEMP
num_stages = int(total_temp_diff / temp_increment)
steps_per_stage = HEATING_STEPS // num_stages

thermo_heat_logger = setup_thermo_logger(atoms, HEATING_LOG_FILE, THERMO_FREQ, start_step=total_steps_so_far)
traj_heat = Trajectory(HEATING_TRAJECTORY_FILE, 'w', atoms)

current_temp = HEATING_START_TEMP
for i in range(num_stages):
    target_temp_stage = current_temp + temp_increment
    print(f"  Heating stage {i+1}/{num_stages}: {current_temp:.0f} K -> {target_temp_stage:.0f} K")

    dyn_heat = Langevin(
        atoms,
        timestep=TIMESTEP_FS * units.fs,
        temperature_K=target_temp_stage,
        friction=0.02,
    )
    
    dyn_heat.attach(thermo_heat_logger, interval=1)
    dyn_heat.attach(traj_heat.write, interval=int(DUMP_FREQ))
    dyn_heat.run(steps_per_stage)
    
    current_temp = target_temp_stage


total_steps_so_far += HEATING_STEPS
write(HEATING_DATA_FILE, atoms, format='lammps-data')
traj_heat.close()
print("--- Gradual Heating Finished ---")

# --- 4.2. Equilibration ---
# --- 4.2.1. Equilibration (NVT) ---

print(f"\n--- Starting Equilibration (NVT) at {HEATING_TARGET_TEMP} K for {HEATING_EQ_NVT_STEPS * TIMESTEP_FS / 1000:.1f} ps ---")
dyn_eq_nvt = Langevin(
    atoms,
    timestep=TIMESTEP_FS * units.fs,
    temperature_K=HEATING_TARGET_TEMP,
    friction=0.02,
)

thermo_eq_nvt_logger = setup_thermo_logger(atoms, HEATING_EQ_NVT_LOG_FILE, THERMO_FREQ, start_step=total_steps_so_far)
traj_eq_nvt = Trajectory(HEATING_EQ_NVT_TRAJECTORY_FILE, 'w', atoms)
dyn_eq_nvt.attach(thermo_eq_nvt_logger, interval=1)
dyn_eq_nvt.attach(traj_eq_nvt.write, interval=int(DUMP_FREQ))
dyn_eq_nvt.run(HEATING_EQ_NVT_STEPS)

total_steps_so_far += HEATING_EQ_NVT_STEPS
write(HEATING_EQ_NVT_DATA_FILE, atoms, format='lammps-data')
traj_eq_nvt.close()
print("--- NVT Equilibration Finished ---")

# ==============================
# === 5. Quenching (NVT) ===
# ==============================

# --- 5.1 Quenching ---
print(f"\n--- Starting Quenching: {HEATING_TARGET_TEMP} K -> {QUENCHING_TARGET_TEMP} K ---")
start_temp_quench = HEATING_TARGET_TEMP
end_temp_quench = QUENCHING_TARGET_TEMP
temp_decrement_quench = 50.0
steps_per_stage_quench = int((temp_decrement_quench / QUENCHING_RATE) * 1000 / TIMESTEP_FS)
num_stages_quench = int((start_temp_quench - end_temp_quench) / temp_decrement_quench)

thermo_quench_logger = setup_thermo_logger(atoms, QUENCHING_LOG_FILE, THERMO_FREQ, start_step=total_steps_so_far)
traj_quench = Trajectory(QUENCHING_TRAJECTORY_FILE, 'w', atoms)

current_temp_quench = start_temp_quench
total_quench_steps = 0
for i in range(num_stages_quench):
    target_temp_quench = current_temp_quench - temp_decrement_quench
    print(f"  Quenching stage {i+1}/{num_stages_quench}: {current_temp_quench:.0f} K -> {target_temp_quench:.0f} K...")

    dyn_quench = Langevin(
        atoms,
        timestep=TIMESTEP_FS * units.fs,
        temperature_K=target_temp_quench,
        friction=0.02
    )
    
    dyn_quench.attach(thermo_quench_logger, interval=1)
    dyn_quench.attach(traj_quench.write, interval=int(DUMP_FREQ))
    dyn_quench.run(steps_per_stage_quench)
    
    current_temp_quench = target_temp_quench
    total_quench_steps += steps_per_stage_quench

total_steps_so_far += total_quench_steps
write(QUENCHING_DATA_FILE, atoms, format='lammps-data')
traj_quench.close()
print("--- Quenching Finished ---")

# ---- 5.2 Quenching Equilibration(NVT) ----
# This section was changed from NPT to NVT to keep the density constant.
print(f"\n--- Starting Quenching Equilibration (NVT) at {QUENCHING_TARGET_TEMP} K for {QUENCHING_EQ_NVT_STEPS * TIMESTEP_FS / 1000:.1f} ps ---")
dyn_quench_eq = Langevin(
    atoms,
    timestep=TIMESTEP_FS * units.fs,
    temperature_K=QUENCHING_TARGET_TEMP,
    friction=0.02,
)

thermo_quench_eq_logger = setup_thermo_logger(atoms, QUENCHING_EQ_LOG_FILE, THERMO_FREQ, start_step=total_steps_so_far)
traj_quench_eq = Trajectory(QUENCHING_EQ_TRAJECTORY_FILE, 'w', atoms)
dyn_quench_eq.attach(thermo_quench_eq_logger, interval=1)
dyn_quench_eq.attach(traj_quench_eq.write, interval=int(DUMP_FREQ))
dyn_quench_eq.run(QUENCHING_EQ_NVT_STEPS)

total_steps_so_far += QUENCHING_EQ_NVT_STEPS
write(QUENCHING_EQ_DATA_FILE, atoms, format='lammps-data')
traj_quench_eq.close()
print("--- Quenching NVT Equilibration Finished ---")

# ===============================
# === 6. Production Run (NVE) ===
# ===============================
# This section was changed to NVE to check for energy conservation in an isolated system.
print(f"\n--- Starting Production Run (NVE) at {QUENCHING_TARGET_TEMP} K for {PROD_NVE_RUN_STEPS * TIMESTEP_FS / 1000:.1f} ps ---")
dyn_prod = VelocityVerlet(
    atoms,
    timestep=TIMESTEP_FS * units.fs,
)

thermo_prod_logger = setup_thermo_logger(atoms, PRODUCT_LOG_FILE, THERMO_FREQ, start_step=total_steps_so_far)
traj_prod = Trajectory(PRODUCT_TRAJECTORY_FILE, 'w', atoms)
dyn_prod.attach(thermo_prod_logger, interval=1)
dyn_prod.attach(traj_prod.write, interval=int(DUMP_FREQ))

dyn_prod.run(PROD_NVE_RUN_STEPS)

total_steps_so_far += PROD_NVE_RUN_STEPS
write(PRODUCT_DATA_FILE, atoms, format='lammps-data')
traj_prod.close()
print("--- Production Run Finished. Simulation Complete. ---")