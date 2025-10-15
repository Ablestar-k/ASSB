# Set environment variable for PyTorch memory management
import os
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
print(f"--- My script actually sees: CUDA_VISIBLE_DEVICES='{os.getenv('CUDA_VISIBLE_DEVICES')}' ---")
import sys
import torch

NUM_CPU_CORES = 1
torch.set_num_threads(NUM_CPU_CORES)
print(f"--- PyTorch is set to use {torch.get_num_threads()} CPU threads ---")

import time

import numpy as np
from ase import Atoms
from ase import units
from ase.io import read, write, Trajectory
from ase.geometry import wrap_positions
from ase.md.andersen import Andersen
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.langevin import Langevin
from ase.md.nose_hoover_chain import NoseHooverChainNVT
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
fmax_criteria = 0.05 # ev/A

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

# Pre-Production (NVT+NVE) 
# NVE(0.2 ns) -> NVT(0.1 ns) * 3
# To solve temperature problem during production run
#PRE_PROD_NVT_RUN_STETS = 50000 # NVT Production for 0.1 ns
#PRE_PROD_NVE_RUN_STEPS = 100000 # NVE Production for 0.2 ns

# Additional Equilibration (NVT)


# Production(NVT)
PROD_NVT_RUN_STEPS = 5000000 # Production for 10 ns 

# I/O Parameters
THERMO_FREQ = 1000  # Frequency for thermodynamic data (Every 2 ps)
DUMP_FREQ = 1000    # Frequency for trajectory data (Every 2 ps)

# --- File and Path Definitions ---
SIMULATION_NAME = f'NTOC_ver4_{ENSEMBLE_INDEX}'
INITIAL_PATH = f'../initial'
DUMP_PATH = f'../dump/dump_{SIMULATION_NAME}'
THERMO_PATH = f'../thermo/thermo_{SIMULATION_NAME}'
DATA_PATH = f'../data/data_{SIMULATION_NAME}'

#INITIAL_STRUCTURE_FILE = f'../dump/dump_NTOC_ver1_{ENSEMBLE_INDEX}/NTOC_ver1_{ENSEMBLE_INDEX}_product_nve.traj'
RESTART_FILE = f'../dump/dump_NTOC_ver4_{ENSEMBLE_INDEX}/NTOC_ver4_{ENSEMBLE_INDEX}_product_nvt.traj'

# --- ML Potential ---
CHECKPOINT_PATH = "./mattersim-v1.0.0-5M.pth"

# Output files
MINIMIZED_DATA_FILE = f'{DATA_PATH}/{SIMULATION_NAME}_minimized.data'
HEATING_DATA_FILE = f'{DATA_PATH}/{SIMULATION_NAME}_heating.data'
HEATING_EQ_NVT_DATA_FILE = f'{DATA_PATH}/{SIMULATION_NAME}_heating_eq_nvt.data'
QUENCHING_DATA_FILE = f'{DATA_PATH}/{SIMULATION_NAME}_quench.data'
QUENCHING_EQ_DATA_FILE = f'{DATA_PATH}/{SIMULATION_NAME}_quench_eq_nvt.data'
PRE_PRODUCT_DATA_FILE = f'{DATA_PATH}/{SIMULATION_NAME}_pre_product.data'
#PRODUCT_DATA_FILE = f'{DATA_PATH}/{SIMULATION_NAME}_product_nve.data'
PRODUCT_DATA_FILE = f'{DATA_PATH}/{SIMULATION_NAME}_product_nvt.data'

HEATING_LOG_FILE = f'{THERMO_PATH}/{SIMULATION_NAME}_heating.thermo'
HEATING_EQ_NVT_LOG_FILE = f'{THERMO_PATH}/{SIMULATION_NAME}_heating_eq_nvt.thermo'
QUENCHING_LOG_FILE = f'{THERMO_PATH}/{SIMULATION_NAME}_quench.thermo'
QUENCHING_EQ_LOG_FILE = f'{THERMO_PATH}/{SIMULATION_NAME}_quench_eq_nvt.thermo'
PRE_PRODUCT_LOG_FILE = f'{THERMO_PATH}/{SIMULATION_NAME}_pre_product.thermo'
#PRODUCT_LOG_FILE = f'{THERMO_PATH}/{SIMULATION_NAME}_product_nve.thermo'
PRODUCT_LOG_FILE = f'{THERMO_PATH}/{SIMULATION_NAME}_product_nvt.thermo'

HEATING_TRAJECTORY_FILE = f'{DUMP_PATH}/{SIMULATION_NAME}_heating.traj'
HEATING_EQ_NVT_TRAJECTORY_FILE = f'{DUMP_PATH}/{SIMULATION_NAME}_heating_eq_nvt.traj'
QUENCHING_TRAJECTORY_FILE = f'{DUMP_PATH}/{SIMULATION_NAME}_quench.traj'
QUENCHING_EQ_TRAJECTORY_FILE = f'{DUMP_PATH}/{SIMULATION_NAME}_quench_eq_nvt.traj'
PRE_PRODUCT_TRAJECTORY_FILE = f'{DUMP_PATH}/{SIMULATION_NAME}_pre_product.traj' 
#PRODUCT_TRAJECTORY_FILE = f'{DUMP_PATH}/{SIMULATION_NAME}_product_nve.traj'
PRODUCT_TRAJECTORY_FILE = f'{DUMP_PATH}/{SIMULATION_NAME}_product_nvt.traj'

BENCHMARK_LOG_FILE = f'{THERMO_PATH}/{SIMULATION_NAME}_benchmark.log'

# --- Hardware Setup ---
if torch.cuda.is_available():
    device = 'cuda'
    print(f"CUDA found. Using GPU: {torch.cuda.get_device_name(0)}")
else:
    device = 'cpu'
    print("CUDA not found. Using CPU.")

# --- Set up log file ---
# For checking Logging and Thermodynamic Data
def setup_thermo_logger(atoms_obj, logfile_handle, thermo_freq, start_step=0):

    if logfile_handle.tell() == 0:
        logfile_handle.write(f"{'Step':>12s} {'Time[ps]':>12s} {'Temp[K]':>10s} "
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

            logfile_handle.write(f"{step:>12d} {time_ps:12.4f} {temp:10.2f} "
                                 f"{pressure_atm:12.4f} {pe:15.6f} {ke:15.6f} "
                                 f"{(pe+ke):15.6f} {vol:15.3f} {density:18.6f}\n")
        
        step_counter[0] += 1

    return log_data

# For benchmarking simulation performance
def log_benchmark(stage_name, steps, elapsed_seconds):
    """Calculates and logs the simulation performance."""
    if elapsed_seconds <= 1e-6: 
        return

    sim_time_ps = steps * TIMESTEP_FS / 1000.0

    ps_per_day = (sim_time_ps / elapsed_seconds) * 86400 

    if not os.path.exists(os.path.dirname(BENCHMARK_LOG_FILE)):
        os.makedirs(os.path.dirname(BENCHMARK_LOG_FILE), exist_ok=True)

    if not os.path.exists(BENCHMARK_LOG_FILE):
        with open(BENCHMARK_LOG_FILE, 'w') as f:
            f.write(f"{'Simulation Stage':<35s} {'Wall Time[s]':>18s} {'Performance[ps/day]':>25s}\n")
            f.write(f"{'-'*35:<35s} {'-'*18:>18s} {'-'*25:>25s}\n")

    with open(BENCHMARK_LOG_FILE, 'a') as f:
        f.write(f"{stage_name:<35s} {elapsed_seconds:>18.2f} {ps_per_day:>25.2f}\n")

    print(f"--- Benchmark for {stage_name}: {ps_per_day:.2f} ps/day ---")


def get_last_step_from_log(logfile_path):
    if not os.path.exists(logfile_path) or os.path.getsize(logfile_path) == 0:
        return None

    try:
        with open(logfile_path, 'r') as f:
            lines = f.readlines()
    except IOError:
        return None 

    for line in reversed(lines):
        line = line.strip()
        if not line:  
            continue
        
        try:
            last_step = int(line.split()[0])
            return last_step 
        except (ValueError, IndexError):
            continue
            
    return None 


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

print(f"Reading initial structure from: '{RESTART_FILE}'")
if not os.path.exists(RESTART_FILE):
    raise FileNotFoundError(f"Error: The initial structure file '{RESTART_FILE}' was not found.")
# index = -1 : restart from the last frame
unit_cell = read(RESTART_FILE, index = -1)

print(f"Successfully loaded {len(unit_cell)} atoms from the unit cell file.")

# box_size = [48.78, 48.78, 48.78]
# unit_cell.set_cell(box_size)

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

#print(f"\n--- Starting Energy Minimization (Steps={MIN_STEPS}, f_max={fmax_criteria:.2e} eV/A) ---")
#dyn_opt = FIRE(atoms, maxstep=0.2, logfile=f'{THERMO_PATH}/{SIMULATION_NAME}_minimization.log')
#dyn_opt.run(fmax=fmax_criteria, steps=MIN_STEPS)
#write(MINIMIZED_DATA_FILE, atoms, format='lammps-data')
#print("--- Energy Minimization Finished ---")
#print(f"Final box dimensions after minimization: {atoms.get_cell().lengths()}")
#
## ==========================
## === 4. Heating Process ===
## ==========================
#
## --- 4.1. Gradual Heating (NVT)---
#total_steps_so_far = 0
#
#print(f"\n--- Starting Gradual Heating (NVT) ---")
#MaxwellBoltzmannDistribution(atoms, temperature_K=HEATING_START_TEMP)
#print(f"Initial velocities set to {HEATING_START_TEMP} K.")
#
#print(f"\n--- Heating gradually: {HEATING_START_TEMP} K -> {HEATING_TARGET_TEMP} K for {HEATING_STEPS * TIMESTEP_FS / 1000:.1f} ps ---")
#temp_increment = 50.0
#total_temp_diff = HEATING_TARGET_TEMP - HEATING_START_TEMP
#num_stages = int(total_temp_diff / temp_increment)
#steps_per_stage = HEATING_STEPS // num_stages
#
#thermo_heat_logger = setup_thermo_logger(atoms, HEATING_LOG_FILE, THERMO_FREQ, start_step=total_steps_so_far)
#traj_heat = Trajectory(HEATING_TRAJECTORY_FILE, 'w', atoms)
#
#start_time_heating = time.time()
#
#current_temp = HEATING_START_TEMP
#for i in range(num_stages):
#    target_temp_stage = current_temp + temp_increment
#    print(f"  Heating stage {i+1}/{num_stages}: {current_temp:.0f} K -> {target_temp_stage:.0f} K")
#
#    dyn_heat = Langevin(
#        atoms,
#        timestep=TIMESTEP_FS * units.fs,
#        temperature_K=target_temp_stage,
#        friction=0.02,
#    )
#    
#    dyn_heat.attach(thermo_heat_logger, interval=1)
#    dyn_heat.attach(traj_heat.write, interval=int(DUMP_FREQ))
#    dyn_heat.run(steps_per_stage)
#    
#    current_temp = target_temp_stage
#
#end_time_heating = time.time()
#log_benchmark('4.1 Gradual Heating (NVT)', HEATING_STEPS, end_time_heating - start_time_heating)
#
#total_steps_so_far += HEATING_STEPS
#write(HEATING_DATA_FILE, atoms, format='lammps-data')
#traj_heat.close()
#print("--- Gradual Heating Finished ---")
#
## --- 4.2. Equilibration ---
## --- 4.2.1. Equilibration (NVT) ---
#
#print(f"\n--- Starting Equilibration (NVT) at {HEATING_TARGET_TEMP} K for {HEATING_EQ_NVT_STEPS * TIMESTEP_FS / 1000:.1f} ps ---")
#dyn_eq_nvt = Langevin(
#    atoms,
#    timestep=TIMESTEP_FS * units.fs,
#    temperature_K=HEATING_TARGET_TEMP,
#    friction=0.02,
#)
#
#thermo_eq_nvt_logger = setup_thermo_logger(atoms, HEATING_EQ_NVT_LOG_FILE, THERMO_FREQ, start_step=total_steps_so_far)
#traj_eq_nvt = Trajectory(HEATING_EQ_NVT_TRAJECTORY_FILE, 'w', atoms)
#dyn_eq_nvt.attach(thermo_eq_nvt_logger, interval=1)
#dyn_eq_nvt.attach(traj_eq_nvt.write, interval=int(DUMP_FREQ))
#
#start_time_eq_nvt = time.time()
#dyn_eq_nvt.run(HEATING_EQ_NVT_STEPS)
#end_time_eq_nvt = time.time()
#log_benchmark('4.2 Heating Equilibration (NVT)', HEATING_EQ_NVT_STEPS, end_time_eq_nvt - start_time_eq_nvt)
#
#total_steps_so_far += HEATING_EQ_NVT_STEPS
#write(HEATING_EQ_NVT_DATA_FILE, atoms, format='lammps-data')
#traj_eq_nvt.close()
#print("--- NVT Equilibration Finished ---")
#
## ==============================
## === 5. Quenching (NVT) ===
## ==============================
#
## --- 5.1 Quenching ---
#print(f"\n--- Starting Quenching: {HEATING_TARGET_TEMP} K -> {QUENCHING_TARGET_TEMP} K ---")
#start_temp_quench = HEATING_TARGET_TEMP
#end_temp_quench = QUENCHING_TARGET_TEMP
#temp_decrement_quench = 50.0
#steps_per_stage_quench = int((temp_decrement_quench / QUENCHING_RATE) * 1000 / TIMESTEP_FS)
#num_stages_quench = int((start_temp_quench - end_temp_quench) / temp_decrement_quench)
#
#thermo_quench_logger = setup_thermo_logger(atoms, QUENCHING_LOG_FILE, THERMO_FREQ, start_step=total_steps_so_far)
#traj_quench = Trajectory(QUENCHING_TRAJECTORY_FILE, 'w', atoms)
#
#start_time_quench = time.time()
#
#current_temp_quench = start_temp_quench
#total_quench_steps = 0
#for i in range(num_stages_quench):
#    target_temp_quench = current_temp_quench - temp_decrement_quench
#    print(f"  Quenching stage {i+1}/{num_stages_quench}: {current_temp_quench:.0f} K -> {target_temp_quench:.0f} K...")
#
#    dyn_quench = Langevin(
#        atoms,
#        timestep=TIMESTEP_FS * units.fs,
#        temperature_K=target_temp_quench,
#        friction=0.02
#    )
#    
#    dyn_quench.attach(thermo_quench_logger, interval=1)
#    dyn_quench.attach(traj_quench.write, interval=int(DUMP_FREQ))
#    dyn_quench.run(steps_per_stage_quench)
#    
#    current_temp_quench = target_temp_quench
#    total_quench_steps += steps_per_stage_quench
#
#end_time_quench = time.time()
#log_benchmark('5.1 Quenching (NVT)', total_quench_steps, end_time_quench - start_time_quench)
#
#total_steps_so_far += total_quench_steps
#write(QUENCHING_DATA_FILE, atoms, format='lammps-data')
#traj_quench.close()
#print("--- Quenching Finished ---")
#
## ---- 5.2 Quenching Equilibration(NVT) ----
#print(f"\n--- Starting Quenching Equilibration (NVT) at {QUENCHING_TARGET_TEMP} K for {QUENCHING_EQ_NVT_STEPS * TIMESTEP_FS / 1000:.1f} ps ---")
#dyn_quench_eq = Langevin(
#    atoms,
#    timestep=TIMESTEP_FS * units.fs,
#    temperature_K=QUENCHING_TARGET_TEMP,
#    friction=0.02,
#)
#
#thermo_quench_eq_logger = setup_thermo_logger(atoms, QUENCHING_EQ_LOG_FILE, THERMO_FREQ, start_step=total_steps_so_far)
#traj_quench_eq = Trajectory(QUENCHING_EQ_TRAJECTORY_FILE, 'w', atoms)
#dyn_quench_eq.attach(thermo_quench_eq_logger, interval=1)
#dyn_quench_eq.attach(traj_quench_eq.write, interval=int(DUMP_FREQ))
#
#start_time_quench_eq = time.time()
#dyn_quench_eq.run(QUENCHING_EQ_NVT_STEPS)
#end_time_quench_eq = time.time()
#log_benchmark('5.2 Quenching Equilibration (NVT)', QUENCHING_EQ_NVT_STEPS, end_time_quench_eq - start_time_quench_eq)
#
#total_steps_so_far += QUENCHING_EQ_NVT_STEPS
#write(QUENCHING_EQ_DATA_FILE, atoms, format='lammps-data')
#traj_quench_eq.close()
#print("--- Quenching NVT Equilibration Finished ---")

# ==============================================================================
# === 6. Pre-Production Run (NVT+NVE) ===
# ==============================================================================
#total_steps_so_far = 790000
#
#print(f"\n Apply Andersen thermostat for 3 times to stablized temperature.")
#
#traj_pre_prod = Trajectory(PRE_PRODUCT_TRAJECTORY_FILE, 'w', atoms)
#with open(PRE_PRODUCT_LOG_FILE, 'w') as logfile:
#    for i in range(3):
#        print(f"\n Start {i+1}th NVE")
#        dyn_pre_nve_prod = VelocityVerlet(
#            atoms,
#            timestep=TIMESTEP_FS * units.fs,
#        )
#
#        thermo_pre_prod_logger = setup_thermo_logger(atoms, logfile, THERMO_FREQ, start_step=total_steps_so_far)
#        dyn_pre_nve_prod.attach(thermo_pre_prod_logger, interval=1)
#        dyn_pre_nve_prod.attach(traj_pre_prod.write, interval=int(DUMP_FREQ))
#
#        start_time_prod = time.time()
#        dyn_pre_nve_prod.run(PRE_PROD_NVE_RUN_STEPS)
#        end_time_prod = time.time()
#        log_benchmark(f'6-{i+1}. NVE Pre_Production', PRE_PROD_NVE_RUN_STEPS, end_time_prod - start_time_prod)
#
#        total_steps_so_far += PRE_PROD_NVE_RUN_STEPS
#
#        # -------------------------------------------------------------------
#        
#        print(f"\n Start {i+1}th NVT")
#        dyn_pre_nvt_prod = NVTBerendsen(
#            atoms, timestep = TIMESTEP_FS * units.fs,
#            temperature_K = QUENCHING_TARGET_TEMP,
#            taut = 100 * TIMESTEP_FS * units.fs
#        )
#
#        thermo_pre_prod_logger = setup_thermo_logger(atoms, logfile, THERMO_FREQ, start_step=total_steps_so_far)
#        dyn_pre_nvt_prod.attach(thermo_pre_prod_logger, interval=1)
#        dyn_pre_nvt_prod.attach(traj_pre_prod.write, interval=int(DUMP_FREQ))
#
#        start_time_quench_eq = time.time()
#        dyn_pre_nvt_prod.run(PRE_PROD_NVT_RUN_STETS)
#        end_time_quench_eq = time.time()
#        log_benchmark(f'6-{i+1}. NVT Pre_Production', PRE_PROD_NVT_RUN_STETS, end_time_quench_eq - start_time_quench_eq)
#
#        total_steps_so_far += PRE_PROD_NVT_RUN_STETS
#
#traj_pre_prod.close()
#write(PRE_PRODUCT_DATA_FILE, atoms, format='lammps-data')
#
#print("--- Pre Production Run Finished. ---")


# =======================================
# === 7. Production Run (NVT) ===
# =======================================
total_steps_so_far = 790000

last_step_found = get_last_step_from_log(PRODUCT_LOG_FILE)

if last_step_found is not None:
    total_steps_so_far = last_step_found
    print(f"Log file found. Restarting simulation from step {total_steps_so_far}.")
    atoms = read(PRODUCT_TRAJECTORY_FILE, index=-1)
    atoms.calc = MatterSimCalculator(potential=potential) 
else:
    #total_steps_so_far = PRE_PROD_NVT_RUN_STEPS
    total_steps_so_far = 790000
    print(f"No valid log found. Starting new production run from step {total_steps_so_far}.")


print(f"\n--- Starting Production Run (NVT) at {QUENCHING_TARGET_TEMP} K for {PROD_NVT_RUN_STEPS * TIMESTEP_FS / 1000:.1f} ps ---")
dyn_prod = NoseHooverChainNVT(
    atoms,
    timestep=TIMESTEP_FS * units.fs,
    temperature_K = QUENCHING_TARGET_TEMP,
    tdamp = 200 * units.fs,
    tchain = 3,
    tloop = 2
)


traj_prod = Trajectory(PRODUCT_TRAJECTORY_FILE, 'a', atoms)
with open(PRODUCT_LOG_FILE, 'a') as logfile:
    thermo_prod_logger = setup_thermo_logger(atoms, logfile, THERMO_FREQ, start_step=total_steps_so_far)
    dyn_prod.attach(thermo_prod_logger, interval=1)
    dyn_prod.attach(traj_prod.write, interval=int(DUMP_FREQ))

    start_time_prod = time.time()
    dyn_prod.run(PROD_NVT_RUN_STEPS)
    end_time_prod = time.time()
    log_benchmark('7. Production (NVT)', PROD_NVT_RUN_STEPS, end_time_prod - start_time_prod)

total_steps_so_far += PROD_NVT_RUN_STEPS
traj_prod.close()
write(PRODUCT_DATA_FILE, atoms, format='lammps-data')

print("--- Production Run Finished. Simulation Complete. ---")

