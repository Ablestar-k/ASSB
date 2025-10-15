import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

file_path = "../result/ensemble_product.txt"
column_names = [
    'Step', 'Time(ps)', 'Temp', 'Temp_std', 'P', 'P_std',
    'PE', 'PE_std', 'KE', 'KE_std', 'Etot', 'Etot_std',
    'Vol', 'Vol_std', 'Density', 'Density_std'
]

try:
    df = pd.read_csv(
        file_path,
        comment='#',
        delim_whitespace=True,
        header=None,
        names=column_names
    )
    print(f"Successfully loaded data from '{file_path}'.")

except FileNotFoundError:
    print(f"Error: The file '{file_path}' was not found. Please check the file path.")
    exit()

output_path = "../result"
os.makedirs(output_path, exist_ok=True)
plt.style.use('seaborn-v0_8-whitegrid')

# Temperature
avg_temp = df['Temp'].mean()


plt.figure(figsize=(10, 7))
plt.plot(df['Time(ps)'], df['Temp'], color='r', label='Temperature')
plt.fill_between(
    df['Time(ps)'],
    df['Temp'] - df['Temp_std'],
    df['Temp'] + df['Temp_std'],
    color='r',
    alpha=0.2,
    label='Std Dev'
)

#plt.axvline(x=1780, color='blue', linestyle='--', linewidth=1.0)
#plt.axvline(x=1880, color='blue', linestyle='--', linewidth=1.0)
#plt.axvline(x=2080, color='blue', linestyle='--', linewidth=1.0)
#plt.axvline(x=2180, color='blue', linestyle='--', linewidth=1.0)
#plt.axvline(x=2380, color='blue', linestyle='--', linewidth=1.0)
#plt.axvline(x=2480, color='blue', linestyle='--', linewidth=1.0)

plt.axhline(y=df['Temp'][0], color='gray', linestyle='--', linewidth=1.0, label=f'Initial: {df["Temp"][0]:.2f} +- {df["Temp_std"][0]:.2f}K')
plt.axhline(y=avg_temp, color='black', linestyle='--', label=f'Average: {avg_temp:.2f} K')
plt.axhline(y=300, color='green', linestyle='--', label='Target Temp: 300 K')

plt.title('Temperature vs. Time', fontsize=16)
plt.xlabel('Time (ps)', fontsize=12)
plt.ylabel('Temperature (K)', fontsize=12)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(output_path, 'temperature_avg.png'), dpi=300)
plt.close()
print(f"Temperature plot saved to 'temperature_avg.png'. Average Temp: {avg_temp:.2f} K")


# Pressure
avg_pressure = df['P'].mean()

plt.figure(figsize=(10, 7))
plt.plot(df['Time(ps)'], df['P'], color='b', label='Pressure')
plt.fill_between(
    df['Time(ps)'],
    df['P'] - df['P_std'],
    df['P'] + df['P_std'],
    color='b',
    alpha=0.2,
    label='Std Dev'
)
plt.axhline(y=avg_pressure, color='black', linestyle='--', label=f'Average: {avg_pressure:.2f} atm')

plt.title('Pressure vs. Time', fontsize=16)
plt.xlabel('Time (ps)', fontsize=12)
plt.ylabel('Pressure (atm)', fontsize=12)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(output_path, 'pressure_avg.png'), dpi=300)
plt.close()
print(f"Pressure plot saved to 'pressure_avg.png'. Average Pressure: {avg_pressure:.2f} atm")


# Total Energy
avg_energy = df['Etot'].mean()

plt.figure(figsize=(10, 7))
plt.plot(df['Time(ps)'], df['Etot'], color='g', label='Total Energy')
plt.fill_between(
    df['Time(ps)'],
    df['Etot'] - df['Etot_std'],
    df['Etot'] + df['Etot_std'],
    color='g',
    alpha=0.2,
    label='Std Dev'
)
plt.axhline(y=avg_energy, color='black', linestyle='--', label=f'Average: {avg_energy:.2f} eV')

plt.title('Total Energy vs. Time', fontsize=16)
plt.xlabel('Time (ps)', fontsize=12)
plt.ylabel('Total Energy (eV)', fontsize=12)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(output_path, 'energy_avg.png'), dpi=300)
plt.close()
print(f"Total Energy plot saved to 'energy_avg.png'. Average Energy: {avg_energy:.2f} eV")

print("\nAll plots have been saved successfully as separate files.")

# Total PE
avg_PE_energy = df['PE'].mean()

plt.figure(figsize=(10, 7))
plt.plot(df['Time(ps)'], df['PE'], color='g', label='Potential Energy')
plt.fill_between(
    df['Time(ps)'],
    df['PE'] - df['PE_std'],
    df['PE'] + df['PE_std'],
    color='g',
    alpha=0.2,
    label='Std Dev'
)
plt.axhline(y=avg_PE_energy, color='black', linestyle='--', label=f'Average: {avg_PE_energy:.2f} eV')

plt.title('PE vs. Time', fontsize=16)
plt.xlabel('Time (ps)', fontsize=12)
plt.ylabel('PE (eV)', fontsize=12)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(output_path, 'potential_energy_avg.png'), dpi=300)
plt.close()
print(f"Total Energy plot saved to 'potential_energy_avg.png'. Average Energy: {avg_PE_energy:.2f} eV")

print("\nAll plots have been saved successfully as separate files.")

# Kinetic Energy
avg_KE_energy = df['KE'].mean()

plt.figure(figsize=(10, 7))
plt.plot(df['Time(ps)'], df['KE'], color='g', label='Kinetic Energy')
plt.fill_between(
    df['Time(ps)'],
    df['KE'] - df['KE_std'],
    df['KE'] + df['KE_std'],
    color='g',
    alpha=0.2,
    label='Std Dev'
)
plt.axhline(y=avg_KE_energy, color='black', linestyle='--', label=f'Average: {avg_KE_energy:.2f} eV')

plt.title('Kinetic vs. Time', fontsize=16)
plt.xlabel('Time (ps)', fontsize=12)
plt.ylabel('Kinetic (eV)', fontsize=12)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(output_path, 'kinetic_energy_avg.png'), dpi=300)
plt.close()
print(f"Total Energy plot saved to 'kinetic_energy_avg.png'. Average Energy: {avg_KE_energy:.2f} eV")

print("\nAll plots have been saved successfully as separate files.")