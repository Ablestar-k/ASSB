import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

file_path = "../result/ensemble_avg.txt"
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

plt.style.use('seaborn-v0_8-whitegrid')

fig, axes = plt.subplots(3, 1, figsize=(12, 18), sharex=True)

# 1) Temperature vs. Time plot
axes[0].plot(df['Time(ps)'], df['Temp'], color='r', label='Average Temperature')
axes[0].fill_between(
    df['Time(ps)'],
    df['Temp'] - df['Temp_std'],
    df['Temp'] + df['Temp_std'],
    color='r',
    alpha=0.2,  
    label='Temperature Std Dev'
)
axes[0].set_title('Temperature vs. Time', fontsize=16)
axes[0].set_ylabel('Temperature (K)', fontsize=12)
axes[0].legend()
axes[0].grid(True)

# 2) Pressure vs. Time plot
axes[1].plot(df['Time(ps)'], df['P'], color='b', label='Average Pressure')
axes[1].fill_between(
    df['Time(ps)'],
    df['P'] - df['P_std'],
    df['P'] + df['P_std'],
    color='b',
    alpha=0.2,  
    label='Pressure Std Dev'
)
axes[1].set_title('Pressure vs. Time', fontsize=16)
axes[1].set_ylabel('Pressure (bar)', fontsize=12)
axes[1].legend()
axes[1].grid(True)

# 3) Total Energy vs. Time plot
axes[2].plot(df['Time(ps)'], df['Etot'], color='g', label='Average Total Energy')
axes[2].fill_between(
    df['Time(ps)'],
    df['Etot'] - df['Etot_std'],
    df['Etot'] + df['Etot_std'],
    color='g',
    alpha=0.2,  
    label='Total Energy Std Dev'
)
axes[2].set_title('Total Energy vs. Time', fontsize=16)
axes[2].set_xlabel('Time (ps)', fontsize=12)  
axes[2].set_ylabel('Total Energy (eV)', fontsize=12)
axes[2].legend()
axes[2].grid(True)


plt.tight_layout()

plt.savefig('simulation_overview_with_std.png')

print("\nThe combined plot has been saved as 'simulation_overview_with_std.png'.")