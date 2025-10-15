import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

target_species = 'Na'

file_path = f"../../result/{target_species}_MSD_product_ensemble_average_mda.dat"

column_names = [
    'Time(ps)', 'MSD_mean', 'MSD_std'
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

output_path = "../../result"

plt.style.use('seaborn-v0_8-whitegrid')

plt.figure(figsize=(10, 7))
plt.plot(df['Time(ps)'], df['MSD_mean'], color='r', label='MSD')
plt.fill_between(
    df['Time(ps)'],
    df['MSD_mean'] - df['MSD_std'],
    df['MSD_mean'] + df['MSD_std'],
    color='r',
    alpha=0.2,
    label='Std Dev'
)

plt.title(f'{target_species} Time vs MSD during production', fontsize=16)
plt.xlabel('Time (ps)', fontsize=12)
plt.ylabel('MSD', fontsize=12)
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(output_path, f'{target_species}_MSD_mda_avg.png'), dpi=300)

plt.close()
