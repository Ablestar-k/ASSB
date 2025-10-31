import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

ATOM_NAME = 'Na'
INPUT_FILE_PATH = "../../result"
INPUT_FILENAME = f'{ATOM_NAME}_MSD_all_product.dat'

INPUT_FILE = os.path.join(INPUT_FILE_PATH, INPUT_FILENAME)

try:
    df = pd.read_csv(INPUT_FILE,
                     sep='\s+',
                     comment='#',
                     names=['Time', 'MSD_avg', 'MSD_std'])

    time_lags = df['Time'].values
    msd_mean = df['MSD_avg'].values
    print(f"Successfully loaded data from '{INPUT_FILE}'.")

except FileNotFoundError:
    print(f"ERROR: File not found at '{INPUT_FILE}'. Please check the file name and location.")
    exit()
except Exception as e:
    print(f"ERROR: An error occurred while reading the file: {e}")
    exit()

valid_indices = np.where(time_lags > 0)
log_time = np.log10(time_lags[valid_indices])
log_msd = np.log10(msd_mean[valid_indices])

alpha = np.diff(log_msd) / np.diff(log_time)

time_for_alpha = time_lags[valid_indices][1:]

plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(10, 7))

ax.plot(time_for_alpha, alpha, color='green', linestyle='-', linewidth=2.5, label=r'Anomalous Exponent $\alpha(t)$')

ax.axhline(1, color='gray', linestyle='--', label=r'$\alpha=1$ (Normal Diffusion)')

ax.set_xlabel('Time (ps)', fontsize=14)
ax.set_ylabel(r'$\alpha$', fontsize=14)
ax.set_title(rf'"{ATOM_NAME}" Time vs $\alpha$', fontsize=16, fontweight='bold')
ax.set_xlim(0.05, 50)
ax.set_ylim(0, 1)
ax.set_xscale('log')
ax.legend(fontsize=12)
ax.tick_params(axis='both', which='major', labelsize=12)

fig.tight_layout()

output_fig_filename = f'{ATOM_NAME}_MSD_alpha.png'

output_path = "../../result" 
output_fig_file = os.path.join(output_path, output_fig_filename)
plt.savefig(output_fig_file, dpi=300)

print(f"\nPlot successfully saved as '{output_fig_file}'")

try:
    df_alpha = pd.DataFrame({
        'Time(ps)': time_for_alpha,
        'Alpha': alpha
    })

    output_data_filename = f'{ATOM_NAME}_MSD_alpha_values.dat'
    output_data_file = os.path.join(output_path, output_data_filename)

    df_alpha.to_csv(output_data_file, 
                    sep='\t', 
                    index=False, 
                    float_format='%.8e',
                    header=['# Time(ps)', 'Alpha'])
    
    print(f"Alpha data successfully saved as '{output_data_file}'")

except Exception as e:
    print(f"ERROR: An error occurred while saving the alpha data: {e}")

plt.show()