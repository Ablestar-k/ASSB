import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

ATOM_NAME = 'Na'
INPUT_FILE_PATH = "../result"
INPUT_FILENAME = f'{ATOM_NAME}_MSD_ensemble_average.dat'

INPUT_FILE = os.path.join(INPUT_FILE_PATH, INPUT_FILENAME)

WINDOW_SIZE = 10  

try:
    df = pd.read_csv(INPUT_FILE,
                     sep='\s+',
                     comment='#',
                     names=['Time', 'MSD_avg', 'MSD_std'])
    time_lags = df['Time'].values
    msd_mean = df['MSD_avg'].values
    print(f"Successfully loaded data from '{INPUT_FILE}'.")
except FileNotFoundError:
    print(f"ERROR: File not found at '{INPUT_FILE}'.")
    exit()
except Exception as e:
    print(f"ERROR: An error occurred while reading the file: {e}")
    exit()

valid_indices = np.where(time_lags > 0)
log_time = np.log10(time_lags[valid_indices])
log_msd = np.log10(msd_mean[valid_indices])


alpha_noisy = np.diff(log_msd) / np.diff(log_time)
time_noisy = time_lags[valid_indices][1:]


alpha_smooth = []
time_centers = []


for i in range(len(log_time) - WINDOW_SIZE + 1):
    t_window = log_time[i : i + WINDOW_SIZE]
    msd_window = log_msd[i : i + WINDOW_SIZE]
    
    slope, intercept = np.polyfit(t_window, msd_window, 1)
    alpha_smooth.append(slope)
    
    center_index = i + WINDOW_SIZE // 2
    time_centers.append(time_lags[valid_indices][center_index])

plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(12, 8))

ax.plot(time_noisy, alpha_noisy, color='mediumseagreen', linestyle='-', linewidth=1.5, alpha=0.6, label=r'Local Slope (Noisy)')

ax.plot(time_centers, alpha_smooth, color='royalblue', linestyle='-', linewidth=3.0, label=f'Sliding Window (size={WINDOW_SIZE}, Smoothed)')

ax.axhline(1, color='gray', linestyle='--', linewidth=1.5, label=r'$\alpha=1$ (Normal Diffusion)')

ax.set_xlabel('Time (ps)', fontsize=14)
ax.set_ylabel(r'Anomalous Exponent $\alpha$', fontsize=14)
ax.set_title(r'Comparison of Methods for Calculating $\alpha(t)$', fontsize=16, fontweight='bold')
ax.set_ylim(-0.5, 2.5)
ax.set_xscale('log') 
ax.legend(fontsize=12)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.grid(True, which="both", linestyle='--', alpha=0.6)

fig.tight_layout()

output_fig_filename = f'{ATOM_NAME}_MSD_alpha_window.png'
output_fig_path = "../result"
output_fig_file = os.path.join(output_fig_path, output_fig_filename)
plt.savefig(output_fig_file, dpi=300, bbox_inches='tight')

plt.savefig('alpha_comparison.png', dpi=300)
print("\nPlot successfully saved as 'alpha_comparison.png'")
plt.show()