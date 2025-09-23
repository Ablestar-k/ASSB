import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
import os

ATOM_NAME = 'Na'
INPUT_FILE_PATH = '../result/MSD'
INPUT_FILENAME = f'{ATOM_NAME}_MSD_ensemble_average.dat'

INPUT_FILE = os.path.join(INPUT_FILE_PATH, INPUT_FILENAME)

FIT_START_TIME_PS = 600.0
FIT_END_TIME_PS = 800.0

DIMENSIONALITY = 3

try:

    df = pd.read_csv(INPUT_FILE,
                     sep='\s+',        
                     comment='#',        
                     names=['Time', 'MSD_avg', 'MSD_std'])

    time_lags = df['Time'].values
    msd_mean = df['MSD_avg'].values
    msd_std = df['MSD_std'].values
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


fit_indices = np.where((time_lags >= FIT_START_TIME_PS) & (time_lags <= FIT_END_TIME_PS))

if len(fit_indices[0]) < 2:
    print("Warning: Not enough data points in the selected range for linear fitting.")
    D_cm2_per_s, slope, intercept = 0, 0, 0
else:
    time_fit = time_lags[fit_indices]
    msd_fit = msd_mean[fit_indices]
    
    lin_model = linregress(time_fit, msd_fit)
    slope = lin_model.slope
    intercept = lin_model.intercept
    
    D_A2_per_ps = slope / (2 * DIMENSIONALITY)
    D_cm2_per_s = D_A2_per_ps * 1e-4
    
    print("\n--- Diffusion Coefficient Analysis ---")
    print(f"Linear fit range: {FIT_START_TIME_PS} ps to {FIT_END_TIME_PS} ps")
    print(f"Slope (2dD): {slope:.4f} Å^2/ps")
    print(f"R-squared value: {lin_model.rvalue**2:.4f}")
    print(f"Diffusion Coefficient (D): {D_cm2_per_s:.4e} cm^2/s")
    print("------------------------------------")


fig, ax1 = plt.subplots(figsize=(10, 7))
plt.style.use('seaborn-v0_8-whitegrid')

ax1.plot(time_lags, msd_mean, label='MSD', color='royalblue', linewidth=2.5)
ax1.fill_between(time_lags, msd_mean - msd_std, msd_mean + msd_std, color='cornflowerblue', alpha=0.3, label='Std Dev')
ax1.set_xlabel('Time (ps)', fontsize=14)
ax1.set_ylabel('MSD ($Å^2$)', fontsize=14, color='royalblue')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.tick_params(axis='y', labelcolor='royalblue')

if len(fit_indices[0]) >= 2:
    fit_line = slope * time_fit + intercept
    ax1.plot(time_fit, fit_line, 'r--', linewidth=2, label=f'Linear Fit (D = {D_cm2_per_s:.2e} cm²/s)')

ax1.legend(loc='upper left')

#ax2 = ax1.twinx()
#ax2.plot(time_lags[valid_indices][1:], alpha, color='green', linestyle=':', linewidth=2, label=r'Anomalous Exponent $\alpha(t)$')
#ax2.set_ylabel(r'Slope $\alpha$', fontsize=14, color='green')
#ax2.tick_params(axis='y', labelcolor='green')
#ax2.set_ylim(0, 2)
#ax2.axhline(1, color='gray', linestyle='--', label=r'$\alpha=1$ (Normal Diffusion)')
#ax2.legend(loc='lower right')
#
#plt.title('MSD and Anomalous Diffusion Analysis', fontsize=16, fontweight='bold')
#fig.tight_layout()

output_path = "../result"
output_fig_filename = f'MSD_{ATOM_NAME}_diffusion_analysis.png'
output_fig_file = os.path.join(output_path, output_fig_filename)

plt.savefig(output_fig_file, dpi=300)
print(f"\nPlot successfully saved as '{output_fig_file}'")
plt.show()