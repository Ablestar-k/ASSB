import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

INPUT_FILENAME = 'MSD_ensemble_average.dat'
WINDOW_SIZE = 10

try:
    df = pd.read_csv(INPUT_FILENAME, sep='\s+', comment='#', names=['Time', 'MSD_avg', 'MSD_std'])
    time_lags = df['Time'].values
    msd_mean = df['MSD_avg'].values
    print(f"Successfully loaded data from '{INPUT_FILENAME}'.")
except FileNotFoundError:
    print(f"ERROR: File not found at '{INPUT_FILENAME}'.")
    exit()

valid_indices = np.where(time_lags > 0)
log_time = np.log10(time_lags[valid_indices])
log_msd = np.log10(msd_mean[valid_indices])

alpha_smooth = []
time_centers = []
for i in range(len(log_time) - WINDOW_SIZE + 1):
    t_window = log_time[i : i + WINDOW_SIZE]
    msd_window = log_msd[i : i + WINDOW_SIZE]
    slope, _ = np.polyfit(t_window, msd_window, 1)
    alpha_smooth.append(slope)
    center_index = i + WINDOW_SIZE // 2
    time_centers.append(time_lags[valid_indices][center_index])


time_centers = np.array(time_centers)
alpha_smooth = np.array(alpha_smooth)

REGIME_TIME_RANGES = {
    'Ballistic': (0.01, 0.1),
    'Sub-diffusive': (1, 100),
    'Fickian': (1000, 10000)
}


print("\n--- Quantitative Analysis of Diffusion Regimes ---")
regime_results = {}
for regime_name, (start_time, end_time) in REGIME_TIME_RANGES.items():
    mask = (time_centers >= start_time) & (time_centers <= end_time)
    

    if np.any(mask):
        alpha_in_regime = alpha_smooth[mask]
        avg_alpha = np.mean(alpha_in_regime)
        std_alpha = np.std(alpha_in_regime)
        regime_results[regime_name] = {'avg': avg_alpha, 'std': std_alpha, 'range': (start_time, end_time)}
        print(f"'{regime_name}' regime ({start_time}-{end_time} ps):")
        print(f"  Average alpha = {avg_alpha:.3f} ± {std_alpha:.3f}")
    else:
        print(f"'{regime_name}' regime ({start_time}-{end_time} ps): No data points found in this range.")


plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(12, 8))


ax.plot(time_centers, alpha_smooth, color='royalblue', linestyle='-', linewidth=2.5, label=r'Anomalous Exponent $\alpha(t)$')


colors = {'Ballistic': 'firebrick', 'Sub-diffusive': 'darkorange', 'Fickian': 'forestgreen'}
for name, results in regime_results.items():
    start, end = results['range']
    color = colors.get(name, 'gray') 
    ax.axvspan(start, end, color=color, alpha=0.15, label=f'{name} Regime')
    
    text_x_pos = np.sqrt(start * end)
    text_y_pos = results['avg']
    ax.text(text_x_pos, text_y_pos + 0.1, f'$\\alpha \\approx {text_y_pos:.2f}$',
            horizontalalignment='center', fontsize=12, fontweight='bold', color=color)


ax.axhline(1, color='gray', linestyle='--', linewidth=1.5, label=r'$\alpha=1$ (Normal Diffusion)')

ax.set_xlabel('Time (ps)', fontsize=14)
ax.set_ylabel(r'Anomalous Exponent $\alpha$', fontsize=14)
ax.set_title(r'Diffusion Regimes Identified from $\alpha(t)$', fontsize=16, fontweight='bold')
ax.set_ylim(-0.5, 2.5)
ax.set_xscale('log')
ax.legend(fontsize=12, loc='upper right')
ax.tick_params(axis='both', which='major', labelsize=12)
ax.grid(True, which="both", linestyle='--', alpha=0.6)

fig.tight_layout()
plt.savefig('alpha_regime_analysis.png', dpi=300)
print("\nPlot with identified regimes saved as 'alpha_regime_analysis.png'")
plt.show()