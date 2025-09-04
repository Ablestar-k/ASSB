import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import statsmodels.api as sm

file_path = "../result/quench_pot_temp.txt"
output_path = "../result"

os.makedirs(output_path, exist_ok=True)

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

    df = df[df['PE_std'] > 0]
    print(f"Successfully loaded and processed data from '{file_path}'.")
except FileNotFoundError:
    print(f"Error: The file '{file_path}' was not found. Please check the file path.")
    exit()


fit_range_1 = (df['Temp'] >= 1000) & (df['Temp'] <= 1200)
df1 = df[fit_range_1].copy()


fit_range_2 = (df['Temp'] >= 300) & (df['Temp'] <= 500)
df2 = df[fit_range_2].copy()



weights1 = 1 / (df1['PE_std']**2)
X1 = sm.add_constant(df1['Temp'])
wls_model1 = sm.WLS(df1['PE'], X1, weights=weights1)
results1 = wls_model1.fit()
intercept1, slope1 = results1.params

weights2 = 1 / (df2['PE_std']**2)
X2 = sm.add_constant(df2['Temp']) 
wls_model2 = sm.WLS(df2['PE'], X2, weights=weights2)
results2 = wls_model2.fit()
intercept2, slope2 = results2.params

print("\n--- Weighted Linear Regression Results ---")
print(f"Range 1 (1000K - 1200K): Slope={slope1:.4f}, Intercept={intercept1:.4f}")
print(f"Range 2 (300K - 500K):  Slope={slope2:.4f}, Intercept={intercept2:.4f}")


if (slope1 - slope2) != 0:
    intersect_temp = (intercept2 - intercept1) / (slope1 - slope2)
    intersect_pe = slope1 * intersect_temp + intercept1
    print(f"\nIntersection Point (Glass Transition Temperature, Tg):")
    print(f"  - Temperature: {intersect_temp:.2f} K")
    print(f"  - Potential Energy: {intersect_pe:.2f} eV")
else:
    print("\nSlopes are parallel, cannot calculate intersection.")
    intersect_temp, intersect_pe = None, None

# --- Plotting ---
plt.style.use('seaborn-v0_8-whitegrid')
plt.figure(figsize=(12, 8))


plt.plot(df['Temp'], df['PE'], 'o', markersize=3, color='gray', label='All Data')


plt.fill_between(
    df['Temp'],
    df['PE'] - df['PE_std'],
    df['PE'] + df['PE_std'],
    color='gray',
    alpha=0.2,
    label='Standard Deviation'
)

fit_x1 = np.array([df1['Temp'].min(), df['Temp'].max()])
fit_y1 = slope1 * fit_x1 + intercept1
plt.plot(fit_x1, fit_y1, color='blue', linestyle='--', linewidth=2, label=f'Fitting 1')


fit_x2 = np.array([df['Temp'].min(), df2['Temp'].max()])
fit_y2 = slope2 * fit_x2 + intercept2
plt.plot(fit_x2, fit_y2, color='green', linestyle='--', linewidth=2, label=f'Fitting 2')

if intersect_temp is not None:
    plt.plot(intersect_temp, intersect_pe, 'r*', markersize=15, label=f'Intersection (Tg): {intersect_temp:.2f} K')
    plt.annotate(
        f'Tg = {intersect_temp:.2f} K\nPE = {intersect_pe:.2f} eV',
        xy=(intersect_temp, intersect_pe),
        xytext=(intersect_temp + 50, intersect_pe - 0.05),
        arrowprops=dict(facecolor='red', shrink=0.05, width=1, headwidth=8),
        fontsize=12,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=1, alpha=0.7)
    )

plt.title('Temperature vs. Potential Energy Linear Fits', fontsize=18)
plt.xlabel('Temperature (K)', fontsize=14)
plt.ylabel('Potential Energy (eV)', fontsize=14)
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()

save_path = os.path.join(output_path, 'T_vs_PE_linear_fit.png')
plt.savefig(save_path, dpi=300)
plt.close()

print(f"\nGraph successfully saved to '{save_path}'.")