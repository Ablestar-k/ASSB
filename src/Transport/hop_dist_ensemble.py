import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, linregress
import os

NUM_ENSEMBLES = 5


def load_hop_data(filepath: str) -> np.ndarray:
    if not os.path.exists(filepath):
        print(f"ERROR: File not found at '{filepath}'")
        return None
    
    all_hop_values = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.strip().startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) > 1:
                hop_values = [float(val) for val in parts[1:]]
                all_hop_values.extend(hop_values)
                
    return np.array(all_hop_values)


def analyze_ensemble_statistics(all_data: list, base_path: str):
    print(f"\n{'='*15} Analyzing Ensemble Statistics & Fitting {'='*15}")
    
    if not all_data:
        print("No data provided for ensemble analysis.")
        return

    min_h = 0.0
    max_h = 5.0
    h_grid = np.linspace(min_h, max_h, 4000)

    all_ccdfs = np.array([[np.sum(data > h) / len(data) for h in h_grid] for data in all_data])
    
    mean_ccdf = np.mean(all_ccdfs, axis=0)
    std_ccdf = np.std(all_ccdfs, axis=0)


    pdf_approx = -np.gradient(mean_ccdf, h_grid)
    valid_range_mask = (mean_ccdf < 0.95) & (mean_ccdf > 1e-8)
    if np.any(valid_range_mask):
        peak_index = np.argmax(pdf_approx[valid_range_mask])
        h_at_peak = h_grid[valid_range_mask][peak_index]
        print(f"[Analysis] Estimated h at PDF Peak: {h_at_peak:.4f} Å^2")
    else:
        print("[Warning] Could not determine a valid PDF peak. Using a fallback.")
        h_at_peak = h_grid[len(h_grid) // 4]

    tail_mask = (h_grid > h_at_peak + 1.5) & (mean_ccdf > 1e-8) 
    h_tail = h_grid[tail_mask]
    ccdf_tail = mean_ccdf[tail_mask]
    
    h_departure = None
    fit_line_log_values = None 

    if len(h_tail) > 10:
        ccdf_tail_log = np.log(ccdf_tail)
        
        slope, intercept, r_value, _, _ = linregress(h_tail, ccdf_tail_log)
        fit_line_log_values = slope * h_grid + intercept
        print(f"[Analysis] Exponential fit completed on the tail (R^2 = {r_value**2:.4f}).")

        residuals = abs(ccdf_tail_log - (slope * h_tail + intercept)) 
        residual_std = np.std(residuals)
        threshold = 0.1 * residual_std  

        departure_indices = np.where(residuals > threshold)[0]
        if len(departure_indices) > 0:
            first_departure_index = departure_indices[0]
            h_departure = h_tail[first_departure_index]
            print(f"[Result] Departure from exponential fit found at h = {h_departure:.4f} Å^2")
        else:
            print("[Info] No significant departure point found in the tail.")
            
    else:
        print("[Info] Insufficient tail data to perform log-fitting.")

    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(10, 7))

    ax.semilogy(h_grid, mean_ccdf, 'b-', lw=2, label='Mean CCDF')
    ax.fill_between(h_grid, mean_ccdf - std_ccdf, mean_ccdf + std_ccdf, 
                     color='b', alpha=0.2, label='Standard Deviation')

    if fit_line_log_values is not None:
        ax.semilogy(h_grid, np.exp(fit_line_log_values), 'r--', lw=2, label=f'Exponential Fit (Extended)')
    
    ax.set_title('Ensemble Average CCDF with Exponential Fit', fontsize=16)
    ax.set_xlabel('Hop value, $h$ (Å^2)', fontsize=12)
    ax.set_ylabel('Probability, $P(H > h)$', fontsize=12)
    ax.set_xlim(0, 5)
    ax.set_ylim(bottom=1e-5, top=1.1)
    ax.legend()
    
    plt.tight_layout()
    save_path = f'{base_path}ensemble_fit_analysis_extended.png'
    plt.show()
    plt.savefig(save_path, dpi=300)
    plt.close(fig)
    print(f"Ensemble fit analysis plot saved to '{save_path}'")

if __name__ == '__main__':
    base_path = '../../result/hop/'
    os.makedirs(base_path, exist_ok=True)

    all_hop_data = []
    for i in range(1, NUM_ENSEMBLES + 1):
        print(f"\n--- Loading file number {i} ---")
        input_file = f'{base_path}{i}_hop.dat'
        
        hop_data = load_hop_data(input_file)
        
        if hop_data is not None:
            all_hop_data.append(hop_data)

    if all_hop_data:
        analyze_ensemble_statistics(all_hop_data, base_path)
    
    print("\nAll processing finished.")