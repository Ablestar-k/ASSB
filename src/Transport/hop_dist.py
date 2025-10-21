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

def analyze_and_save_distributions(data: np.ndarray, file_id: int, base_path: str):
    if data is None or len(data) < 20: 
        print("Not enough data for analysis.")
        return

    # --- 1. Calculate Distributions (PDF, CDF, CCDF) ---
    # PDF using Kernel Density Estimation
    kde = gaussian_kde(data)
    x_range_pdf = np.linspace(data.min(), data.max(), 1000)
    pdf_values = kde(x_range_pdf)
    
    # Empirical CDF and CCDF
    sorted_data = np.sort(data)
    y_values_cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
    y_values_ccdf = 1 - y_values_cdf

    # --- 2. Save Distribution Data ---
    # Define output paths
    pdf_path = f'{base_path}{file_id}_hop_pdf.dat'
    cdf_path = f'{base_path}{file_id}_hop_cdf.dat'
    ccdf_path = f'{base_path}{file_id}_hop_ccdf.dat'
    
    np.savetxt(pdf_path, np.vstack((x_range_pdf, pdf_values)).T, header='Hop_Value_(h)\tPDF', fmt='%.8f', delimiter='\t')
    np.savetxt(cdf_path, np.vstack((sorted_data, y_values_cdf)).T, header='Hop_Value_(h)\tCDF', fmt='%.8f', delimiter='\t')
    np.savetxt(ccdf_path, np.vstack((sorted_data, y_values_ccdf)).T, header='Hop_Value_(h)\tCCDF_(1-CDF)', fmt='%.8f', delimiter='\t')
    print(f"Distribution data saved for dataset {file_id}.")

    # --- 3. Find PDF Peak (Steepest drop in CCDF) ---
    peak_index = np.argmax(pdf_values)
    h_at_peak = x_range_pdf[peak_index]
    print(f"\n[Analysis Result 1] h at PDF Peak: {h_at_peak:.4f} Å^2")
    
    # --- 4. Find Departure Point from Log-Linear Fit ---
    tail_mask = sorted_data > h_at_peak + 0.25
    h_tail = sorted_data[tail_mask]
    ccdf_tail = y_values_ccdf[tail_mask]

    valid_mask = ccdf_tail > 0
    h_tail_valid = h_tail[valid_mask]
    ccdf_tail_log = np.log(ccdf_tail[valid_mask])

    h_departure = None
    fit_line_log = None
    if len(h_tail_valid) > 10: 
        slope, intercept, _, _, _ = linregress(h_tail_valid, ccdf_tail_log)
        fit_line_log = slope * h_tail_valid + intercept
        
        residuals = abs(ccdf_tail_log - fit_line_log)
 
        residual_std = np.std(residuals)
        threshold = 2.0 * residual_std

        # Find the first index where the residual exceeds the threshold
        departure_indices = np.where(residuals > threshold)[0]
        if len(departure_indices) > 0:
            first_departure_index = departure_indices[0]
            h_departure = h_tail_valid[first_departure_index]
            print(f"[Analysis Result 2] h at Log-Fit Departure: {h_departure:.4f} Å^2 (Threshold: {threshold:.4f})")
    else:
        print("[Info] Insufficient tail data to perform log-fitting.")

    # --- 5. Visualization ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle(f'Hop Value Distribution Analysis (Dataset {file_id})', fontsize=16)

    # Plot 1: PDF with key points
    ax1.plot(x_range_pdf, pdf_values, 'r-', lw=2, label='KDE-PDF')
    ax1.hist(data, bins=100, density=True, alpha=0.5, label='Histogram')
    ax1.axvline(h_at_peak, color='k', linestyle='--', label=f'PDF Peak h = {h_at_peak:.2f} Å^2')
    if h_departure is not None:
        ax1.axvline(h_departure, color='b', linestyle=':', lw=2, label=f'Departure h = {h_departure:.2f} Å^2')
    ax1.set_title('Probability Density Function (PDF)', fontsize=14)
    ax1.set_xlim(0,6)
    ax1.set_xlabel('Hop value, $h$ (Å^2)', fontsize=12)
    ax1.set_ylabel('Probability Density', fontsize=12)
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.6)

    # Plot 2: CCDF on a log scale with fit
    ax2.semilogy(sorted_data, y_values_ccdf, '.', ms=2, color='green', label='CCDF Data')
    if h_departure is not None and fit_line_log is not None:
        ax2.semilogy(h_tail_valid, np.exp(fit_line_log), 'r-', label='Exponential Fit')
        ax2.axvline(h_departure, color='b', linestyle=':', lw=2, label=f'Departure h = {h_departure:.2f} Å^2')
    ax2.set_title('CCDF (Log Scale) and Exponential Fit', fontsize=14)
    ax2.set_xlabel('Hop value, $h$ (Å^2)', fontsize=12)
    ax2.set_ylabel('Probability, $P(H > h)$', fontsize=12)
    ax2.set_yscale('log')
    ax2.legend()
    ax2.grid(True, which="both", linestyle='--', alpha=0.6)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(f'{base_path}{file_id}_hop_logfit_analysis.png', dpi=300)
    #plt.show()

if __name__ == '__main__':
    base_path = '../../result/hop/'

    os.makedirs(base_path, exist_ok=True)

    for i in range(1, NUM_ENSEMBLES + 1):
        print(f"\n{'='*15} Processing file number {i} {'='*15}")
        input_file = f'{base_path}{i}_hop.dat'
        
        hop_data = load_hop_data(input_file)
        
        if hop_data is not None:
            analyze_and_save_distributions(
                data=hop_data,
                file_id=i,
                base_path=base_path
            )
        print(f"{'='*15} Finished processing file number {i} {'='*15}")