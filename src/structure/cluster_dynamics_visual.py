import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


BASE_PATH = '../../result/cluster_dynamics/'

DELTA_T_STEP_C_CODE = 1

TRAJ_INTERVAL_PS = 2.0

ENSEMBLE_START = 1
ENSEMBLE_END = 5

LIFETIME_HIST_BINS = 10000

# These must match the prefixes of the .dat files (e.g., "1_tcf_poly.dat")
PREFIX_TCF_POLY = "tcf_intermit_poly"
PREFIX_TCF_ISO = "tcf_intermit_iso"
PREFIX_TCF_SIZE = "tcf_size_corr"
PREFIX_LIFE_ISO_POLY = "life_iso_poly"
PREFIX_LIFE_POLY_ISO = "life_poly_iso"


PLOT_TCF_POLY = "ensemble_avg_tcf_poly.png"
PLOT_TCF_ISO = "ensemble_avg_tcf_iso.png"
PLOT_TCF_COMBINED = "ensemble_avg_tcf_combined.png"
PLOT_TCF_SIZE = "ensemble_avg_tcf_size_corr.png"
PLOT_LIFE_ISO_POLY = "ensemble_sum_life_iso_poly.png"
PLOT_LIFE_POLY_ISO = "ensemble_sum_life_poly_iso.png"


# --- Internal Setup ---
ENSEMBLE_RANGE = range(ENSEMBLE_START, ENSEMBLE_END + 1)
NUM_ENSEMBLES = len(ENSEMBLE_RANGE)

if not os.path.exists(BASE_PATH):
    print(f"Warning: Directory {BASE_PATH} does not exist. The script might fail.")


# --- 2. Ensemble Averaging Function (for TCFs) ---

def process_and_average_time_data(ensemble_range, base_path, file_prefix, num_ensembles):
    """
    Reads TCF data (.dat) from all ensembles, calculates the mean and standard error (SEM).
    (Used for tcf_poly, tcf_iso, tcf_size_corr)
    """
    print(f"--- Starting Ensemble TCF Processing ({file_prefix}) ---")
    raw_data_list = []
    
    for i in ensemble_range:
        filename = os.path.join(base_path, f'{i}_{file_prefix}.dat')
        if not os.path.exists(filename):
            print(f"  File not found (Skipping): {filename}")
            continue
        try:
            # C code output format: dt(step), TCF, <h(0)h(t)>, <h(0)>
            data = pd.read_csv(filename, sep='\s+', comment='#', 
                               names=['dt_frame', 'TCF', 'h0ht', 'h0'])
            if data.empty:
                print(f"  File is empty: {filename}")
                continue
            
            raw_data_list.append(data)
            print(f"  Loaded successfully: {filename} ({len(data)} data points)")
        except Exception as e:
            print(f"  Error reading file {filename}: {e}")

    if len(raw_data_list) == 0:
        print(f"--- Ensemble TCF Processing FAILED: No data loaded ({file_prefix}) ---")
        return None, None
    elif len(raw_data_list) != num_ensembles:
        print(f"  Warning: Requested ensemble count ({num_ensembles}) differs from loaded count ({len(raw_data_list)}).")

    # --- Find common time axis length ---
    base_t_frames = raw_data_list[0]['dt_frame'].values
    min_common_length = len(base_t_frames)
    
    for data in raw_data_list[1:]:
        current_t_frames = data['dt_frame'].values
        new_len = 0
        for k in range(min(min_common_length, len(current_t_frames))):
            if base_t_frames[k] == current_t_frames[k]:
                new_len = k + 1
            else:
                break
        if new_len < min_common_length:
            min_common_length = new_len
            
    print(f"  Averaging all ensembles up to common length {min_common_length} (dt = {base_t_frames[min_common_length-1]} frames).")
    
    # Truncate to common length
    t_values_frames = base_t_frames[:min_common_length]
    truncated_df_list = [data.iloc[:min_common_length] for data in raw_data_list]

    # Stack the 'TCF' column (column 1) for ensemble averaging
    tcf_stack = np.stack([df['TCF'].values for df in truncated_df_list], axis=0)
    
    # Calculate ensemble average and standard error (SEM)
    avg_tcf = np.mean(tcf_stack, axis=0)
    sem_tcf = np.std(tcf_stack, axis=0) / np.sqrt(tcf_stack.shape[0])
    
    print(f"  Ensemble averaging complete. (N={tcf_stack.shape[0]})")

    # --- Save average/SEM data ---
    t_values_ps = t_values_frames * TRAJ_INTERVAL_PS
    
    save_df = pd.DataFrame({
        'dt_frame': t_values_frames,
        'dt_ps': t_values_ps,
        'TCF_Avg': avg_tcf,
        'TCF_SEM': sem_tcf
    })
    
    avg_output_file = os.path.join(base_path, f'ensemble_avg_{file_prefix}.dat')
    
    header = f"# Ensemble Average TCF for {file_prefix}\n"
    header += f"# Averaged over {tcf_stack.shape[0]} ensembles.\n"
    header += f"# Truncated to common length {min_common_length}.\n"
    header += f"# dt_frame dt_ps TCF_Avg TCF_SEM"
    
    save_df.to_csv(avg_output_file, sep=' ', index=False, float_format='%.8e', header=header)

    print(f"  Successfully saved ensemble average: {avg_output_file}")
    print(f"--- Ensemble TCF Processing FINISHED ({file_prefix}) ---\n")
    
    return avg_output_file, save_df

# --- 3. Ensemble Averaging Function (for Lifetimes) ---

def process_and_average_lifetime_dist(ensemble_range, base_path, file_prefix, num_ensembles, hist_bins):
    """
    Reads Lifetime distribution (.dat) from all ensembles, sums them into one histogram, and normalizes.
    (Used for life_iso_poly, life_poly_iso)
    """
    print(f"--- Starting Ensemble Lifetime Processing ({file_prefix}) ---")
    
    total_counts_hist = np.zeros(hist_bins, dtype=np.int64)
    total_events = 0
    files_found = 0

    for i in ensemble_range:
        filename = os.path.join(base_path, f'{i}_{file_prefix}.dat')
        if not os.path.exists(filename):
            print(f"  File not found (Skipping): {filename}")
            continue
        try:
            data = pd.read_csv(filename, sep='\s+', comment='#', 
                               names=['lifetime_frame', 'Count'])
            if data.empty:
                print(f"  File is empty: {filename}")
                continue
            
            files_found += 1
            for _, row in data.iterrows():
                lifetime = int(row['lifetime_frame'])
                count = int(row['Count'])
                if lifetime < hist_bins:
                    total_counts_hist[lifetime] += count
                    total_events += count
            
            print(f"  Loaded successfully: {filename} (Added {data['Count'].sum()} events)")
        except Exception as e:
            print(f"  Error reading file {filename}: {e}")

    if total_events == 0:
        print(f"--- Ensemble Lifetime Processing FAILED: No valid events found ({file_prefix}) ---")
        return None, None

    print(f"  Found {total_events} total events across {files_found} ensembles.")

    # --- Save summed/normalized data ---
    lifetime_frames = np.arange(hist_bins)
    lifetime_ps = lifetime_frames * TRAJ_INTERVAL_PS
    
    if total_events > 0:
        probability = total_counts_hist / total_events
    else:
        probability = np.zeros(hist_bins)

    save_df = pd.DataFrame({
        'lifetime_frame': lifetime_frames,
        'lifetime_ps': lifetime_ps,
        'Count_Sum': total_counts_hist,
        'Probability': probability
    })

    avg_output_file = os.path.join(base_path, f'ensemble_sum_{file_prefix}.dat')
    
    header = f"# Ensemble Summed Lifetime Distribution for {file_prefix}\n"
    header += f"# Summed over {files_found} ensembles (Total Events: {total_events}).\n"
    header += f"# lifetime_frame lifetime_ps Count_Sum Probability"
    
    save_df_sparse = save_df[save_df['Count_Sum'] > 0]
    save_df_sparse.to_csv(avg_output_file, sep=' ', index=False, float_format='%.8e', header=header)

    print(f"  Successfully saved ensemble sum distribution: {avg_output_file}")
    print(f"--- Ensemble Lifetime Processing FINISHED ({file_prefix}) ---\n")
    
    return avg_output_file, save_df_sparse

# --- 4. Plotting Function (for TCFs) ---

def plot_tcf_data(avg_data_df, title, output_filename, y_label="TCF (Ensemble Avg ± SEM)"):
    """Saves a single TCF dataset as a (log-y) plot."""
    print(f"  Plotting started: {output_filename}")
    plt.figure(figsize=(12, 8))
    
    t_ps = avg_data_df['dt_ps'].values
    tcf_avg = avg_data_df['TCF_Avg'].values
    tcf_sem = avg_data_df['TCF_SEM'].values

    plt.plot(t_ps, tcf_avg, label=f'Ensemble Avg (N={NUM_ENSEMBLES})', color='blue', linewidth=2)
    plt.fill_between(t_ps, 
                     tcf_avg - tcf_sem, 
                     tcf_avg + tcf_sem, 
                     color='blue', alpha=0.2, label='SEM')

    plt.title(title)
    plt.xlabel('Time dt (ps)')
    plt.ylabel(y_label)
    plt.xscale('log')
    #plt.yscale('log')
    plt.legend()
    plt.ylim(bottom=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(BASE_PATH, output_filename))
    plt.close()
    print(f"  Plot saved successfully: {output_filename}\n")

# --- 5. Plotting Function (for Lifetimes) ---

def plot_lifetime_data(avg_data_df, title, output_filename):
    """Saves Lifetime distribution data as a (log-y) plot."""
    print(f"  Plotting started: {output_filename}")
    plt.figure(figsize=(12, 8))
    
    t_ps = avg_data_df['lifetime_ps'].values
    prob = avg_data_df['Probability'].values
    
    t_ps_plot = t_ps[prob > 0]
    prob_plot = prob[prob > 0]

    if len(prob_plot) == 0:
        print("  Warning: No data to plot for Lifetime.")
        plt.close()
        return

    plt.plot(t_ps_plot, prob_plot, 'o-', markersize=3, label=f'Prob. Dist. (N={NUM_ENSEMBLES} total)')
    
    plt.title(title)
    plt.xlabel('Lifetime $\\tau$ (ps)')
    plt.ylabel('Probability P($\\tau$) (Log Scale)')
    plt.yscale('log')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(BASE_PATH, output_filename))
    plt.close()
    print(f"  Plot saved successfully: {output_filename}\n")

# --- 6. Plotting Function (Combined TCFs) ---

def plot_combined_tcf(poly_df, iso_df, title, output_filename):
    """Plots Polymeric TCF and Isolated TCF on the same graph."""
    print(f"  Combined plotting started: {output_filename}")
    plt.figure(figsize=(12, 8))

    plt.plot(poly_df['dt_ps'], poly_df['TCF_Avg'], 
             label=f'Poly-Ta Pair (N={NUM_ENSEMBLES})', color='blue', linewidth=2)
    plt.fill_between(poly_df['dt_ps'], 
                     poly_df['TCF_Avg'] - poly_df['TCF_SEM'], 
                     poly_df['TCF_Avg'] + poly_df['TCF_SEM'], 
                     color='blue', alpha=0.15)

    plt.plot(iso_df['dt_ps'], iso_df['TCF_Avg'], 
             label=f'Iso-Ta State (N={NUM_ENSEMBLES})', color='red', linewidth=2)
    plt.fill_between(iso_df['dt_ps'], 
                     iso_df['TCF_Avg'] - iso_df['TCF_SEM'], 
                     iso_df['TCF_Avg'] + iso_df['TCF_SEM'], 
                     color='red', alpha=0.15)

    plt.title(title)
    plt.xlabel('Time dt (ps)')
    plt.ylabel('TCF (Ensemble Avg ± SEM)')
    plt.ylim(bottom=0.7)
    plt.xscale('log')
    #plt.yscale('log')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(BASE_PATH, output_filename))
    plt.close()
    print(f"  Combined plot saved successfully: {output_filename}\n")


# --- 7. Main Script Execution ---
if __name__ == "__main__":
    print("--- Ensemble Dynamics Analysis Script ---")
    print(f"Ensemble range: {ENSEMBLE_START} to {ENSEMBLE_END} (Total {NUM_ENSEMBLES})")
    print(f"Data Path: {BASE_PATH}")
    print(f"Trajectory Interval: {TRAJ_INTERVAL_PS} ps/frame\n")
    
    poly_tcf_df = None
    iso_tcf_df = None

    # --- 1. Intermittent TCF (Poly) ---
    poly_avg_file, poly_tcf_df = process_and_average_time_data(
        ENSEMBLE_RANGE, BASE_PATH, PREFIX_TCF_POLY, NUM_ENSEMBLES
    )
    if poly_avg_file:
        plot_tcf_data(poly_tcf_df, 
                      "Intermittent TCF (Polymeric Ta-Ta Pairs)",
                      PLOT_TCF_POLY,
                      "C_poly(t)")

    # --- 2. Intermittent TCF (Iso) ---
    iso_avg_file, iso_tcf_df = process_and_average_time_data(
        ENSEMBLE_RANGE, BASE_PATH, PREFIX_TCF_ISO, NUM_ENSEMBLES
    )
    if iso_avg_file:
        plot_tcf_data(iso_tcf_df, 
                      "Intermittent TCF (Isolated Ta State)",
                      PLOT_TCF_ISO,
                      "C_iso(t)")

    # --- 3. Combined TCF Plot ---
    if poly_tcf_df is not None and iso_tcf_df is not None:
        plot_combined_tcf(poly_tcf_df, iso_tcf_df,
                          "Intermittent TCF Comparison",
                          PLOT_TCF_COMBINED)

    # --- 4. Cluster Size Autocorrelation ---
    size_avg_file, size_tcf_df = process_and_average_time_data(
        ENSEMBLE_RANGE, BASE_PATH, PREFIX_TCF_SIZE, NUM_ENSEMBLES
    )
    if size_avg_file:
        plot_tcf_data(size_tcf_df, 
                      "Ta Cluster Size Autocorrelation",
                      PLOT_TCF_SIZE,
                      "C_N(t)")

    # --- 5. Exchange Dynamics (Iso -> Poly) ---
    iso_poly_avg_file, iso_poly_life_df = process_and_average_lifetime_dist(
        ENSEMBLE_RANGE, BASE_PATH, PREFIX_LIFE_ISO_POLY, NUM_ENSEMBLES, LIFETIME_HIST_BINS
    )
    if iso_poly_avg_file:
        plot_lifetime_data(iso_poly_life_df,
                           "Lifetime Distribution (Isolated State)",
                           PLOT_LIFE_ISO_POLY)

    # --- 6. Exchange Dynamics (Poly -> Iso) ---
    poly_iso_avg_file, poly_iso_life_df = process_and_average_lifetime_dist(
        ENSEMBLE_RANGE, BASE_PATH, PREFIX_LIFE_POLY_ISO, NUM_ENSEMBLES, LIFETIME_HIST_BINS
    )
    if poly_iso_avg_file:
        plot_lifetime_data(poly_iso_life_df,
                           "Lifetime Distribution (Polymeric State)",
                           PLOT_LIFE_POLY_ISO)

    print("--- All dynamics analysis and plotting tasks are complete. ---")