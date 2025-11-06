import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.colors as mcolors

MAX_R = 10.0
try:
    NUMBINS = 1000
except ValueError:
    print("Warning: NUMBINS 1000  has unusual characters. Defaulting to 1000.")
    NUMBINS = 1000
    
TRAJ_INTERVAL_PS = 2.0

ENSEMBLE_START = 1
ENSEMBLE_END = 5
ENSEMBLE_RANGE = range(ENSEMBLE_START, ENSEMBLE_END + 1)
NUM_ENSEMBLES = len(ENSEMBLE_RANGE)

BASE_PATH = '../../result/structure/'
if not os.path.exists(BASE_PATH):
    print(f"Warning: Directory {BASE_PATH} does not exist. Script may fail.")

# For Gs(r,t)
TARGET_TIMES_PS = [100] # in ps

VERSION = "4"

CLUSTER_DIST_PREFIX = f"cluster_dist_ver{VERSION}"
NA_POLY_GSRT_PREFIX = f"vhf_poly_ver{VERSION}"
NA_ISO_GSRT_PREFIX  = f"vhf_iso_ver{VERSION}"
TA_POLY_MSD_PREFIX  = f"msd_ta_poly_ver{VERSION}"
TA_ISO_MSD_PREFIX   = f"msd_ta_iso_ver{VERSION}"
NA_NBO_GSRT_PREFIX  = f"vhf_na_nbo_ver{VERSION}" 
NA_BO_GSRT_PREFIX   = f"vhf_na_bo_ver{VERSION}"  
TA_O_COORD_PREFIX   = f"ta_o_coord_ver{VERSION}"
TA_BRIDGE_PREFIX    = f"ta_bridge_ver{VERSION}"
O_TYPE_DIST_PREFIX  = f"o_type_dist_ver{VERSION}"


# --- 1. Cluster Size Distribution ---
def process_and_plot_ensemble_cluster_dist(ensemble_range, base_path):
    print("--- 1. Processing Ensemble Cluster Size Distribution ---")
    all_cluster_data = []
    
    for i in ensemble_range:
        filename = os.path.join(base_path, f'{i}_{CLUSTER_DIST_PREFIX}.dat')
        if not os.path.exists(filename):
            print(f"  File not found (skipping): {filename}")
            continue
            
        try:
            data = pd.read_csv(filename, sep='\s+', comment='#',
                               names=['Timestep', 'Cluster_Size'])
            if not data.empty:
                all_cluster_data.append(data)
                print(f"  Load success: {filename} ({len(data)} data points)")
            else:
                print(f"  File empty: {filename}")
        except Exception as e:
            print(f"  File read error {filename}: {e}")

    if not all_cluster_data:
        print("--- 1. Cluster Size processing failed: No data loaded ---\n")
        return

    full_ensemble_data = pd.concat(all_cluster_data, ignore_index=True)
    output_all_file = os.path.join(base_path, 'ensemble_all_cluster_dist.dat')
    full_ensemble_data.to_csv(output_all_file, sep=' ', index=False)
    print(f"  Ensemble cluster data merged and saved: {output_all_file}")

    sizes = full_ensemble_data['Cluster_Size']
    print(f"  Total {len(sizes)} clusters (all ensembles, all timesteps).")
    print(f"  Min: {sizes.min()}, Max: {sizes.max()}, Mean: {sizes.mean():.2f}")

    plt.figure(figsize=(10, 6))
    if sizes.max() > 0 and sizes.min() > 0:
        min_val = np.floor(np.log10(sizes.min()))
        max_val = np.ceil(np.log10(sizes.max()))
        bins = np.logspace(min_val, max_val, 100) if max_val > min_val else 20
        plt.xscale('log')
        plt.yscale('log')
    else:
        bins = 100
        print("  Warning: Cluster sizes include 0 or negative values. Plotting with linear scale.")

    plt.hist(sizes, bins=bins, edgecolor='black', alpha=0.7, density=True)
    plt.title(f'Ensemble Cluster Size Distribution (N={NUM_ENSEMBLES}, All Timesteps)', fontsize=16)
    plt.xlabel('Cluster Size (Number of Atoms)', fontsize=12)
    plt.ylabel('Probability Density', fontsize=12)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)

    output_filename = 'ensemble_cluster_size_distribution.png'
    output_filename = os.path.join(base_path, output_filename)
    plt.tight_layout()
    plt.savefig(output_filename)
    print(f"  Success! Ensemble plot saved: {output_filename}")
    print("--- 1. Cluster Size processing complete ---\n")
    plt.close()


# --- 2. Ensemble Average MSD ---
def process_and_average_msd(ensemble_range, base_path, file_prefix, num_ensembles):
    print(f"--- 2. Processing Ensemble MSD ({file_prefix}) ---")
    raw_data_list = []
    min_common_length = -1

    for i in ensemble_range:
        filename = os.path.join(base_path, f'{i}_{file_prefix}.dat')
        if not os.path.exists(filename):
            print(f"  File not found (skipping): {filename}")
            continue

        try:
            data = pd.read_csv(filename, sep='\s+', comment='#', header=None)
            if data.empty:
                print(f"  File empty: {filename}")
                continue
            
            # Data: [dt, MSD, Count]
            raw_data_list.append(data.values)
            print(f"  Load success: {filename} (Shape: {data.values.shape})")
            
            current_len = len(data)
            if min_common_length == -1 or current_len < min_common_length:
                min_common_length = current_len

        except Exception as e:
            print(f"  File read error {filename}: {e}")

    if len(raw_data_list) == 0:
        print(f"--- 2. MSD processing failed: No data loaded ({file_prefix}) ---")
        return None
    elif len(raw_data_list) != num_ensembles:
        print(f"  Warning: Requested {num_ensembles} ensembles, but {len(raw_data_list)} were loaded.")

    if min_common_length == -1 or min_common_length == 0:
        print(f"--- 2. MSD processing failed: No common data length ({file_prefix}) ---")
        return None
        
    print(f"  All ensembles will be truncated to common length: {min_common_length}")

    # Truncate all datasets and stack them
    # We stack columns 1 (MSD) and 2 (Count)
    truncated_data_list = [data[:min_common_length, 1:] for data in raw_data_list]
    ensemble_stack = np.stack(truncated_data_list, axis=0) # Shape: (ensemble, time, 2)
    
    # Calculate average and SEM
    avg_data = np.mean(ensemble_stack, axis=0)
    sem_data = np.std(ensemble_stack, axis=0) / np.sqrt(ensemble_stack.shape[0])
    
    # Get the time column (dt) from the first file
    t_values_frames = raw_data_list[0][:min_common_length, 0]

    print(f"  Ensemble average calculated. (Avg Shape: {avg_data.shape}) (N={ensemble_stack.shape[0]})")

    # --- Save Average/SEM data ---
    avg_output_file = os.path.join(base_path, f'ensemble_avg_{file_prefix}.dat')
    sem_output_file = os.path.join(base_path, f'ensemble_sem_{file_prefix}.dat')
    
    # Save format: [dt, avg_MSD, avg_Count]
    avg_save_data = np.hstack((t_values_frames.reshape(-1, 1), avg_data))
    # Save format: [dt, sem_MSD, sem_Count]
    sem_save_data = np.hstack((t_values_frames.reshape(-1, 1), sem_data))
    
    header = f"# Ensemble Average MSD for {file_prefix}\n"
    header += f"# Averaged over {ensemble_stack.shape[0]} ensembles.\n"
    header += f"# Truncated to common length {min_common_length}.\n"
    header += f"# Col 0: dt (frames), Col 1: Avg_MSD (A^2), Col 2: Avg_Count"
    
    np.savetxt(avg_output_file, avg_save_data, fmt='%.8e', delimiter=' ', header=header, comments='')
    
    header = header.replace("Average", "Standard Error (SEM)").replace("Avg_", "SEM_")
    np.savetxt(sem_output_file, sem_save_data, fmt='%.8e', delimiter=' ', header=header, comments='')

    print(f"  Success! Ensemble average saved: {avg_output_file}")
    print(f"  Success! Ensemble SEM saved: {sem_output_file}")
    print(f"--- 2. MSD processing complete ({file_prefix}) ---\n")
    
    return {'avg': avg_output_file, 'sem': sem_output_file}


# --- 3. Plot Combined MSD ---
def plot_combined_msd(poly_files, iso_files, output_filename, traj_interval_ps):
    """
    Plots the averaged Polymeric and Isolated MSD data on two graphs:
    1. Linear scale (MSD vs. Time)
    2. Log-log scale (log(MSD) vs. log(Time))
    """
    print(f"--- 3. Plotting Combined Ensemble MSD ---")
    if not poly_files or not iso_files:
        print("  Error: Missing Polymeric or Isolated MSD files. Skipping plot.")
        print(f"--- 3. MSD Plotting Failed ---\n")
        return

    try:
        # Load Poly data
        poly_avg_data = pd.read_csv(poly_files['avg'], sep='\s+', comment='#', header=None).values
        poly_sem_data = pd.read_csv(poly_files['sem'], sep='\s+', comment='#', header=None).values
        t_frames = poly_avg_data[:, 0]
        t_ps = t_frames * traj_interval_ps
        
        poly_msd_avg = poly_avg_data[:, 1]
        poly_msd_sem = poly_sem_data[:, 1]
        
        # Load Iso data
        iso_avg_data = pd.read_csv(iso_files['avg'], sep='\s+', comment='#', header=None).values
        iso_sem_data = pd.read_csv(iso_files['sem'], sep='\s+', comment='#', header=None).values
        iso_t_frames = iso_avg_data[:, 0]
        
        iso_msd_avg = iso_avg_data[:, 1]
        iso_msd_sem = iso_sem_data[:, 1]

        # Check for time axis consistency
        if not np.array_equal(t_frames, iso_t_frames):
            print("  Error: Time axes for Poly and Iso MSD are different. Skipping plot.")
            print(f"--- 3. MSD Plotting Failed ---\n")
            return
            
    except Exception as e:
        print(f"  Error loading MSD data files: {e}")
        print(f"--- 3. MSD Plotting Failed ---\n")
        return

    # --- Plot 1: Linear Scale ---
    plt.figure(figsize=(10, 7))
    
    # Plot Poly
    plt.plot(t_ps, poly_msd_avg, label='Ta (Polymeric)', color='blue', linewidth=2, marker='none')
    plt.fill_between(t_ps, 
                     poly_msd_avg - poly_msd_sem, 
                     poly_msd_avg + poly_msd_sem, 
                     color='blue', alpha=0.2)
                     
    # Plot Iso
    plt.plot(t_ps, iso_msd_avg, label='Ta (Isolated)', color='red', linewidth=2, marker='none')
    plt.fill_between(t_ps, 
                     iso_msd_avg - iso_msd_sem, 
                     iso_msd_avg + iso_msd_sem, 
                     color='red', alpha=0.2)
    
    plt.title(f'Ensemble Avg Ta MSD (N={NUM_ENSEMBLES})', fontsize=16)
    plt.xlabel('Time (ps)', fontsize=14)
    plt.ylabel('MSD (Å²)', fontsize=14)
    plt.legend(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    
    linear_output_file = output_filename.replace(".png", "_linear.png")
    plt.tight_layout()
    plt.savefig(linear_output_file)
    print(f"  Success! Linear MSD plot saved: {linear_output_file}")
    plt.close()

    # --- Plot 2: Log-Log Scale ---
    plt.figure(figsize=(10, 7))

    # Filter out zero/negative values for log plot
    non_zero_mask_poly = (t_ps > 0) & (poly_msd_avg > 0)
    non_zero_mask_iso = (t_ps > 0) & (iso_msd_avg > 0)

    # Plot Poly
    plt.plot(t_ps[non_zero_mask_poly], poly_msd_avg[non_zero_mask_poly], 
             label='Ta (Polymeric)', color='blue', linewidth=2)
    
    # Plot Iso
    plt.plot(t_ps[non_zero_mask_iso], iso_msd_avg[non_zero_mask_iso], 
             label='Ta (Isolated)', color='red', linewidth=2)

    plt.title(f'Ensemble Avg Ta MSD (Log-Log Plot, N={NUM_ENSEMBLES})', fontsize=16)
    plt.xlabel('Time (ps)', fontsize=14)
    plt.ylabel('MSD (Å²)', fontsize=14)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(fontsize=12)
    plt.grid(True, which='both', linestyle='--', alpha=0.6)
    
    log_output_file = output_filename.replace(".png", "_loglog.png")
    plt.tight_layout()
    plt.savefig(log_output_file)
    print(f"  Success! Log-Log MSD plot saved: {log_output_file}")
    print(f"--- 3. MSD Plotting Complete ---\n")
    plt.close()


# --- 4. Ta-O Coordination Number Distribution ---
def process_and_plot_ta_o_coord(ensemble_range, base_path):
    """
    Loads all Ta-O coordination data, concatenates, and plots the
    overall distribution as a bar chart.
    """
    print("--- 4. Processing Ensemble Ta-O Coordination ---")
    all_ta_coord_data = []
    
    for i in ensemble_range:
        filename = os.path.join(base_path, f'{i}_{TA_O_COORD_PREFIX}.dat')
        if not os.path.exists(filename):
            print(f"  File not found (skipping): {filename}")
            continue
            
        try:
            data = pd.read_csv(filename, sep='\s+', comment='#',
                               names=['Timestep', 'Ta_Atom_Index', 'Ta_O_Coordination_Number'])
            if not data.empty:
                all_ta_coord_data.append(data)
                print(f"  Load success: {filename} ({len(data)} data points)")
            else:
                print(f"  File empty: {filename}")
        except Exception as e:
            print(f"  File read error {filename}: {e}")

    if not all_ta_coord_data:
        print("--- 4. Ta-O Coordination processing failed: No data loaded ---\n")
        return

    full_ensemble_data = pd.concat(all_ta_coord_data, ignore_index=True)
    coord_counts = full_ensemble_data['Ta_O_Coordination_Number']
    
    print(f"  Total {len(coord_counts)} Ta-atom observations (all ensembles, all timesteps).")
    
    # Calculate normalized counts for the bar plot
    cn_distribution = coord_counts.value_counts(normalize=True).sort_index()
    
    plt.figure(figsize=(10, 7))
    bars = plt.bar(cn_distribution.index, cn_distribution.values, 
                   edgecolor='black', alpha=0.8, color='c')
    
    plt.title(f'Ensemble Ta-O Coordination Distribution (N={NUM_ENSEMBLES}, All Timesteps)', fontsize=16)
    plt.xlabel('Ta-O Coordination Number (CN)', fontsize=14)
    plt.ylabel('Probability (Normalized Frequency)', fontsize=14)
    # X축 눈금을 정수로 강제
    if not cn_distribution.empty:
        plt.xticks(range(int(cn_distribution.index.min()), int(cn_distribution.index.max()) + 1))
    
    plt.grid(True, axis='y', linestyle='--', linewidth=0.5, alpha=0.7)
    
    # Add percentage labels on top of bars
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2.0, yval + 0.005,
                 f'{yval*100:.1f}%', ha='center', va='bottom')

    output_filename = 'ensemble_ta_o_coord_distribution.png'
    output_filename = os.path.join(base_path, output_filename)
    plt.tight_layout()
    plt.savefig(output_filename)
    print(f"  Success! Ensemble Ta-O CN plot saved: {output_filename}")
    print("--- 4. Ta-O Coordination processing complete ---\n")
    plt.close()


# --- 5. Ta Bridging Analysis (Polymeric Ta) ---
def process_and_plot_ta_bridging(ensemble_range, base_path):
    """
    Loads all Ta bridging data, focusing on Polymeric Ta atoms.
    Plots 1D histograms of Ta-O-Ta and Ta-Cl-Ta bridge counts.
    Plots a 2D histogram correlating O-bridges vs. Cl-bridges.
    """
    print("--- 5. Processing Ensemble Ta Bridging ---")
    all_ta_bridge_data = []
    
    for i in ensemble_range:
        filename = os.path.join(base_path, f'{i}_{TA_BRIDGE_PREFIX}.dat')
        if not os.path.exists(filename):
            print(f"  File not found (skipping): {filename}")
            continue
            
        try:
            data = pd.read_csv(filename, sep='\s+', comment='#',
                               names=['Timestep', 'Ta_Atom_Index', 'Ta_O_Ta_Bridges', 
                                      'Ta_Cl_Ta_Bridges', 'Ta_Classification'])
            if not data.empty:
                all_ta_bridge_data.append(data)
                print(f"  Load success: {filename} ({len(data)} data points)")
            else:
                print(f"  File empty: {filename}")
        except Exception as e:
            print(f"  File read error {filename}: {e}")

    if not all_ta_bridge_data:
        print("--- 5. Ta Bridging processing failed: No data loaded ---\n")
        return

    full_data = pd.concat(all_ta_bridge_data, ignore_index=True)
    
    # Filter for POLYMERIC Ta atoms only (Classification == 0)
    poly_data = full_data[full_data['Ta_Classification'] == 0].copy()
    
    if poly_data.empty:
        print("  Warning: No Polymeric Ta atoms (Classification 0) found in data. Skipping plots.")
        print("--- 5. Ta Bridging processing complete (no data) ---\n")
        return
        
    print(f"  Analyzing {len(poly_data)} Polymeric Ta-atom observations...")

    # --- Plot 1: 1D Histograms ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle(f'Ensemble Bridge Counts for Polymeric Ta (N={NUM_ENSEMBLES}, All Timesteps)', fontsize=16)

    # Plot Ta-O-Ta Bridges
    o_bridges = poly_data['Ta_O_Ta_Bridges']
    o_dist = o_bridges.value_counts(normalize=True).sort_index()
    if not o_dist.empty:
        ax1.bar(o_dist.index, o_dist.values, edgecolor='black', alpha=0.8, color='royalblue')
        ax1.set_xticks(range(int(o_dist.index.min()), int(o_dist.index.max()) + 1))
    ax1.set_title('Ta-O-Ta Bridges', fontsize=14)
    ax1.set_xlabel('Count of Ta-O-Ta Bridges per Ta', fontsize=12)
    ax1.set_ylabel('Probability (Normalized Frequency)', fontsize=12)
    ax1.grid(True, axis='y', linestyle='--', alpha=0.6)

    # Plot Ta-Cl-Ta Bridges
    cl_bridges = poly_data['Ta_Cl_Ta_Bridges']
    cl_dist = cl_bridges.value_counts(normalize=True).sort_index()
    if not cl_dist.empty:
        ax2.bar(cl_dist.index, cl_dist.values, edgecolor='black', alpha=0.8, color='seagreen')
        ax2.set_xticks(range(int(cl_dist.index.min()), int(cl_dist.index.max()) + 1))
    ax2.set_title('Ta-Cl-Ta Bridges', fontsize=14)
    ax2.set_xlabel('Count of Ta-Cl-Ta Bridges per Ta', fontsize=12)
    ax2.set_ylabel('Probability (Normalized Frequency)', fontsize=12)
    ax2.grid(True, axis='y', linestyle='--', alpha=0.6)
    
    output_filename_1d = 'ensemble_ta_bridge_1d_hists.png'
    output_filename_1d = os.path.join(base_path, output_filename_1d)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(output_filename_1d)
    print(f"  Success! Ensemble 1D Bridge Histograms saved: {output_filename_1d}")
    plt.close(fig)

    # --- Plot 2: 2D Correlation Heatmap ---
    plt.figure(figsize=(10, 8))
    
    x_data = poly_data['Ta_Cl_Ta_Bridges']
    y_data = poly_data['Ta_O_Ta_Bridges']
    
    if x_data.empty or y_data.empty:
        print("  Warning: No data for 2D Bridge Correlation plot. Skipping.")
        plt.close()
        return

    max_x = int(x_data.max())
    max_y = int(y_data.max())
    
    bins_x = np.arange(-0.5, max_x + 1.5, 1)
    bins_y = np.arange(-0.5, max_y + 1.5, 1)
    
    h, xedges, yedges, img = plt.hist2d(
        x_data, y_data, 
        bins=[bins_x, bins_y],
        cmap='viridis', 
        norm=mcolors.LogNorm() 
    )
    
    plt.colorbar(img, label='Observation Count (Log Scale)')
    plt.title(f'Correlation of O vs. Cl Bridges for Polymeric Ta (N={NUM_ENSEMBLES})', fontsize=16)
    plt.xlabel('Ta-Cl-Ta Bridges', fontsize=14)
    plt.ylabel('Ta-O-Ta Bridges', fontsize=14)
    
    plt.xticks(np.arange(0, max_x + 1, 1))
    plt.yticks(np.arange(0, max_y + 1, 1))
    
    output_filename_2d = 'ensemble_ta_bridge_2d_correlation.png'
    output_filename_2d = os.path.join(base_path, output_filename_2d)
    plt.tight_layout()
    plt.savefig(output_filename_2d)
    print(f"  Success! Ensemble 2D Bridge Correlation plot saved: {output_filename_2d}")
    print("--- 5. Ta Bridging processing complete ---\n")
    plt.close()


# --- 6. (### NEW FUNCTION ###) Oxygen Type Distribution (NBO/BO/UO) ---
def process_and_plot_oxygen_types(ensemble_range, base_path):
    """
    Loads all Oxygen Type classification data, concatenates, and plots the
    overall distribution (NBO/BO/UO) as a bar chart.
    """
    print("--- 6. Processing Ensemble Oxygen Type Distribution (NBO/BO/UO) ---")
    all_o_type_data = []
    
    for i in ensemble_range:
        filename = os.path.join(base_path, f'{i}_{O_TYPE_DIST_PREFIX}.dat')
        if not os.path.exists(filename):
            print(f"  File not found (skipping): {filename}")
            continue
            
        try:
            # C-code Format: Timestep, O_Global_Idx, O_Local_Idx, O_Type
            data = pd.read_csv(filename, sep='\s+', comment='#',
                               names=['Timestep', 'O_Global_Idx', 'O_Local_Idx', 'O_Type'])
            if not data.empty:
                all_o_type_data.append(data)
                print(f"  Load success: {filename} ({len(data)} data points)")
            else:
                print(f"  File empty: {filename}")
        except Exception as e:
            print(f"  File read error {filename}: {e}")

    if not all_o_type_data:
        print("--- 6. Oxygen Type processing failed: No data loaded ---\n")
        return

    full_data = pd.concat(all_o_type_data, ignore_index=True)
    o_type_counts = full_data['O_Type']
    
    print(f"  Total {len(o_type_counts)} Oxygen-atom observations (all ensembles, all timesteps).")
    
    # Calculate normalized counts for the bar plot
    # C-Code Defines: 0=NBO, 1=BO, 2=UO
    type_distribution = o_type_counts.value_counts(normalize=True).sort_index()
    
    # Define labels for the types
    type_labels = {
        0: 'NBO\n(1-Coord)',
        1: 'BO\n(2+-Coord)',
        2: 'UO\n(0-Coord)'
    }
    # Map index to labels, drop any types not found
    plot_labels = type_distribution.index.map(lambda x: type_labels.get(x, f'Unknown ({x})'))
    
    plt.figure(figsize=(10, 7))
    # Assign specific colors
    colors = ['tomato' if 'NBO' in l else 'forestgreen' if 'BO' in l else 'skyblue' for l in plot_labels]
    
    bars = plt.bar(plot_labels, type_distribution.values, 
                   edgecolor='black', alpha=0.8, color=colors)
    
    plt.title(f'Ensemble Oxygen Type Distribution (N={NUM_ENSEMBLES}, All Timesteps)', fontsize=16)
    plt.xlabel('Oxygen Type (Coordination to Ta)', fontsize=14)
    plt.ylabel('Probability (Normalized Frequency)', fontsize=14)
    plt.grid(True, axis='y', linestyle='--', linewidth=0.5, alpha=0.7)
    
    # Add percentage labels
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2.0, yval + 0.005,
                 f'{yval*100:.1f}%', ha='center', va='bottom', fontsize=12)

    output_filename = 'ensemble_oxygen_type_distribution.png'
    output_filename = os.path.join(base_path, output_filename)
    plt.tight_layout()
    plt.savefig(output_filename)
    print(f"  Success! Ensemble Oxygen Type plot saved: {output_filename}")
    print("--- 6. Oxygen Type processing complete ---\n")
    plt.close()


# --- 7. Ensemble Average Gs(r,t) ---
def process_and_average_van_hove(ensemble_range, base_path, file_prefix, num_ensembles):
    """
    Reads all ensemble Van Hove data, calculates average and SEM, and saves them.
    (This is a generic function for any Gs(r,t) data)
    """
    print(f"--- 7. Processing Ensemble Gs(r,t) ({file_prefix}) ---")
    ensemble_data_list = []
    t_values_frames = None
    min_common_length = -1 
    
    raw_data_list = []  
    
    for i in ensemble_range:
        filename = os.path.join(base_path, f'{i}_{file_prefix}.dat')
        if not os.path.exists(filename):
            print(f"  File not found (skipping): {filename}")
            continue

        try:
            data = pd.read_csv(filename, sep='\s+', comment='#', header=None)
            if data.empty:
                print(f"  File empty: {filename}")
                continue

            current_t_frames = data.iloc[:, 0].values
            gs_data = data.iloc[:, 1:].values
            
            raw_data_list.append((current_t_frames, gs_data, filename))
            print(f"  Load success: {filename} (Shape: {gs_data.shape})")

        except Exception as e:
            print(f"  File read error {filename}: {e}")

    if len(raw_data_list) == 0:
        print(f"--- 7. Gs(r,t) processing failed: No data loaded ({file_prefix}) ---")
        return None
    elif len(raw_data_list) != num_ensembles:
        print(f"  Warning: Requested {num_ensembles} ensembles, but {len(raw_data_list)} were loaded.")
    
    # Set base time axis (from first file)
    t_values_frames, base_gs_data, base_filename = raw_data_list[0]
    min_common_length = len(t_values_frames)
    ensemble_data_list.append(base_gs_data)
    
    # Find common length by comparing with other files
    for (current_t_frames, gs_data, filename) in raw_data_list[1:]:
        new_len = 0
        for k in range(min(min_common_length, len(current_t_frames))):
            if t_values_frames[k] == current_t_frames[k]:
                new_len = k + 1
            else:
                print(f"  Warning: {filename} differs from {base_filename} at line {k} (dt={current_t_frames[k]}).")
                break
        
        if new_len == 0:
            print(f"  CRITICAL WARNING: {filename} has no common time axis with {base_filename}. Skipping.")
        else:
            if new_len < min_common_length:
                print(f"  Info: Common time axis length reduced from {min_common_length} to {new_len}.")
                min_common_length = new_len
            ensemble_data_list.append(gs_data)
            
    if min_common_length == -1 or min_common_length == 0:
         print(f"--- 7. Gs(r,t) processing failed: No common time axis found ({file_prefix}) ---")
         return None

    print(f"  All ensembles will be truncated to common length {min_common_length} (dt = {t_values_frames[min_common_length-1]}) for averaging.")
    
    t_values_frames = t_values_frames[:min_common_length]
    
    truncated_gs_list = [gs_data[:min_common_length, :] for gs_data in ensemble_data_list if gs_data.shape[0] >= min_common_length]
    
    ensemble_stack = np.stack(truncated_gs_list, axis=0) # Shape: (ensemble, time, radius)
    
    avg_gs_data = np.mean(ensemble_stack, axis=0)
    sem_gs_data = np.std(ensemble_stack, axis=0) / np.sqrt(ensemble_stack.shape[0])
    
    print(f"  Ensemble average calculated. (Avg Shape: {avg_gs_data.shape}) (N={ensemble_stack.shape[0]})")

    # --- Save Average/SEM data ---
    avg_output_file = os.path.join(base_path, f'ensemble_avg_{file_prefix}.dat')
    sem_output_file = os.path.join(base_path, f'ensemble_sem_{file_prefix}.dat')
    
    avg_save_data = np.hstack((t_values_frames.reshape(-1, 1), avg_gs_data))
    sem_save_data = np.hstack((t_values_frames.reshape(-1, 1), sem_gs_data))
    
    header = f"# Ensemble Average Gs(r,t) for {file_prefix}\n"
    header += f"# Averaged over {ensemble_stack.shape[0]} ensembles.\n"
    header += f"# Truncated to common length {min_common_length}.\n"
    header += f"# Col 0: dt (frames), Col 1...: Gs(r,t) bins"
    
    np.savetxt(avg_output_file, avg_save_data, fmt='%.8e', delimiter=' ', header=header, comments='')
    
    header = header.replace("Average", "Standard Error (SEM)")
    np.savetxt(sem_output_file, sem_save_data, fmt='%.8e', delimiter=' ', header=header, comments='')

    print(f"  Success! Ensemble average saved: {avg_output_file}")
    print(f"  Success! Ensemble SEM saved: {sem_output_file}")
    print(f"--- 7. Gs(r,t) processing complete ({file_prefix}) ---\n")
    
    return {'avg': avg_output_file, 'sem': sem_output_file}


# --- 8. Plot Individual Gs(r,t) (2D Map + 1D Lines) ---
def plot_individual_van_hove(data_files, plot_title, output_prefix, 
                             times_to_plot_ps, traj_interval_ps, max_r, num_bins_from_c):
    """
    Reads ensemble average Gs(r,t) and SEM to create:
    1. A 2D heatmap of Gs(r,t) vs. r and t.
    2. A 1D line plot of Gs(r,t) vs. r at specific times, with error bars.
    """
    print(f"--- 8. Plotting Individual Ensemble Gs(r,t) ({output_prefix}) ---")
    
    if not data_files:
        print(f"  Error: No data files provided. Skipping plot.")
        print("--- 8. Individual Gs(r,t) Plotting Failed ---\n")
        return
        
    avg_filename = data_files['avg']
    sem_filename = data_files['sem']

    if not os.path.exists(avg_filename) or not os.path.exists(sem_filename):
        print(f"  Error: {avg_filename} or {sem_filename} not found. Skipping plot.")
        print("--- 8. Individual Gs(r,t) Plotting Failed ---\n")
        return

    try:
        avg_data_pd = pd.read_csv(avg_filename, sep='\s+', comment='#', header=None)
        t_values_frames = avg_data_pd.iloc[:, 0].values
        avg_gs_data = avg_data_pd.iloc[:, 1:].values
        sem_gs_data = pd.read_csv(sem_filename, sep='\s+', comment='#', header=None).iloc[:, 1:].values
        t_values = t_values_frames * traj_interval_ps
        
        if avg_gs_data.shape[1] != num_bins_from_c:
            print(f"  Warning: Python NUMBINS ({num_bins_from_c}) does not match data columns ({avg_gs_data.shape[1]}).")
            current_numbins = avg_gs_data.shape[1]
        else:
            current_numbins = num_bins_from_c
            
        binsize = max_r / current_numbins
        r_values = np.linspace(binsize / 2.0, max_r - binsize / 2.0, current_numbins)

        print(f"  Loaded {avg_gs_data.shape[0]} timesteps, {avg_gs_data.shape[1]} radius bins.")
    except Exception as e:
        print(f"  File read error {avg_filename} or {sem_filename}: {e}\n")
        return

    # --- 2D Heatmap Plot ---
    plt.figure(figsize=(12, 8))
    r_edges = np.linspace(0, max_r, current_numbins + 1)

    if len(t_values) > 1:
        dt_step = t_values[1] - t_values[0]
        t_edges = np.append(t_values, t_values[-1] + dt_step)
        t_edges = t_edges - (dt_step / 2.0)
        t_edges[t_edges < 0] = 0.0
    else: 
        dt_step = t_values[0] if len(t_values) > 0 and t_values[0] > 0 else traj_interval_ps
        t_edges = np.array([max(0, t_values[0] - (dt_step / 2.0)), t_values[0] + (dt_step / 2.0)])

    gs_data_transposed = avg_gs_data.T
    vmax = np.percentile(gs_data_transposed[np.isfinite(gs_data_transposed)], 99)
    if vmax <= 0: vmax = 1.0
    non_zero_finite_data = gs_data_transposed[(gs_data_transposed > 0) & np.isfinite(gs_data_transposed)]
    vmin = np.min(non_zero_finite_data) if len(non_zero_finite_data) > 0 else vmax / 1e6
    if vmin >= vmax: vmin = vmax / 1e6

    plt.pcolormesh(t_edges, r_edges, gs_data_transposed,
                   cmap='viridis', shading='flat',
                   norm=mcolors.LogNorm(vmin=vmin, vmax=vmax))

    plt.colorbar(label='Ensemble Avg Gs(r, t) [Log Scale]')
    plt.title(f'Ensemble Avg Van Hove Gs(r, t) - {plot_title} (N={NUM_ENSEMBLES})', fontsize=16)
    plt.xlabel('Time dt (ps)', fontsize=12)
    plt.ylabel('Distance r (Å)', fontsize=12)
    plt.ylim(0, 10.0)
    if t_edges.max() > 0:
        plt.xlim(0, t_edges.max())

    output_filename_2d = f'ensemble_avg_{output_prefix}_2d_heatmap.png'
    output_filename_2d = os.path.join(BASE_PATH, output_filename_2d)
    plt.tight_layout()
    plt.savefig(output_filename_2d)
    print(f"  Success! Ensemble 2D map saved: {output_filename_2d}")
    plt.close()

    # --- 1D Line Plot ---
    plt.figure(figsize=(10, 7))
    found_indices = []
    for requested_t_ps in times_to_plot_ps:
        closest_idx = np.abs(t_values - requested_t_ps).argmin()
        found_indices.append(closest_idx)
    indices_to_plot = sorted(list(set(found_indices)))

    if len(indices_to_plot) == 0:
        print("  (1D Plot) No time data found to plot. Skipping 1D plot.")
        plt.close()
        return

    print(f"  (1D Plot) Plotting data for time indices (Frame / ps):")
    colors = plt.cm.viridis(np.linspace(0, 1, len(indices_to_plot)))

    for i, idx in enumerate(indices_to_plot):
        t_ps = t_values[idx]
        t_frame = t_values_frames[idx]
        gs_at_t_avg = avg_gs_data[idx, :]
        gs_at_t_sem = sem_gs_data[idx, :]
        color = colors[i]
        print(f"    - Index {idx}: {t_frame} frames / {t_ps:.2f} ps")
        
        plt.plot(r_values, gs_at_t_avg, label=f'dt = {t_ps:.1f} ps', color=color, linewidth=2, marker='none')
        plt.fill_between(r_values, 
                         gs_at_t_avg - gs_at_t_sem, 
                         gs_at_t_avg + gs_at_t_sem, 
                         color=color, alpha=0.2)

    plt.title(f'Ensemble Avg Gs(r, t) vs r - {plot_title} (N={NUM_ENSEMBLES})', fontsize=16)
    plt.xlabel('Distance r (Å)', fontsize=12)
    plt.ylabel('Ensemble Avg Gs(r, t) (± SEM)', fontsize=12)
    plt.legend()
    plt.ylim(bottom=0)
    plt.xlim(0, 10.0)

    output_filename_1d = f'ensemble_avg_{output_prefix}_1d_lines_with_error.png'
    output_filename_1d = os.path.join(BASE_PATH, output_filename_1d)
    plt.tight_layout()
    plt.savefig(output_filename_1d)
    print(f"  Success! Ensemble 1D line plot (with error) saved: {output_filename_1d}\n")
    plt.close()


# --- 9. Plot Combined Gs(r,t) - Version 1 (4*pi*r^2, Log Scale) ---
def plot_combined_1d_gsrt_4pir2_log(data_files_1, data_files_2, 
                                  label_1, label_2, 
                                  plot_title_suffix, output_filename, 
                                  times_to_plot_ps, traj_interval_ps, max_r, num_bins_from_c,
                                  cmap_name_1='Blues', cmap_name_2='Reds'):
    """
    Combines two Gs(r,t) datasets onto a single 1D plot.
    Y-axis is 4*pi*r^2 * Gs(r,t) on a LOG scale.
    """
    print(f"--- 9. Plotting Combined Gs(r,t) [4*pi*r^2, Log] ({plot_title_suffix}) ---")
    
    # --- Load Data 1 ---
    try:
        avg_pd_1 = pd.read_csv(data_files_1['avg'], sep='\s+', comment='#', header=None)
        t_values_frames = avg_pd_1.iloc[:, 0].values
        avg_gs_data_1 = avg_pd_1.iloc[:, 1:].values
        sem_gs_data_1 = pd.read_csv(data_files_1['sem'], sep='\s+', comment='#', header=None).iloc[:, 1:].values
    except Exception as e:
        print(f"  Error: Failed to load {label_1} data files ({e}). Skipping combined plot.")
        return

    # --- Load Data 2 ---
    try:
        avg_pd_2 = pd.read_csv(data_files_2['avg'], sep='\s+', comment='#', header=None)
        t_frames_2 = avg_pd_2.iloc[:, 0].values
        avg_gs_data_2 = avg_pd_2.iloc[:, 1:].values
        sem_gs_data_2 = pd.read_csv(data_files_2['sem'], sep='\s+', comment='#', header=None).iloc[:, 1:].values
        
        if not np.array_equal(t_values_frames, t_frames_2):
             print(f"  Error: Time axes for {label_1} and {label_2} are different. Skipping plot.")
             return
    except Exception as e:
        print(f"  Error: Failed to load {label_2} data files ({e}). Skipping combined plot.")
        return

    # --- Calculate r and t axes ---
    t_values = t_values_frames * traj_interval_ps
    current_numbins = avg_gs_data_1.shape[1]
    if current_numbins != num_bins_from_c:
         print(f"  Warning: Python NUMBINS ({num_bins_from_c}) does not match data columns ({current_numbins}).")
    
    binsize = max_r / current_numbins
    r_values = np.linspace(binsize / 2.0, max_r - binsize / 2.0, current_numbins)
    
    four_pi_r_sq = 4.0 * np.pi * r_values**2
    
    plt.figure(figsize=(12, 8))
    
    found_indices = []
    for requested_t_ps in times_to_plot_ps:
        closest_idx = np.abs(t_values - requested_t_ps).argmin()
        found_indices.append(closest_idx)
    indices_to_plot = sorted(list(set(found_indices)))

    if len(indices_to_plot) == 0:
        print("  (Combined Plot) No time data found to plot.")
        plt.close()
        return

    print(f"  (Combined Plot) Plotting data for time indices (Frame / ps):")

    cmap_1 = plt.cm.get_cmap(cmap_name_1, len(indices_to_plot) + 2)
    colors_1 = [cmap_1(i) for i in np.linspace(0.4, 1.0, len(indices_to_plot))] 

    cmap_2 = plt.cm.get_cmap(cmap_name_2, len(indices_to_plot) + 2)
    colors_2 = [cmap_2(i) for i in np.linspace(0.4, 1.0, len(indices_to_plot))]

    all_positive_y_values = [] 

    for i, idx in enumerate(indices_to_plot):
        t_ps = t_values[idx]
        t_frame = t_values_frames[idx]
        
        color_1_current = colors_1[i]
        color_2_current = colors_2[i]

        print(f"    - Index {idx}: {t_frame} frames / {t_ps:.2f} ps")

        # Plot Data 1
        gs_avg_1 = avg_gs_data_1[idx, :]
        gs_sem_1 = sem_gs_data_1[idx, :]
        y_avg_1 = four_pi_r_sq * gs_avg_1
        y_sem_1 = four_pi_r_sq * gs_sem_1  
        
        plt.plot(r_values, y_avg_1, 
               label=f'{label_1} (dt = {t_ps:.1f} ps)', 
               color=color_1_current, linestyle='-', linewidth=2, marker='none')
        plt.fill_between(r_values, 
                         y_avg_1 - y_sem_1, 
                         y_avg_1 + y_sem_1, 
                         color=color_1_current, alpha=0.15)

        # Plot Data 2
        gs_avg_2 = avg_gs_data_2[idx, :]
        gs_sem_2 = sem_gs_data_2[idx, :]
        y_avg_2 = four_pi_r_sq * gs_avg_2
        y_sem_2 = four_pi_r_sq * gs_sem_2  
        
        plt.plot(r_values, y_avg_2, 
               label=f'{label_2} (dt = {t_ps:.1f} ps)', 
               color=color_2_current, linestyle='-', linewidth=2, marker='none')
        plt.fill_between(r_values, 
                         y_avg_2 - y_sem_2, 
                         y_avg_2 + y_sem_2, 
                         color=color_2_current, alpha=0.15, hatch='//') # Use hatch to differentiate

        all_positive_y_values.extend(y_avg_1[y_avg_1 > 0])
        all_positive_y_values.extend(y_avg_2[y_avg_2 > 0])

    plt.title(f'Ensemble Avg $4\\pi r^2 Gs(r, t)$ - {plot_title_suffix}', fontsize=20)
    plt.xlabel('Distance r (Å)', fontsize=20)
    plt.ylabel('$4\\pi r^2 Gs(r, t)$ (Log Scale)', fontsize=20)
    
    plt.yscale('log')
    
    plt.xlim(0, 10.0)
    if all_positive_y_values:
        min_y = np.min(all_positive_y_values) / 10.0
        max_y = np.max(all_positive_y_values) * 10.0
        plt.ylim(max(1e-3, min_y), max_y)
    else:
        plt.ylim(1e-3, 1) # Default fallback

    plt.legend(ncol=2, fontsize=16) 
    plt.tick_params(axis='both', which='major', labelsize=14) 

    plt.tight_layout()
    plt.savefig(output_filename)
    print(f"  Success! Combined Gs(r,t) plot (4pi*r^2, Log) saved: {output_filename}\n")
    plt.close()


# --- 10. (### NEW FUNCTION ###) Plot Combined Gs(r,t) - Version 2 (Original, Linear Scale) ---
def plot_combined_1d_gsrt_original(data_files_1, data_files_2, 
                                  label_1, label_2, 
                                  plot_title_suffix, output_filename, 
                                  times_to_plot_ps, traj_interval_ps, max_r, num_bins_from_c,
                                  cmap_name_1='Blues', cmap_name_2='Reds'):
    """
    Combines two Gs(r,t) datasets onto a single 1D plot.
    Y-axis is the ORIGINAL Gs(r,t) value on a LINEAR scale.
    """
    print(f"--- 10. Plotting Combined Gs(r,t) [Original, Linear] ({plot_title_suffix}) ---")
    
    # --- Load Data 1 ---
    try:
        avg_pd_1 = pd.read_csv(data_files_1['avg'], sep='\s+', comment='#', header=None)
        t_values_frames = avg_pd_1.iloc[:, 0].values
        avg_gs_data_1 = avg_pd_1.iloc[:, 1:].values
        sem_gs_data_1 = pd.read_csv(data_files_1['sem'], sep='\s+', comment='#', header=None).iloc[:, 1:].values
    except Exception as e:
        print(f"  Error: Failed to load {label_1} data files ({e}). Skipping combined plot.")
        return

    # --- Load Data 2 ---
    try:
        avg_pd_2 = pd.read_csv(data_files_2['avg'], sep='\s+', comment='#', header=None)
        t_frames_2 = avg_pd_2.iloc[:, 0].values
        avg_gs_data_2 = avg_pd_2.iloc[:, 1:].values
        sem_gs_data_2 = pd.read_csv(data_files_2['sem'], sep='\s+', comment='#', header=None).iloc[:, 1:].values
        
        if not np.array_equal(t_values_frames, t_frames_2):
             print(f"  Error: Time axes for {label_1} and {label_2} are different. Skipping plot.")
             return
    except Exception as e:
        print(f"  Error: Failed to load {label_2} data files ({e}). Skipping combined plot.")
        return

    # --- Calculate r and t axes ---
    t_values = t_values_frames * traj_interval_ps
    current_numbins = avg_gs_data_1.shape[1]
    if current_numbins != num_bins_from_c:
         print(f"  Warning: Python NUMBINS ({num_bins_from_c}) does not match data columns ({current_numbins}).")
    
    binsize = max_r / current_numbins
    r_values = np.linspace(binsize / 2.0, max_r - binsize / 2.0, current_numbins)
    
    plt.figure(figsize=(12, 8))
    
    found_indices = []
    for requested_t_ps in times_to_plot_ps:
        closest_idx = np.abs(t_values - requested_t_ps).argmin()
        found_indices.append(closest_idx)
    indices_to_plot = sorted(list(set(found_indices)))

    if len(indices_to_plot) == 0:
        print("  (Combined Plot) No time data found to plot.")
        plt.close()
        return

    print(f"  (Combined Plot) Plotting data for time indices (Frame / ps):")

    cmap_1 = plt.cm.get_cmap(cmap_name_1, len(indices_to_plot) + 2)
    colors_1 = [cmap_1(i) for i in np.linspace(0.4, 1.0, len(indices_to_plot))] 

    cmap_2 = plt.cm.get_cmap(cmap_name_2, len(indices_to_plot) + 2)
    colors_2 = [cmap_2(i) for i in np.linspace(0.4, 1.0, len(indices_to_plot))]

    for i, idx in enumerate(indices_to_plot):
        t_ps = t_values[idx]
        t_frame = t_values_frames[idx]
        
        color_1_current = colors_1[i]
        color_2_current = colors_2[i]

        print(f"    - Index {idx}: {t_frame} frames / {t_ps:.2f} ps")

        # Plot Data 1
        y_avg_1 = avg_gs_data_1[idx, :]
        y_sem_1 = sem_gs_data_1[idx, :]
        
        plt.plot(r_values, y_avg_1, 
               label=f'{label_1} (dt = {t_ps:.1f} ps)', 
               color=color_1_current, linestyle='-', linewidth=2, marker='none')
        plt.fill_between(r_values, 
                         y_avg_1 - y_sem_1, 
                         y_avg_1 + y_sem_1, 
                         color=color_1_current, alpha=0.15)

        # Plot Data 2
        y_avg_2 = avg_gs_data_2[idx, :]
        y_sem_2 = sem_gs_data_2[idx, :]
        
        plt.plot(r_values, y_avg_2, 
               label=f'{label_2} (dt = {t_ps:.1f} ps)', 
               color=color_2_current, linestyle='-', linewidth=2, marker='none')
        plt.fill_between(r_values, 
                         y_avg_2 - y_sem_2, 
                         y_avg_2 + y_sem_2, 
                         color=color_2_current, alpha=0.15, hatch='//') # Use hatch to differentiate

    plt.title(f'Ensemble Avg Original Gs(r, t) - {plot_title_suffix}', fontsize=20)
    plt.xlabel('Distance r (Å)', fontsize=20)
    plt.ylabel('Gs(r, t) (Linear Scale)', fontsize=20)
    
    plt.yscale('linear') # Use linear scale
    
    plt.xlim(0, 10.0)
    plt.ylim(bottom=0) # Start y-axis at 0

    plt.legend(ncol=2, fontsize=16) 
    plt.tick_params(axis='both', which='major', labelsize=14) 

    plt.tight_layout()
    plt.savefig(output_filename)
    print(f"  Success! Combined Gs(r,t) plot (Original, Linear) saved: {output_filename}\n")
    plt.close()


# --- Main Script Execution ---
if __name__ == "__main__":
    print("--- Ensemble Averaging and Visualization Script START ---")
    print(f"Ensemble Range: {ENSEMBLE_START} to {ENSEMBLE_END} (Total {NUM_ENSEMBLES})")
    print(f"Data Path: {BASE_PATH}\n")

#    # 1. Process Cluster Distribution
#    process_and_plot_ensemble_cluster_dist(ENSEMBLE_RANGE, BASE_PATH)
#
#    # 2. Process Ta MSD Data
#    ta_poly_files = process_and_average_msd(
#        ENSEMBLE_RANGE, BASE_PATH, TA_POLY_MSD_PREFIX, NUM_ENSEMBLES
#    )
#    ta_iso_files = process_and_average_msd(
#        ENSEMBLE_RANGE, BASE_PATH, TA_ISO_MSD_PREFIX, NUM_ENSEMBLES
#    )
#
#    # 3. Plot Combined Ta MSD
#    if ta_poly_files and ta_iso_files:
#        plot_combined_msd(
#            ta_poly_files, ta_iso_files,
#            output_filename=os.path.join(BASE_PATH, "ensemble_avg_msd_ta_combined.png"),
#            traj_interval_ps=TRAJ_INTERVAL_PS
#        )
#
#    # 4. Process Ta-O Coordination
#    process_and_plot_ta_o_coord(ENSEMBLE_RANGE, BASE_PATH)
#    
#    # 5. Process Ta Bridging
#    process_and_plot_ta_bridging(ENSEMBLE_RANGE, BASE_PATH)
#    
#    # 6. [신규 추가] Process Oxygen Type (NBO/BO/UO)
#    process_and_plot_oxygen_types(ENSEMBLE_RANGE, BASE_PATH)

    # 7. (was 4) Process Na+ Gs(r,t) (Poly vs. Iso)
    na_poly_files = process_and_average_van_hove(
        ENSEMBLE_RANGE, BASE_PATH, NA_POLY_GSRT_PREFIX, NUM_ENSEMBLES
    )
    na_iso_files = process_and_average_van_hove(
        ENSEMBLE_RANGE, BASE_PATH, NA_ISO_GSRT_PREFIX, NUM_ENSEMBLES
    )
    
    # 8. (was 5) Process Na+ Gs(r,t) (NBO vs. BO)
    na_nbo_files = process_and_average_van_hove(
        ENSEMBLE_RANGE, BASE_PATH, NA_NBO_GSRT_PREFIX, NUM_ENSEMBLES
    )
    na_bo_files = process_and_average_van_hove(
        ENSEMBLE_RANGE, BASE_PATH, NA_BO_GSRT_PREFIX, NUM_ENSEMBLES
    )

    # 9. (was 6) Plot Individual Gs(r,t) (Poly)
    if na_poly_files:
        plot_individual_van_hove(
            na_poly_files,
            plot_title="Na+ (Polymeric)",
            output_prefix="vhf_na_poly",
            times_to_plot_ps=TARGET_TIMES_PS,
            traj_interval_ps=TRAJ_INTERVAL_PS, max_r=MAX_R, num_bins_from_c=NUMBINS
        )

    # 10. (was 7) Plot Individual Gs(r,t) (Iso)
    if na_iso_files:
        plot_individual_van_hove(
            na_iso_files,
            plot_title="Na+ (Isolated)",
            output_prefix="vhf_na_iso",
            times_to_plot_ps=TARGET_TIMES_PS,
            traj_interval_ps=TRAJ_INTERVAL_PS, max_r=MAX_R, num_bins_from_c=NUMBINS
        )

    # 11. (was 8) Plot Individual Gs(r,t) (NBO)
    if na_nbo_files:
        plot_individual_van_hove(
            na_nbo_files,
            plot_title="Na+ (near NBO/UO)",
            output_prefix="vhf_na_nbo",
            times_to_plot_ps=TARGET_TIMES_PS,
            traj_interval_ps=TRAJ_INTERVAL_PS, max_r=MAX_R, num_bins_from_c=NUMBINS
        )

    # 12. (was 9) Plot Individual Gs(r,t) (BO)
    if na_bo_files:
        plot_individual_van_hove(
            na_bo_files,
            plot_title="Na+ (near BO)",
            output_prefix="vhf_na_bo",
            times_to_plot_ps=TARGET_TIMES_PS,
            traj_interval_ps=TRAJ_INTERVAL_PS, max_r=MAX_R, num_bins_from_c=NUMBINS
        )

    # 13. (was 10) Plot Combined Gs(r,t) (Poly vs. Iso)
    if na_poly_files and na_iso_files:
        plot_combined_1d_gsrt_4pir2_log(
            na_poly_files, na_iso_files,
            label_1="Na+ (Poly)", label_2="Na+ (Iso)",
            plot_title_suffix="Na+ (Poly vs. Iso)",
            output_filename=os.path.join(BASE_PATH, "ensemble_avg_gsrt_na_poly_iso_4pir2_log.png"),
            times_to_plot_ps=TARGET_TIMES_PS,
            traj_interval_ps=TRAJ_INTERVAL_PS, max_r=MAX_R, num_bins_from_c=NUMBINS,
            cmap_name_1='Blues', cmap_name_2='Reds'
        )
        plot_combined_1d_gsrt_original(
            na_poly_files, na_iso_files,
            label_1="Na+ (Poly)", label_2="Na+ (Iso)",
            plot_title_suffix="Na+ (Poly vs. Iso)",
            output_filename=os.path.join(BASE_PATH, "ensemble_avg_gsrt_na_poly_iso_original_linear.png"),
            times_to_plot_ps=TARGET_TIMES_PS,
            traj_interval_ps=TRAJ_INTERVAL_PS, max_r=MAX_R, num_bins_from_c=NUMBINS,
            cmap_name_1='Blues', cmap_name_2='Reds'
        )

    # 14. (was 11) Plot Combined Gs(r,t) (NBO vs. BO)
    if na_nbo_files and na_bo_files:
        plot_combined_1d_gsrt_4pir2_log(
            na_nbo_files, na_bo_files,
            label_1="Na+ (NBO/UO)", label_2="Na+ (BO)",
            plot_title_suffix="Na+ (NBO/UO vs. BO)",
            output_filename=os.path.join(BASE_PATH, "ensemble_avg_gsrt_na_nbo_bo_4pir2_log.png"),
            times_to_plot_ps=TARGET_TIMES_PS,
            traj_interval_ps=TRAJ_INTERVAL_PS, max_r=MAX_R, num_bins_from_c=NUMBINS,
            cmap_name_1='Greens', cmap_name_2='Oranges'
        )
        plot_combined_1d_gsrt_original(
            na_nbo_files, na_bo_files,
            label_1="Na+ (NBO/UO)", label_2="Na+ (BO)",
            plot_title_suffix="Na+ (NBO/UO vs. BO)",
            output_filename=os.path.join(BASE_PATH, "ensemble_avg_gsrt_na_nbo_bo_original_linear.png"),
            times_to_plot_ps=TARGET_TIMES_PS,
            traj_interval_ps=TRAJ_INTERVAL_PS, max_r=MAX_R, num_bins_from_c=NUMBINS,
            cmap_name_1='Greens', cmap_name_2='Oranges'
        )


    if na_poly_files and na_bo_files:
        print("\n>>> EXECUTING NEW PLOT 15: Na+ (Poly) vs. Na+ (BO) [4pir2 Log] <<<")
        plot_combined_1d_gsrt_4pir2_log(
            na_poly_files, na_bo_files,
            label_1="Na+ (Poly)", label_2="Na+ (BO)",
            plot_title_suffix="Na+ (Poly vs. BO)",
            output_filename=os.path.join(BASE_PATH, "ensemble_avg_gsrt_na_poly_bo_4pir2_log.png"),
            times_to_plot_ps=TARGET_TIMES_PS,
            traj_interval_ps=TRAJ_INTERVAL_PS, max_r=MAX_R, num_bins_from_c=NUMBINS,
            cmap_name_1='Blues', cmap_name_2='Oranges'
        )

    print("--- All Ensemble Averaging and Visualization Tasks COMPLETE ---")