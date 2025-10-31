import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import seaborn as sns
import pandas as pd

def load_gsrt_data(filepath: str):
    """
    Loads Gs(r, t) data from a specified file path.
    The file is expected to have data blocks for each time step,
    prefixed by a '# t_delta = ... ps' header.
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"ERROR: File not found: '{filepath}'")
        return None, None, None

    all_t = []
    all_r = []
    all_gs_data_blocks = []
    
    current_gs_block = []
    current_r_block = []
    r_values_captured = False # Flag to read 'r' values only from the first block

    for line in lines:
        line = line.strip()
        
        # Check for the start of a new time block
        if line.startswith('# t_delta'):
            # If we have a block, save it
            if current_gs_block:
                all_gs_data_blocks.append(current_gs_block)
                # If this was the first block, save the r values
                if not r_values_captured:
                    all_r = np.array(current_r_block)
                    r_values_captured = True
            
            # Start new block
            current_gs_block = []
            if not r_values_captured:
                current_r_block = []
            
            # Parse the time value
            try:
                value_str = line.split('=')[-1]
                value_str = value_str.replace('ps', '')
                current_t_val = float(value_str.strip())
                all_t.append(current_t_val)
            except (IndexError, ValueError):
                print(f"WARNING: Could not parse time value: {line}")
                
        # Read data lines (non-comment, non-empty)
        elif line and not line.startswith('#'):
            try:
                parts = line.split()
                if len(parts) >= 2:
                    if not r_values_captured:
                        current_r_block.append(float(parts[0])) # r
                    current_gs_block.append(float(parts[1]))  # Gs(r,t)
                else:
                    print(f"WARNING: Skipping malformed data line: {line}")
            except (IndexError, ValueError):
                 print(f"WARNING: Could not parse data line: {line}")

    # Append the very last block
    if current_gs_block:
        all_gs_data_blocks.append(current_gs_block)
        if not r_values_captured:
            all_r = np.array(current_r_block)
            r_values_captured = True

    if not all_t or all_r.size == 0:
        print("ERROR: Failed to read valid data blocks.")
        return None, None, None
        
    t_values = np.array(all_t)
    gs_matrix = np.array(all_gs_data_blocks)
    
    # Validation check
    if gs_matrix.shape != (len(t_values), len(all_r)):
        print(f"ERROR: Data shape mismatch.")
        print(f"       Gs(r,t) matrix: {gs_matrix.shape}")
        print(f"       Time axis (t): {len(t_values)}, Distance axis (r): {len(all_r)}")
        # Try to fix common mismatch (e.g., last block was shorter)
        if gs_matrix.shape[0] == len(t_values) -1:
             print("       Mismatch: 1 fewer Gs block than t values. Dropping last t value.")
             t_values = t_values[:-1]
        elif gs_matrix.shape[0] -1 == len(t_values):
             print("       Mismatch: 1 more Gs block than t values. Dropping last Gs block.")
             gs_matrix = gs_matrix[:-1]
        else:
             print("       Cannot resolve mismatch. Aborting.")
             return None, None, None

    return all_r, t_values, gs_matrix

def plot_gsrt_heatmap_seaborn(r_values, t_values, gs_matrix, output_path, species_name):
    """
    Visualizes Gs(r, t) using a Seaborn heatmap with a robust LogNorm scale.
    (Original Version)
    """
    print(f"Generating Gs(r,t) Seaborn heatmap at '{output_path}'...")
    
    # 1. Create DataFrame
    try:
        # Transpose gs_matrix so that r is index (rows) and t is columns (cols)
        df_gsrt = pd.DataFrame(gs_matrix.T, index=r_values, columns=t_values)
    except ValueError as e:
        print(f"ERROR: Failed to create DataFrame. Dimension mismatch (r, t, gs_matrix): {e}")
        return

    # 2. Set up Robust LogNorm
    # Find all positive values to avoid log(0)
    positive_gs = df_gsrt.values[df_gsrt.values > 1e-12] 
    if positive_gs.size == 0:
        print("WARNING: No Gs(r,t) data > 0 found. Plotting with linear scale.")
        norm = None
    else:
        # Use percentiles for a robust min/max, avoiding extreme outliers
        robust_vmin = np.percentile(positive_gs, 1.0)
        robust_vmax = np.percentile(positive_gs, 99.9)
        
        if robust_vmin >= robust_vmax:
             robust_vmin = positive_gs.min()
             robust_vmax = positive_gs.max()

        print(f"  Adjusting color scale (Robust LogNorm, 1%-99.9%): vmin={robust_vmin:.2e}, vmax={robust_vmax:.2e}")
        norm = LogNorm(vmin=robust_vmin, vmax=robust_vmax)

    # 3. Calculate tick label spacing
    x_ticks_skip = max(1, df_gsrt.shape[1] // 10) # Show ~10 labels on x-axis
    y_ticks_skip = max(1, df_gsrt.shape[0] // 10) # Show ~10 labels on y-axis

    # 4. Plotting
    sns.set_theme(style="white", context="paper") 
    fig, ax = plt.subplots(figsize=(10, 8))
    
    sns.heatmap(
        df_gsrt,
        ax=ax,
        cmap="rocket_r",    # Use a perceptually uniform colormap (reversed)
        norm=norm,          # Apply the log normalization
        cbar_kws={'label': r'Gs(r, t) ($\AA^{-3}$)'}, # Colorbar label
        xticklabels=x_ticks_skip,  # Sparse x-ticks
        yticklabels=y_ticks_skip,  # Sparse y-ticks
    )

    ax.invert_yaxis() # Put r=0 at the bottom
    
    # Format ticks
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

    # Labels and Title
    ax.set_xlabel('Time, t (ps)', fontsize=12)
    ax.set_ylabel('Distance, r ($\AA$)', fontsize=12)
    ax.set_title(f'Van Hove Self-Correlation, $G_s(r, t)$ - Species: {species_name}', fontsize=14)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close(fig)
    print("Gs(r,t) Seaborn heatmap generation complete.")

def plot_gsrt_heatmap_4pir2_seaborn(r_values, t_values, gs_matrix, output_path, species_name):
    """
    Visualizes 4*pi*r^2 * Gs(r, t) using a Seaborn heatmap with a robust LogNorm scale.
    (New Version)
    """
    print(f"Generating 4*pi*r^2 Gs(r, t) Seaborn heatmap at '{output_path}'...")

    # 1. Calculate 4*pi*r^2 * Gs(r, t)
    print(f"  Calculating 4*pi*r^2 * Gs(r, t)...")
    # r_values shape is (n_r,)
    # gs_matrix shape is (n_t, n_r)
    if r_values.shape[0] != gs_matrix.shape[1]:
        print(f"ERROR: r_values dimension ({r_values.shape[0]}) does not match gs_matrix dim 1 ({gs_matrix.shape[1]})")
        return
        
    # Create the 4*pi*r^2 prefactor, shape (n_r,)
    prefactor = 4 * np.pi * (r_values**2) 
    
    # Broadcast prefactor (n_r,) across gs_matrix (n_t, n_r)
    gs_r_weighted = gs_matrix * prefactor 
    # gs_r_weighted shape is (n_t, n_r)

    # 2. Create DataFrame
    try:
        # Transpose to (n_r, n_t) for the DataFrame
        df_gsrt = pd.DataFrame(gs_r_weighted.T, index=r_values, columns=t_values)
    except ValueError as e:
        print(f"ERROR: Failed to create DataFrame. Dimension mismatch: {e}")
        return

    # 3. Set up Robust LogNorm (same logic, applied to new data)
    positive_gs = df_gsrt.values[df_gsrt.values > 1e-12] 
    if positive_gs.size == 0:
        print("WARNING: No 4*pi*r^2 Gs(r,t) data > 0 found. Plotting with linear scale.")
        norm = None
    else:
        robust_vmin = np.percentile(positive_gs, 1.0)
        robust_vmax = np.percentile(positive_gs, 99.9)
        
        if robust_vmin >= robust_vmax:
             robust_vmin = positive_gs.min()
             robust_vmax = positive_gs.max()

        print(f"  Adjusting color scale (Robust LogNorm, 1%-99.9%): vmin={robust_vmin:.2e}, vmax={robust_vmax:.2e}")
        norm = LogNorm(vmin=robust_vmin, vmax=robust_vmax)

    # 4. Calculate tick label spacing
    x_ticks_skip = max(1, df_gsrt.shape[1] // 10) 
    y_ticks_skip = max(1, df_gsrt.shape[0] // 10)

    # 5. Plotting
    sns.set_theme(style="white", context="paper") 
    fig, ax = plt.subplots(figsize=(10, 8))
    
    sns.heatmap(
        df_gsrt,
        ax=ax,
        cmap="rocket_r", 
        norm=norm,
        # *** CHANGED CBAR LABEL ***
        cbar_kws={'label': r'4$\pi$r$^2$ Gs(r, t) ($\AA^{-1}$)'}, 
        xticklabels=x_ticks_skip,
        yticklabels=y_ticks_skip,
    )

    ax.invert_yaxis()
    
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

    ax.set_xlabel('Time, t (ps)', fontsize=12)
    ax.set_ylabel('Distance, r ($\AA$)', fontsize=12)
    # *** CHANGED TITLE ***
    ax.set_title(f'Van Hove Self-Correlation (Radial), $4\pi r^2 G_s(r, t)$ - Species: {species_name}', fontsize=14)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close(fig)
    print("4*pi*r^2 Seaborn heatmap generation complete.")


def main():
    SPECIES = 'Na' 
    BASE_PATH = '../../result/'
    INPUT_FILE = os.path.join(BASE_PATH, f'gsrt_all_{SPECIES}.dat')
    
    # *** Define TWO output files ***
    OUTPUT_FILE_GS = os.path.join(BASE_PATH, f'gsrt_heatmap_{SPECIES}.png')
    OUTPUT_FILE_4PIR2 = os.path.join(BASE_PATH, f'gsrt_heatmap_4pir2_{SPECIES}.png')

    print(f"Reading data file: {INPUT_FILE}")
    r, t, gs_matrix = load_gsrt_data(INPUT_FILE)
    
    if r is None:
        print("Data loading failed. Exiting program.")
        return

    print(f"Data loaded successfully: {gs_matrix.shape[0]} time steps, {gs_matrix.shape[1]} distance bins")
    
    # *** Call BOTH plotting functions ***
    
    # 1. Plot original Gs(r, t)
    plot_gsrt_heatmap_seaborn(r, t, gs_matrix, OUTPUT_FILE_GS, SPECIES)
    
    # 2. Plot 4*pi*r^2 * Gs(r, t)
    plot_gsrt_heatmap_4pir2_seaborn(r, t, gs_matrix, OUTPUT_FILE_4PIR2, SPECIES)

if __name__ == "__main__":
    main()