import numpy as np
import matplotlib.pyplot as plt
import os


NUM_ENSEMBLES = 5

TIME_PER_FRAME_PS = 2.0 

def load_msd_data(file_list):
    all_msd_data = []
    time_data = None
    
    for f_path in file_list:
        try:
            data = np.loadtxt(f_path, comments='#', usecols=(0, 1))
            
            if data.size == 0:
                print(f"  [Warning] File is empty: {f_path}. Skipping.")
                continue

            if time_data is None:
                time_data = data[:, 0] * TIME_PER_FRAME_PS
                
            # Store MSD data
            all_msd_data.append(data[:, 1])
            
        except OSError as e:
            print(f"  [Error] Could not read file: {e}")
        except IndexError as e:
            print(f"  [Error] File format error in {f_path}: {e}. Is the file empty?")
            
    if not all_msd_data:
        print("  [Fatal] No data was loaded. Exiting.")
        return None, None
        
    print(f"  Successfully loaded {len(all_msd_data)} files.")

    try:
        min_len = min(len(msd) for msd in all_msd_data)
    except ValueError:
        print("  [Fatal] No valid data was loaded at all.")
        return None, None
        
    if time_data is not None:
        time_data = time_data[:min_len]
    msd_array_truncated = [msd[:min_len] for msd in all_msd_data]
    
    return time_data, np.array(msd_array_truncated)

def calculate_average_msd(msd_array):
    if msd_array is None or msd_array.size == 0:
        return None, None
        
    n_ensembles = msd_array.shape[0]
    
    avg_msd = np.mean(msd_array, axis=0)
    

    sem_msd = np.std(msd_array, axis=0) / np.sqrt(n_ensembles)
    
    print(f"  Calculated average over {n_ensembles} ensembles.")
    return avg_msd, sem_msd

def plot_msd(traces, output_filename, title):

    #plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(14, 9))

    print(f"  Generating plot: {output_filename}")
    
    for trace in traces:
        t = trace['time']
        msd = trace['msd']
        sem = trace['sem']
        label = trace['label']
        
        ax.plot(t, msd, label=label, color=trace['color'], linewidth=2)
        
        if sem is not None:
            ax.fill_between(t, msd - sem, msd + sem, 
                            color=trace['color'], alpha=0.2, 
                            label=f"{label} (SEM)")
    
    ax.set_xlabel("Time (ps)") 
    ax.set_ylabel(r"MSD ($\AA^2$)") 
    ax.set_title(title)
    ax.legend()
    

    plt.tight_layout()
    
    plt.savefig(output_filename, dpi=300)
    plt.close(fig)

def main():
    base_path = "../../result/structure"
    
    poly_template = "{index}_msd_ta_poly_ver2.dat"
    iso_template = "{index}_msd_ta_iso_ver2.dat"
    
    os.makedirs(base_path, exist_ok=True)
    
    poly_files = []
    iso_files = []

    for i in range(1, NUM_ENSEMBLES + 1):
        poly_name = poly_template.format(index=i)
        iso_name = iso_template.format(index=i)
        
        poly_files.append(os.path.join(base_path, poly_name))
        iso_files.append(os.path.join(base_path, iso_name))
    
    print("--- MSD Plotter Initialized ---")
    print(f"Base path: {os.path.abspath(base_path)}")
    print(f"Ensembles: {NUM_ENSEMBLES}")
    print(f"Time per frame: {TIME_PER_FRAME_PS} ps")
    
    print("\nTarget files for 'Polymeric':")
    for f in poly_files:
        print(f"  - {f}")
    
    print("\nTarget files for 'Isolated':")
    for f in iso_files:
        print(f"  - {f}")

    
    # --- 1. Process Polymeric Ta Data ---
    print("\nProcessing Polymeric Ta MSD data...")
    t_poly, msd_poly_all = load_msd_data(poly_files)
    avg_poly, sem_poly = calculate_average_msd(msd_poly_all)
    
    # --- 2. Process Isolated Ta Data ---
    print("\nProcessing Isolated Ta MSD data...")
    t_iso, msd_iso_all = load_msd_data(iso_files)
    avg_iso, sem_iso = calculate_average_msd(msd_iso_all)
    
    # Check if data loading was successful
    data_loaded_poly = t_poly is not None and avg_poly is not None
    data_loaded_iso = t_iso is not None and avg_iso is not None

    if not data_loaded_poly and not data_loaded_iso:
        print("\n[Fatal] No data was loaded for either type. Exiting.")
        return

    print("\nGenerating plots...")
    
    # --- 3. Generate Plots (Save to base_path) ---
    
    traces_to_plot = []

    if data_loaded_poly:
        trace_poly = {
            'time': t_poly, 
            'msd': avg_poly, 
            'sem': sem_poly, 
            'label': 'Polymeric Ta', 
            'color': 'crimson'
        }
        plot_msd([trace_poly], 
                 os.path.join(base_path, 'msd_ta_poly_average.png'), 
                 'Ensemble-Averaged MSD (Polymeric Ta)')
        traces_to_plot.append(trace_poly)
    else:
        print("  Skipping Polymeric plot (no data).")

    if data_loaded_iso:
        trace_iso = {
            'time': t_iso, 
            'msd': avg_iso, 
            'sem': sem_iso, 
            'label': 'Isolated Ta', 
            'color': 'royalblue'
        }
        plot_msd([trace_iso], 
                 os.path.join(base_path, 'msd_ta_iso_average.png'), 
                 'Ensemble-Averaged MSD (Isolated Ta)')
        traces_to_plot.append(trace_iso)
    else:
        print("  Skipping Isolated plot (no data).")

    
    # Plot 3: Combined Plot (if both were loaded)
    if traces_to_plot:
        plot_msd(traces_to_plot, 
                 os.path.join(base_path, 'msd_ta_combined_average.png'), 
                 'Ensemble-Averaged MSD (Polymeric vs. Isolated Ta)')
    
    print(f"\nAll tasks completed. Check for .png files in '{base_path}'")

if __name__ == "__main__":
    main()