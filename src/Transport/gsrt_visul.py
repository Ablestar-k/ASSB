import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.colors as mcolors
from cycler import cycler

SPECIES = 'Na'
INPUT_FILE = f'../../result/gsrt_ensemble_average_{SPECIES}.dat'
SELECTED_TIMES_PS = [4.0, 40.0, 400.0, 1000.0, 2000.0, 4000.0]

def parse_gsrt_data(filepath):
    data = {}
    current_time = None
    time_interval = None

    print(f"Parsing file: '{filepath}'...")
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if '# t_delta =' in line and time_interval is None:
                    match = re.search(r'# t_delta = \s*(\d+)', line)
                    if match:
                        first_time_val = float(match.group(1))
                        time_interval = first_time_val
                        print(f"Estimated Time Interval (T_INTERVAL): {time_interval} ps")
                        break

            for line in lines:
                line = line.strip()
                if line.startswith('# t_delta ='):
                    match = re.search(r'# t_delta = \s*([\d.]+)', line)
                    if match:
                        current_time = float(match.group(1))
                        data[current_time] = []
                elif line.startswith('# r (Angstrom)') or not line:
                    continue
                elif current_time is not None:
                    parts = line.split()
                    if len(parts) == 3:
                        data[current_time].append([float(p) for p in parts])

    except FileNotFoundError:
        print(f"Error: File not found at '{filepath}'.")
        return None, None
    except Exception as e:
        print(f"An error occurred while parsing the file: {e}")
        return None, None

    for t, values in data.items():
        if not values: continue
        arr = np.array(values)
        data[t] = {
            'r': arr[:, 0],
            'gs_mean': arr[:, 1],
            'gs_std': arr[:, 2]
        }
    
    sorted_times = sorted(data.keys())
    sorted_data = {t: data[t] for t in sorted_times}

    print("File parsing complete.")
    return sorted_data, time_interval

def plot_gsrt_vs_r(data, selected_times, species_name, figname="../../result/gs_vs_r.png"):
    if not data:
        print("No data available to plot.")
        return

    #plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(14, 9))
    
    num_lines = len(selected_times)
    ax.set_prop_cycle(cycler('color', plt.cm.viridis(np.linspace(0, 1, num_lines))))

    for t in selected_times:
        if t in data:
            r = data[t]['r']
            gs_mean = data[t]['gs_mean']
            gs_std = data[t]['gs_std']
            
            ax.plot(r, gs_mean, label=f't = {t} ps')
            ax.fill_between(r, gs_mean - gs_std, gs_mean + gs_std, alpha=0.2)
        else:
            print(f"Warning: No data found for selected time t = {t} ps.")

    ax.set_xlabel('r (Å)')
    ax.set_ylabel(r'$G_s(r, t) \quad (\AA^{-3})$')
    ax.set_title(f'Self-Part of Van Hove Correlation Function for {species_name}')
    ax.legend(title="Time Delay")
    ax.set_xlim(0, 5)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    plt.savefig(figname, dpi=300)
    print(f"Plot saved successfully: '{figname}'")
    plt.close(fig)


def plot_gsrt_heatmap(data, species_name, figname="../../result/gs_heatmap.png"):
    if not data:
        print("No data available to plot heatmap.")
        return
        
    times = np.array(sorted(data.keys()))
    r_values = data[times[0]]['r']
    
    gs_grid = np.zeros((len(r_values), len(times)))
    
    for i, t in enumerate(times):
        gs_grid[:, i] = data[t]['gs_mean']

    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(12, 8)) 
    
    epsilon = 1e-10
    im = ax.pcolormesh(times, r_values, gs_grid + epsilon, shading='gouraud', cmap='viridis', 
                       norm=mcolors.LogNorm(vmin=np.min(gs_grid[gs_grid>0]), vmax=np.max(gs_grid)))
    
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('r (Å)')
    ax.set_title(f'Heatmap of $G_s(r, t)$ for {species_name}')

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r'$G_s(r, t) \quad (\AA^{-3})$')
    
    plt.tight_layout()
    plt.savefig(figname, dpi=300)
    print(f"Heatmap saved successfully: '{figname}'")
    plt.close(fig) 

def plot_gsrt_4pir2_vs_r_logy(data, selected_times, species_name, figname="../../result/gs_4pir2_logy_vs_r.png"):
    if not data:
        print("No data available to plot.")
        return

    #plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(18, 10))
    
    num_lines = len(selected_times)
    ax.set_prop_cycle(cycler('color', plt.cm.viridis(np.linspace(0, 1, num_lines))))

    for t in selected_times:
        if t in data:
            r = data[t]['r']
            gs_mean = data[t]['gs_mean']
            gs_std = data[t]['gs_std']
            
            four_pi_r2 = 4.0 * np.pi * r**2
            
         
            y_mean = four_pi_r2 * gs_mean
            y_std_band = four_pi_r2 * gs_std 
            
          
            valid_indices = r > 1e-6 
            r_plot = r[valid_indices]
            y_mean_plot = y_mean[valid_indices]
            y_std_band_plot = y_std_band[valid_indices]

            if len(r_plot) > 0:
                ax.plot(r_plot, y_mean_plot, label=f't = {t} ps')
                ax.fill_between(r_plot, y_mean_plot - y_std_band_plot, y_mean_plot + y_std_band_plot, alpha=0.2)
  
        else:
            print(f"Warning: No data found for selected time t = {t} ps.")

    ax.set_xlabel('r (Å)')
    ax.set_ylabel(r'$4 \pi r^2 G_s(r, t)$') 
    ax.set_title(r'$4 \pi r^2 G_s(r, t)$ vs r for ' + f'{species_name}')
    ax.legend()
    #ax.set_xlim(0, 5)
    
    ax.set_yscale('log') 

    plt.tight_layout()
    plt.savefig(figname, dpi=300)
    print(f"Plot saved successfully: '{figname}'")
    plt.close(fig)


if __name__ == '__main__':

    gsrt_data, t_interval = parse_gsrt_data(INPUT_FILE)
    
    if gsrt_data:
        plot_gsrt_vs_r(gsrt_data, SELECTED_TIMES_PS, SPECIES)
        
        plot_gsrt_4pir2_vs_r_logy(gsrt_data, SELECTED_TIMES_PS, SPECIES)

        # plot_gsrt_heatmap(gsrt_data, SPECIES)