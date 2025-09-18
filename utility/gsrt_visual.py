import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import re

def parse_gsrt_data(filename):
    all_t = []
    all_r = []
    all_gs = []
    
    current_t_data = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                if current_t_data:
                    r_vals, gs_vals = zip(*current_t_data)
                    all_r.append(r_vals)
                    all_gs.append(gs_vals)
                    current_t_data = []
                continue

            if line.startswith('# t_delta'):
                match = re.search(r'=\s*(\d+)', line)
                if match:
                    all_t.append(int(match.group(1)))
                continue
            
            if line.startswith('#'):
                continue

            parts = line.split()
            try:
                r = float(parts[0])
                gs = float(parts[1])
                current_t_data.append((r, gs))
            except (ValueError, IndexError):
                print(f"Warning: Could not parse line: {line}")

    if current_t_data:
        r_vals, gs_vals = zip(*current_t_data)
        all_r.append(r_vals)
        all_gs.append(gs_vals)

    if not all_t:
        raise ValueError("No data parsed. Check file format and content.")

    t_coords = np.array(all_t)
    r_coords = np.array(all_r[0])
    gs_grid = np.array(all_gs).T
    
    return r_coords, t_coords, gs_grid

def plot_gsrt_heatmap(r, t, gs, species_name, output_filename="../result/van_hove_gsrt.png"):

    vmin = np.min(gs[gs > 0]) if np.any(gs > 0) else 1e-10
    
    fig, ax = plt.subplots(figsize=(10, 8))

    t_step = t[1] - t[0] if len(t) > 1 else 1
    r_step = r[1] - r[0] if len(r) > 1 else 1
    
    t_edges = np.append(t, t[-1] + t_step)
    r_edges = np.append(r, r[-1] + r_step)

    T, R = np.meshgrid(t_edges, r_edges)
    
    c = ax.pcolormesh(T, R, gs, 
                      shading='auto', 
                      cmap='viridis', 
                      norm=LogNorm(vmin=vmin, vmax=gs.max()))
    
    cbar = fig.colorbar(c, ax=ax)
    cbar.set_label(r'$G_s(r, t) \quad [\AA^{-3}]$', fontsize=14)

    ax.set_xlabel(r'Time (ps))', fontsize=14)
    ax.set_ylabel('Distance, $r$ (Å)', fontsize=14)
    ax.set_title(f'Self-part of Van Hove Correlation Function $G_s(r, t)$ for {species_name}', fontsize=16)

    ax.set_ylim(0, r.max())
    ax.set_xlim(t.min(), t.max())
    
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    print(f"Plot saved to {output_filename}")
    plt.show()


if __name__ == '__main__':
    SPECIES = 'Na'
    FILENAME = f'../result/gsrt_all_{SPECIES}.dat'
    
    try:
        r_coords, t_coords, gs_grid = parse_gsrt_data(FILENAME)
        print("Data parsing complete.")
        print(f"Data grid shape (r, t): {gs_grid.shape}")

        plot_gsrt_heatmap(r_coords, t_coords, gs_grid, SPECIES)
        
    except FileNotFoundError:
        print(f"Error: The file '{FILENAME}' was not found.")
        print("Please make sure the Fortran analysis has been run and the output file is in the same directory.")
    except Exception as e:
        print(f"An error occurred: {e}")