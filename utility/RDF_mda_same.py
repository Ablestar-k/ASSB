import os
import MDAnalysis as mda
from MDAnalysis.analysis import rdf
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

BASE_DIR_NAME = "../dump/dump_NTOC_ver1_{}"
NUM_ENSEMBLES = 5
ATOM_GROUP_TO_ANALYZE = "name Cl"
NBINS = 400
RANGE = (0.0, 20.0)
SKIP_PERCENT = 0.0
STEP_SIZE = 1

def get_box_matrix_from_xyz(filename):
    with open(filename, 'r') as f:
        f.readline()
        second_line = f.readline()
        try:
            lattice_str = second_line.split('Lattice="')[1].split('"')[0]
            vecs = np.array([float(x) for x in lattice_str.split()])
            return vecs.reshape(3, 3)
        except (IndexError, ValueError):
            return None

def convert_matrix_to_dimensions(matrix):
    a, b, c = matrix[0], matrix[1], matrix[2]
    lx, ly, lz = np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c)
    epsilon = 1e-10
    alpha = np.rad2deg(np.arccos(np.dot(b, c) / ((ly * lz) + epsilon)))
    beta = np.rad2deg(np.arccos(np.dot(a, c) / ((lx * lz) + epsilon)))
    gamma = np.rad2deg(np.arccos(np.dot(a, b) / ((lx * ly) + epsilon)))
    return np.array([lx, ly, lz, alpha, beta, gamma])

all_rdfs = []
distance_bins = None

print("Starting analysis using InterRDF with exclusion_block...")

for i in range(1, NUM_ENSEMBLES + 1):
    directory = BASE_DIR_NAME.format(i)
    xyz_filename = f"{i}_product.xyz"
    # =====================
    # Edit when you calculate production run
    xyz_path = os.path.join(directory, xyz_filename)
    
    if not os.path.exists(xyz_path):
        print(f"WARNING: File not found - {xyz_path}. Skipping.")
        continue
        
    print(f"Processing [{i}/{NUM_ENSEMBLES}]: {xyz_path}")
    
    try:
        box_matrix = get_box_matrix_from_xyz(xyz_path)
        if box_matrix is None:
            print(f"  ERROR: Could not read 'Lattice' information. Skipping.")
            continue

        box_dimensions = convert_matrix_to_dimensions(box_matrix)
        print(f"  Box calculated successfully. Dimensions: {np.round(box_dimensions, 2)}")

        u = mda.Universe(xyz_path)
        
        def set_box(ts):
            ts.dimensions = box_dimensions
            return ts
        u.trajectory.add_transformations(set_box)
        
        atoms = u.select_atoms(ATOM_GROUP_TO_ANALYZE)
        
        if len(atoms) < 2:
            print(f"  WARNING: Fewer than 2 atoms found for '{ATOM_GROUP_TO_ANALYZE}'. Skipping RDF calculation.")
            continue
        
        num_frames = len(u.trajectory)
        start_frame = int(num_frames * (SKIP_PERCENT / 100.0))

        rdf_analysis = rdf.InterRDF(atoms, atoms, 
                                    nbins=NBINS, 
                                    range=RANGE, 
                                    exclusion_block=(1, 1))

        rdf_analysis.run(start=start_frame, step=STEP_SIZE)
        all_rdfs.append(rdf_analysis.results.rdf)
        if distance_bins is None:
            distance_bins = rdf_analysis.results.bins

    except Exception as e:
        print(f"  ERROR: An unexpected problem occurred. - {e}")


if all_rdfs:
    print(f"\nAnalysis completed for a total of {len(all_rdfs)} ensembles.")
    rdf_array = np.array(all_rdfs)
    avg_rdf = np.mean(rdf_array, axis=0)
    std_rdf = np.std(rdf_array, axis=0)
    

    atom_name = ATOM_GROUP_TO_ANALYZE.split()[-1]
    output_file_path = "../result"
    output_filename = f'rdf_{atom_name}-{atom_name}.txt'
    output_file = os.path.join(output_file_path, output_filename)

    df = pd.DataFrame({
        'Distance_A': distance_bins,
        'g(r)_avg': avg_rdf,
        'g(r)_std': std_rdf
    })
    df.to_csv(output_file, index=False, float_format='%.6f', sep='\t') 
    print(f"g(r) results have been saved to {output_file}")


    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(distance_bins, avg_rdf, color='royalblue', lw=2.5, label=f'Ensemble Average g({atom_name}-{atom_name})')
    ax.fill_between(distance_bins, avg_rdf - std_rdf, avg_rdf + std_rdf, color='cornflowerblue', alpha=0.3, label='Standard Deviation')
    ax.set_xlabel('Distance, $r$ (Å)', fontsize=14)
    ax.set_ylabel('$g(r)$', fontsize=14)
    ax.set_title('Ensemble-Averaged Radial Distribution Function', fontsize=16)
    ax.legend(fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_xlim(RANGE)
    ax.set_ylim(bottom=0)
    plt.tight_layout()

    plot_filename = f'rdf_{atom_name}-{atom_name}.png'
    output_fig_file = os.path.join(output_file_path, plot_filename)
    plt.savefig(output_fig_file, dpi=300, bbox_inches='tight')
    print(f"Plot successfully saved as '{plot_filename}'")

    plt.show()

else:
    print("\nNo ensembles were analyzed successfully. Exiting script.")