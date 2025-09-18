import MDAnalysis as mda
import MDAnalysis.analysis.msd as msd
from MDAnalysis.transformations import unwrap, nojump
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os


BASE_DIR_NAME = "../dump/dump_NTOC_ver3_{}"
NUM_ENSEMBLES = 5
ATOM_GROUP_TO_ANALYZE = "name Na"
SIMULATION_DT = 2.0

all_msd = []
print("Starting MSD analysis...")

for i in range(1, NUM_ENSEMBLES + 1):
    directory = BASE_DIR_NAME.format(i)

    #xyz_filename = f"{i}_product.xyz"
    xyz_filename = f"{i}_quench_eq.xyz"
    xyz_path = os.path.join(directory, xyz_filename)

    if not os.path.exists(xyz_path):
        print(f"WARNING: File not found - {xyz_path}. Skipping.")
        continue

    print(f"Processing [{i}/{NUM_ENSEMBLES}]: {xyz_path}")

    
    try:
        u = mda.Universe(xyz_path)
        ag = u.select_atoms(ATOM_GROUP_TO_ANALYZE)

        if len(ag) == 0:
            print(f"WARNING: No atoms matched in {xyz_path}. Skipping.")
            continue


        transform = mda.transformations.nojump.Nojump(u.atoms) 
        u.trajectory.add_transformations(transform)

        msd_analysis = msd.EinsteinMSD(ag, msd_type='xyz', fft=True)
        msd_analysis.run()

        all_msd.append(msd_analysis.results.timeseries)

    except Exception as e:
        print(f"ERROR processing {xyz_path}: {e}")

if all_msd:
    print(f"\nMSD analysis completed for {len(all_msd)} ensembles.")
   
    max_len = max(len(m) for m in all_msd)
    msd_padded = [np.pad(m, (0, max_len - len(m)), 'constant', constant_values=np.nan) for m in all_msd]
    msd_array = np.array(msd_padded)
    

    msd_mean = np.nanmean(msd_array, axis=0)
    msd_std = np.nanstd(msd_array, axis=0)
    
    num_lags = len(msd_mean)
    time_lags_fs = np.arange(num_lags) * SIMULATION_DT
    time_lags_ps = time_lags_fs / 1000.0

    output_path = "../result"
    os.makedirs(output_path, exist_ok=True) 

    atom_name = ATOM_GROUP_TO_ANALYZE.split()[-1]
    output_data_filename = f'{atom_name}_MSD_quench_ensemble_average.dat'
    output_data_file = os.path.join(output_path, output_data_filename)

    data_to_save = {
        'Time(ps)': time_lags_ps,
        'MSD_avg(A^2)': msd_mean,
        'MSD_std(A^2)': msd_std
    }
    df_to_save = pd.DataFrame(data_to_save)
    df_to_save.to_csv(output_data_file, sep='\t', index=False, float_format='%.6f')
    print(f"\nMSD data successfully saved to '{output_data_file}'")

    # --- Plotting ---
    plt.figure(figsize=(10, 7))
    plt.plot(time_lags_ps, msd_mean, label='Ensemble Average MSD', color='royalblue', linewidth=2)
    plt.fill_between(time_lags_ps, msd_mean - msd_std, msd_mean + msd_std, color='cornflowerblue', alpha=0.3, label='Standard Deviation')
    
    plt.xlabel('Time Lag (ps)', fontsize=14)
    plt.ylabel('MSD ($Å^2$)', fontsize=14)
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Ensemble-Averaged Mean Squared Displacement (MSD)', fontsize=16, fontweight='bold')
    plt.grid(True, which="both", linestyle='--', alpha=0.6)
    plt.legend(fontsize=12)
    plt.tick_params(axis='both', which='major', labelsize=12)
    

    output_fig_filename = f'MSD_{atom_name}_final.png'
    output_fig_file = os.path.join(output_path, output_fig_filename)
    plt.savefig(output_fig_file, dpi=300, bbox_inches='tight')
    print(f"Plot successfully saved as '{output_fig_file}'")
    
    plt.show()

else:
    print("\nNo ensembles were analyzed successfully. Exiting script.")