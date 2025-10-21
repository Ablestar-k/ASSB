import numpy as np
from ase.io import read, write
import os

BASE_DIR_NAME = "../dump/dump_NTOC_ver4_{}"
BASE_FILE_NAME = "NTOC_ver4_{}_product_nvt.traj"
OUTPUT_FILE_NAME = "{}_product_unwrapped.xyz"  
NUM_ENSEMBLES = 5

def unwrap_trajectory_memory_efficient(trajectory_file, output_file):

    print(f"Reading trajectory from: {trajectory_file}")
    traj = read(trajectory_file, index=':')

    if not traj:
        print("Trajectory is empty.")
        return

    print("Unwrapping frames and writing to file...")
    first_frame = traj[0].copy()
    write(output_file, first_frame, format='extxyz')

    prev_unwrapped_pos = first_frame.get_positions()
    images = np.zeros_like(prev_unwrapped_pos, dtype=int)

    for i in range(1, len(traj)):
        current_atoms = traj[i].copy()

        cell = current_atoms.get_cell()
        inv_cell = np.linalg.inv(cell)

        scaled_pos = np.dot(current_atoms.get_positions(), inv_cell)
        scaled_prev_pos = np.dot(prev_unwrapped_pos, inv_cell)
        
        jump = np.rint(scaled_prev_pos - scaled_pos)
        images += jump.astype(int)
        
        new_unwrapped_pos = current_atoms.get_positions() + np.dot(images, cell)
        current_atoms.set_positions(new_unwrapped_pos)

        write(output_file, current_atoms, format='extxyz', append=True)
        prev_unwrapped_pos = new_unwrapped_pos

    print(f"Successfully wrote unwrapped trajectory to: {output_file}")

print("--- Starting trajectory unwrapping for all ensembles ---")
for i in range(1, NUM_ENSEMBLES + 1):
    dir_path = BASE_DIR_NAME.format(i)
    input_filename = BASE_FILE_NAME.format(i)
    output_filename = OUTPUT_FILE_NAME.format(i) 

    input_file = os.path.join(dir_path, input_filename)
    output_file = os.path.join(dir_path, output_filename)

    if os.path.exists(input_file):
        print(f"\n--- Processing Ensemble #{i} ---")
        unwrap_trajectory_memory_efficient(input_file, output_file)
    else:
        print(f"\nWARNING: Input file not found, skipping: {input_file}")

print("\n--- All processing finished. ---")