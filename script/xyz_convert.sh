#!/bin/bash

for i in {1..5}
do
	echo "Processing ensemble ver1_$i..."

	echo "  Converting heating.traj..."
	ase convert ../dump/dump_NTOC_ver1_$i/NTOC_ver1_heating.traj ${i}_heating.xyz
	
	echo "  Converting heating_eq_nvt.traj..."
	ase convert ../dump/dump_NTOC_ver1_$i/NTOC_ver1_heating_eq_nvt.traj ${i}_heating_eq.xyz
	
	echo "  Converting quench.traj..."
	ase convert ../dump/dump_NTOC_ver1_$i/NTOC_ver1_quench.traj ${i}_quench.xyz
	
	echo "  Converting quench_eq_nvt.traj..."
	ase convert ../dump/dump_NTOC_ver1_$i/NTOC_ver1_quench_eq_nvt.traj ${i}_quench_eq.xyz
	
	echo "  Converting product_nve.traj..."
	ase convert ../dump/dump_NTOC_ver1_$i/NTOC_ver1_product_nve.traj ${i}_product.xyz
done

echo "Finished all conversions."
