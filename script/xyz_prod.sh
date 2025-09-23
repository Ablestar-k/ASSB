#!/bin/bash

for i in {1..5}
do
	echo "Processing ensemble ver3_$i..."
	
	echo "  Converting quench_eq_nvt.traj..."
	ase convert ../dump/dump_NTOC_ver3_$i/NTOC_ver3_${i}_pre_product.traj ../dump/dump_NTOC_ver3_$i/${i}_pre_prodcut.xyz
	
	echo "  Converting product_nve.traj..."
	ase convert ../dump/dump_NTOC_ver3_$i/NTOC_ver3_${i}_product_nve.traj ../dump/dump_NTOC_ver3_$i/${i}_product.xyz
done

echo "Finished all conversions."
