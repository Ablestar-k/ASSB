#!/bin/bash

for i in {1..5}
do
	echo "Processing ensemble ver4_$i..."

	echo "  Converting product_nvt.traj..."
	rm -f ../dump/dump_NTOC_ver4_$i/${i}_product.xyz
	ase convert ../dump/dump_NTOC_ver4_$i/NTOC_ver4_${i}_product_nvt.traj ../dump/dump_NTOC_ver4_$i/${i}_product.xyz
done

echo "Finished all conversions."
