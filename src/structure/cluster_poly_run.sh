#!/bin/bash

for i in {1..5}
do
	echo "Processing ensemble ver4_$i..."

	nohup ./cluster_ver3.x /SSD2/sejung/NTOC/dump/dump_NTOC_ver4_$i/${i}_product.xyz /SSD2/sejung/NTOC/result/structure/${i}_cluster_polymer.dat /SSD2/sejung/NTOC/result/structure/${i}_Na_cluster_polymer.dat > ${i}_clster_polymer.log 2>&1 &
done

echo "Clustering Analysis Finished"
