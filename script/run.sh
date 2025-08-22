#!/bin/bash

CUDA_VISIBLE_DEVICES=1

for i in {1..5}
do
    echo "Running ensemble $i ..."
    nohup python3 NTOC_Ver1.py $i > NTOC_$i.log 2>&1
    echo "Ensemble $i finished. Waiting for 10 seconds before starting the next one..."
    sleep 10
done

echo "All ensembles have been submitted."