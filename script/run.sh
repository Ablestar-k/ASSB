#!/bin/bash

GPU_LIST=(0 3 4 5 6)
.
ENSEMBLE_LIST=(1 2 3 4 5)

CORES_PER_JOB=1


for i in "${!GPU_LIST[@]}"
do
    GPU_ID=${GPU_LIST[i]}
    ENSEMBLE_ID=${ENSEMBLE_LIST[i]}

    CPU_START=$((i * CORES_PER_JOB))
    CPU_END=$((CPU_START + CORES_PER_JOB - 1))
    CPU_LIST="$CPU_START-$CPU_END"

    echo "Running Ensemble ${ENSEMBLE_ID}: GPU=${GPU_ID}, CPU Cores=${CPU_LIST}"

    CUDA_VISIBLE_DEVICES=${GPU_ID} \
    OMP_NUM_THREADS=${CORES_PER_JOB} \
    taskset -c ${CPU_LIST} \
    nohup python3 NTOC_Ver1.py ${ENSEMBLE_ID} > NTOC_${ENSEMBLE_ID}.log 2>&1 &

done
