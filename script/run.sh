#!/bin/bash

GPU_LIST=(0 1 7)

ENSEMBLE_LIST=(1 2 3)

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
    nohup python3 NTOC_Ver4_restart.py ${ENSEMBLE_ID} > NTOC_ver4_${ENSEMBLE_ID}_restart.log 2>&1 &

done
