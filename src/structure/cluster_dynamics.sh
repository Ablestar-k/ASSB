#!/bin/bash

PROGRAM="./cluster_dynamics.x"

DUMP_BASE_DIR="../../dump/dump_NTOC_ver4_"
RESULT_BASE_DIR="../../result/dynamics"

LOG_FILE="${RESULT_BASE_DIR}/run_dynamics_analysis.log"

ENSEMBLE_START=1
ENSEMBLE_END=5

MAX_CORES=5 

mkdir -p "$RESULT_BASE_DIR"
echo "Output directory: $RESULT_BASE_DIR (Checked/Created)."

(
    echo "Starting processing for $ENSEMBLE_END ensembles..."
    echo "Log file: $LOG_FILE"
    echo "Using up to $MAX_CORES cores concurrently."

    for (( i=$ENSEMBLE_START; i<=$ENSEMBLE_END; i++ ))
    do
        (
            echo "------------------------------------------------"
            echo ">>> Starting processing for ensemble $i (PID: $$)..."

            INPUT_FILE="${DUMP_BASE_DIR}${i}/${i}_product_unwrapped.xyz"

            OUT_TCF_POLY="${RESULT_BASE_DIR}/${i}_tcf_poly.dat"
            OUT_TCF_ISO="${RESULT_BASE_DIR}/${i}_tcf_iso.dat"
            OUT_TCF_SIZE="${RESULT_BASE_DIR}/${i}_tcf_size_corr.dat"
            OUT_LIFE_ISO_POLY="${RESULT_BASE_DIR}/${i}_life_iso_poly.dat"
            OUT_LIFE_POLY_ISO="${RESULT_BASE_DIR}/${i}_life_poly_iso.dat"

            if [ ! -f "$INPUT_FILE" ]; then
                echo "Warning: Input file '$INPUT_FILE' for ensemble $i not found."
                echo "Skipping this ensemble."
                exit 1
            fi

            echo "Executing C program: $PROGRAM (Ensemble $i)"
            echo "  Input: $INPUT_FILE"

            "$PROGRAM" "$INPUT_FILE" \
                       "$OUT_TCF_POLY" "$OUT_TCF_ISO" \
                       "$OUT_TCF_SIZE" "$OUT_LIFE_ISO_POLY" "$OUT_LIFE_POLY_ISO"

            echo ">>> Finished processing ensemble $i."

        ) &

        if (( $(jobs -p | wc -l) >= $MAX_CORES )); then
            wait -n
        fi

    done

    wait

    echo "------------------------------------------------"
    echo "All ensemble processing is complete. :D"

) > "$LOG_FILE" 2>&1 &

echo "Dynamics processing (1-$ENSEMBLE_END) started in the background (Main PID: $!)."
echo "Jobs will run in parallel (max $MAX_CORES cores)."
echo "Monitor progress: tail -f $LOG_FILE"
