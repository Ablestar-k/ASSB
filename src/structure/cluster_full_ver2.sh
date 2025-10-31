#!/bin/bash

PROGRAM="./cluster_full_ver2.x"

DUMP_BASE_DIR="../../dump/dump_NTOC_ver4_"
RESULT_BASE_DIR="../../result/structure"

LOG_FILE="${RESULT_BASE_DIR}/run_ensembles_ver2.log"

ENSEMBLE_START=1
ENSEMBLE_END=5

mkdir -p "$RESULT_BASE_DIR"
echo "Output directory: $RESULT_BASE_DIR (Checked/Created)."

(
    echo "Starting processing for $ENSEMBLE_END ensembles..."
    echo "Log file: $LOG_FILE"

    for (( i=$ENSEMBLE_START; i<=$ENSEMBLE_END; i++ ))
    do
        echo "------------------------------------------------"
        echo ">>> Starting processing for ensemble $i..."

        INPUT_FILE="${DUMP_BASE_DIR}${i}/${i}_product_unwrapped.xyz"

        OUT_CLUSTER_SIZE="${RESULT_BASE_DIR}/${i}_cluster_size_ver2.dat"
        OUT_NA_COORD="${RESULT_BASE_DIR}/${i}_na_coord_ver2.dat"
        OUT_CLUSTER_DIST="${RESULT_BASE_DIR}/${i}_cluster_dist_ver2.dat"
        OUT_VHF_POLY="${RESULT_BASE_DIR}/${i}_vhf_poly_ver2.dat"
        OUT_VHF_ISO="${RESULT_BASE_DIR}/${i}_vhf_iso_ver2.dat"
        OUT_MSD_TA_POLY="${RESULT_BASE_DIR}/${i}_msd_ta_poly_ver2.dat"
        OUT_MSD_TA_ISO="${RESULT_BASE_DIR}/${i}_msd_ta_iso_ver2.dat"
        OUT_NA_O_COORD="${RESULT_BASE_DIR}/${i}_na_o_coord_ver2.dat"
        OUT_VHF_NA_NBO="${RESULT_BASE_DIR}/${i}_vhf_na_nbo_ver2.dat"
        OUT_VHF_NA_BO="${RESULT_BASE_DIR}/${i}_vhf_na_bo_ver2.dat"


        if [ ! -f "$INPUT_FILE" ]; then
            echo "Warning: Input file '$INPUT_FILE' for ensemble $i not found."
            echo "Skipping this ensemble."
            continue
        fi

        echo "Executing C program: $PROGRAM"
        echo "  Input: $INPUT_FILE"
        echo "  Output: ${RESULT_BASE_DIR}/${i}_*_ver2.dat" 

        "$PROGRAM" "$INPUT_FILE" \
                   "$OUT_CLUSTER_SIZE" "$OUT_NA_COORD" "$OUT_CLUSTER_DIST" \
                   "$OUT_VHF_POLY" "$OUT_VHF_ISO" \
                   "$OUT_MSD_TA_POLY" "$OUT_MSD_TA_ISO" \
                   "$OUT_NA_O_COORD" "$OUT_VHF_NA_NBO" "$OUT_VHF_NA_BO"

        echo ">>> Finished processing ensemble $i."
    done

    echo "------------------------------------------------"
    echo "All ensemble processing is complete. :D"

) > "$LOG_FILE" 2>&1 &

echo "Ensemble processing (1-$ENSEMBLE_END) started in the background."
echo "PID: $!"
echo "You can monitor the progress by typing the following command in your terminal:"
echo "tail -f $LOG_FILE"