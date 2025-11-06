#!/bin/bash

PROGRAM="./cluster_full_ver5.x" 
ENSEMBLE_START=1
ENSEMBLE_END=5

DUMP_BASE_DIR="../../dump/dump_NTOC_ver4_"
RESULT_BASE_DIR="../../result/structure"

LOG_FILE="${RESULT_BASE_DIR}/run_ensembles_ver5.log" 

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

        OUT_CLUSTER_SIZE="${RESULT_BASE_DIR}/${i}_cluster_size_ver5.dat"  # 2
        OUT_NA_COORD="${RESULT_BASE_DIR}/${i}_na_coord_ver5.dat"      # 3
        OUT_CLUSTER_DIST="${RESULT_BASE_DIR}/${i}_cluster_dist_ver5.dat"  # 4
        OUT_VHF_POLY="${RESULT_BASE_DIR}/${i}_vhf_poly_ver5.dat"      # 5
        OUT_VHF_ISO="${RESULT_BASE_DIR}/${i}_vhf_iso_ver5.dat"       # 6
        OUT_MSD_TA_POLY="${RESULT_BASE_DIR}/${i}_msd_ta_poly_ver5.dat" # 7
        OUT_MSD_TA_ISO="${RESULT_BASE_DIR}/${i}_msd_ta_iso_ver5.dat"  # 8
        
        OUT_NA_ANION_RATIO="${RESULT_BASE_DIR}/${i}_na_anion_ratio_ver5.dat" # 9
        OUT_VHF_NA_NBO_UO="${RESULT_BASE_DIR}/${i}_vhf_na_nbo_uo_ver5.dat"  # 10
        OUT_VHF_NA_BO="${RESULT_BASE_DIR}/${i}_vhf_na_bo_ver5.dat"        # 11
        OUT_VHF_NA_CL="${RESULT_BASE_DIR}/${i}_vhf_na_cl_ver5.dat"        # 12 

        OUT_TA_O_COORD="${RESULT_BASE_DIR}/${i}_ta_o_coord_ver5.dat"   # 13
        OUT_TA_BRIDGE="${RESULT_BASE_DIR}/${i}_ta_bridge_ver5.dat"     # 14
        OUT_O_TYPE_DIST="${RESULT_BASE_DIR}/${i}_o_type_dist_ver5.dat"  # 15

        if [ ! -f "$INPUT_FILE" ]; then
            echo "Warning: Input file '$INPUT_FILE' for ensemble $i not found."
            echo "Skipping this ensemble."
            continue
        fi

        echo "Executing C program: $PROGRAM"
        echo "  Input: $INPUT_FILE"
        echo "  Output: ${RESULT_BASE_DIR}/${i}_*_ver5.dat" 

        "$PROGRAM" "$INPUT_FILE" \
            "$OUT_CLUSTER_SIZE" "$OUT_NA_COORD" "$OUT_CLUSTER_DIST" \
            "$OUT_VHF_POLY" "$OUT_VHF_ISO" \
            "$OUT_MSD_TA_POLY" "$OUT_MSD_TA_ISO" \
            "$OUT_NA_ANION_RATIO" \
            "$OUT_VHF_NA_NBO_UO" "$OUT_VHF_NA_BO" "$OUT_VHF_NA_CL" \
            "$OUT_TA_O_COORD" "$OUT_TA_BRIDGE" "$OUT_O_TYPE_DIST"

        echo ">>> Finished processing ensemble $i."
    done

    echo "------------------------------------------------"
    echo "All ensemble processing is complete. :D"

) > "$LOG_FILE" 2>&1 &

echo "Ensemble processing (1-$ENSEMBLE_END) started in the background."
echo "PID: $!"
echo "You can monitor the progress by typing the following command in your terminal:"
echo "tail -f $LOG_FILE"