#!/bin/bash

# --- 설정 ---
PROGRAM="./cluster_full_ver4.x" # [수정] 새 실행 파일 이름
ENSEMBLE_START=1
ENSEMBLE_END=5

DUMP_BASE_DIR="../../dump/dump_NTOC_ver4_"
RESULT_BASE_DIR="../../result/structure"

LOG_FILE="${RESULT_BASE_DIR}/run_ensembles_ver4.log" # [수정] 새 로그 파일

mkdir -p "$RESULT_BASE_DIR"
echo "Output directory: $RESULT_BASE_DIR (Checked/Created)."

(
    # ... (로그 메시지)

    for (( i=$ENSEMBLE_START; i<=$ENSEMBLE_END; i++ ))
    do
        # ... (앙상블 시작 메시지)
        INPUT_FILE="${DUMP_BASE_DIR}${i}/${i}_product_unwrapped.xyz"

        # [수정] 모든 ver3 -> ver4
        OUT_CLUSTER_SIZE="${RESULT_BASE_DIR}/${i}_cluster_size_ver4.dat"
        OUT_NA_COORD="${RESULT_BASE_DIR}/${i}_na_coord_ver4.dat"
        OUT_CLUSTER_DIST="${RESULT_BASE_DIR}/${i}_cluster_dist_ver4.dat"
        OUT_VHF_POLY="${RESULT_BASE_DIR}/${i}_vhf_poly_ver4.dat"
        OUT_VHF_ISO="${RESULT_BASE_DIR}/${i}_vhf_iso_ver4.dat"
        OUT_MSD_TA_POLY="${RESULT_BASE_DIR}/${i}_msd_ta_poly_ver4.dat"
        OUT_MSD_TA_ISO="${RESULT_BASE_DIR}/${i}_msd_ta_iso_ver4.dat"
        OUT_NA_O_COORD="${RESULT_BASE_DIR}/${i}_na_o_coord_ver4.dat"
        OUT_VHF_NA_NBO="${RESULT_BASE_DIR}/${i}_vhf_na_nbo_ver4.dat"
        OUT_VHF_NA_BO="${RESULT_BASE_DIR}/${i}_vhf_na_bo_ver4.dat"
        OUT_TA_O_COORD="${RESULT_BASE_DIR}/${i}_ta_o_coord_ver4.dat"
        OUT_TA_BRIDGE="${RESULT_BASE_DIR}/${i}_ta_bridge_ver4.dat"
        
        # [신규 추가] 새 출력 파일
        OUT_O_TYPE_DIST="${RESULT_BASE_DIR}/${i}_o_type_dist_ver4.dat"

        # ... (파일 확인)

        echo "Executing C program: $PROGRAM"
        echo "  Input: $INPUT_FILE"
        echo "  Output: ${RESULT_BASE_DIR}/${i}_*_ver4.dat" 

        # [수정] 마지막에 $OUT_O_TYPE_DIST 인자 추가 (총 15개 인자)
        "$PROGRAM" "$INPUT_FILE" \
            "$OUT_CLUSTER_SIZE" "$OUT_NA_COORD" "$OUT_CLUSTER_DIST" \
            "$OUT_VHF_POLY" "$OUT_VHF_ISO" \
            "$OUT_MSD_TA_POLY" "$OUT_MSD_TA_ISO" \
            "$OUT_NA_O_COORD" "$OUT_VHF_NA_NBO" "$OUT_VHF_NA_BO" \
            "$OUT_TA_O_COORD" "$OUT_TA_BRIDGE" \
            "$OUT_O_TYPE_DIST" # <-- 신규 인자

        echo ">>> Finished processing ensemble $i."
    done


) > "$LOG_FILE" 2>&1 &