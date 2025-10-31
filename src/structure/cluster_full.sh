#!/bin/bash
PROGRAM="./cluster_full.x"

DUMP_BASE_DIR="../../dump/dump_NTOC_ver4_"
RESULT_BASE_DIR="../../result/structure"

ENSEMBLE_START=1
ENSEMBLE_END=5

mkdir -p "$RESULT_BASE_DIR"
echo "출력 디렉토리: $RESULT_BASE_DIR (존재 여부 확인/생성 완료)"


echo "총 $ENSEMBLE_END 개의 앙상블 처리를 시작합니다..."

for (( i=$ENSEMBLE_START; i<=$ENSEMBLE_END; i++ ))
do
    echo "------------------------------------------------"
    echo ">>> 앙상블 $i 처리 시작..."

    INPUT_FILE="${DUMP_BASE_DIR}${i}/${i}_product_unwrapped.xyz"

    OUT_CLUSTER_SIZE="${RESULT_BASE_DIR}/${i}_cluster_size.dat"
    OUT_NA_COORD="${RESULT_BASE_DIR}/${i}_na_coord.dat"
    OUT_CLUSTER_DIST="${RESULT_BASE_DIR}/${i}_cluster_dist.dat"
    OUT_VHF_POLY="${RESULT_BASE_DIR}/${i}_vhf_poly.dat"
    OUT_VHF_ISO="${RESULT_BASE_DIR}/${i}_vhf_iso.dat"

    if [ ! -f "$INPUT_FILE" ]; then
        echo "경고: 앙상블 $i 의 입력 파일 '$INPUT_FILE'을(를) 찾을 수 없습니다."
        echo "이 앙상블을 건너뜁니다."
        continue 
    fi

    echo "C 프로그램 실행: $PROGRAM"
    echo "  입력: $INPUT_FILE"
    echo "  출력: ${RESULT_BASE_DIR}/${i}_*.dat"
    
    "$PROGRAM" "$INPUT_FILE" "$OUT_CLUSTER_SIZE" "$OUT_NA_COORD" \
               "$OUT_CLUSTER_DIST" "$OUT_VHF_POLY" "$OUT_VHF_ISO"

    echo ">>> 앙상블 $i 처리 완료."
done

echo "------------------------------------------------"
echo "모든 앙상블 처리가 완료되었습니다. :D"
