#!/bin/bash

OUT_DATA_DIR=$1

ALIGNMENTS_NEEDLE_DIR=${OUT_DATA_DIR}"/alignments_needle"
ACTIVE_ENHANCERS=${OUT_DATA_DIR}"/enh_active"
ALL_ENH_ACTIVE=${ACTIVE_ENHANCERS}/all_enh_active.bed

[ -e ${ALL_ENH_ACTIVE} ] && rm ${ALL_ENH_ACTIVE}

cat ${ACTIVE_ENHANCERS}/*_enh_active.bed > ${ALL_ENH_ACTIVE}

bedtools intersect -wa -a ${ALL_ENH_ACTIVE} \
				   -wb -b ${ALIGNMENTS_NEEDLE_DIR}"/alignment_best_regions.bed" | \
				   sort -k1,1 -k2,2n | uniq  \
				   > ${ALIGNMENTS_NEEDLE_DIR}/alignment_best_enh.bed