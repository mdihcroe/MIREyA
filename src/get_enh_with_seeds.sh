#!/bin/bash

DB_PATH=$1
FORWARD_SEEDS=$2
REVERSE_COMPL_SEEDS=$3
OUT_DATA_DIR=$4
ENHANCERS_BED=$5

[ -d ${OUT_DATA_DIR} ] || mkdir ${OUT_DATA_DIR}

mirnas=$(sort -k1,1 -u ${FORWARD_SEEDS} | cut -f1)
for mirna in ${mirnas}
do
	echo ${mirna}
	rm -f ${OUT_DATA_DIR}/${mirna}_seed_pos_in_genome
	seeds1=$(grep -P "^"${mirna}'\t' ${REVERSE_COMPL_SEEDS} | cut -f2)
	seeds2=$(grep -P "^"${mirna}'\t' ${FORWARD_SEEDS} | cut -f2)
	seeds=$(echo ${seeds1}; echo ${seeds2})
	for seed in ${seeds}
	do
		echo ${seed}
		for FILE_PATH in ${DB_PATH}/chr*.fa
		do
			echo ${FILE_PATH}
			FILE_NAME=${FILE_PATH##*/}
			IFS='.' read -r chr FILE_NAME <<< "$FILE_NAME"

			tail ${FILE_PATH} -n +2 | tr -d '\n' | \
			grep -aob ${seed} | grep -oE '[0-9]+' > ${OUT_DATA_DIR}/${mirna}_tmp.start_positions
			awk -v prefix="$chr" -v postfix=" ${#seed}" '{print prefix"\t"$0"\t"$0+postfix}' \
			${OUT_DATA_DIR}/${mirna}_tmp.start_positions >> ${OUT_DATA_DIR}/${mirna}_seed_pos_in_genome
		done
	done
	rm -f ${OUT_DATA_DIR}/${mirna}_tmp.start_positions
	bedtools intersect -wa -a ${ENHANCERS_BED} -b ${OUT_DATA_DIR}/${mirna}_seed_pos_in_genome | sort -k1,1 -k2,2n | uniq \
            > ${OUT_DATA_DIR}/${mirna}_enh_with_seeds.bed
done
