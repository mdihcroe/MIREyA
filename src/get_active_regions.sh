#!/bin/bash

interval=14

OUT_DIR=$1
GENOME_FASTAS=$2
OUT_DIR_ENH_ACTIVE=${OUT_DIR}/enh_active
OUT_DIR_ACTIVE_REGIONS=${OUT_DIR_ENH_ACTIVE}/active_regions
ALL_ENH_ACTIVE=${OUT_DIR}/all_enh_active.bed
OUT_DIR_FASTA=${OUT_DIR_ACTIVE_REGIONS}/fasta
ALL_GENOME_FASTA=${OUT_DIR}/all_chr.fa


[ -d ${OUT_DIR_ACTIVE_REGIONS} ] || mkdir ${OUT_DIR_ACTIVE_REGIONS}
[ -e ${ALL_ENH_ACTIVE} ] && rm ${ALL_ENH_ACTIVE}
[ -d ${OUT_DIR_FASTA} ] || mkdir ${OUT_DIR_FASTA}

for enh_per_mirna in ${OUT_DIR_ENH_ACTIVE}/*enh_active.bed
do
	FILE_NAME=${enh_per_mirna##*/}
	mirna_id="$( cut -d '_' -f 1 <<< "$FILE_NAME" )"
	bedtools intersect -wa \
	-a ${OUT_DIR}/${mirna_id}_seed_pos_in_genome \
	-b ${enh_per_mirna} | sort -k1,1 -k2,2n | uniq | \
	awk -v interval=${interval} '{print $1"\t"($2-interval)"\t"($3+interval)}' \
		> ${OUT_DIR_ACTIVE_REGIONS}/${mirna_id}_active_regions.bed
done

[ -f ${ALL_GENOME_FASTA} ] && echo "Fasta file with whole genome of interest is found. Continue..." || cat ${GENOME_FASTAS}/* > ${ALL_GENOME_FASTA}

for res_per_mirna in ${OUT_DIR_ACTIVE_REGIONS}/*_active_regions.bed
do
	FILE_NAME=${res_per_mirna##*/}
	mirna_id="$( cut -d '_' -f 1 <<< "$FILE_NAME" )"
	bedtools getfasta -fi ${ALL_GENOME_FASTA} -bed ${res_per_mirna} \
					  -fo ${OUT_DIR_FASTA}/${mirna_id}_active_regions.fasta
done