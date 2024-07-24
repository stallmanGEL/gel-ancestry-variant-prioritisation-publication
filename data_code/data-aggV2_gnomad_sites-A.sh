#!/bin/bash

##-------------------------------------------------------------------------
## Configs + packages

source ../config/config.cfg

module purge
module load plink/2.0

###########################################################################
###########################################################################
###                                                                     ###
###                           DATA CREATE                               ###
###                                                                     ###
###########################################################################
###########################################################################

### AGGV2
CHUNK_LIST=../config/aggV2.file_list
DATA_DIR=${OUT_DIR}/aggV2_gnomad_sites_data

### Reformat ancestry alleles TSV (from gnomAD) to site IDs
sed 's/\t\["/_/g' ${GNOMAD_ANC_ALLELES} | sed 's/","/_/g' | sed 's/"\]//g' | tail -n+2 > ${TMP_DIR}/gnomad_alleles.ids

## Run per file chunk
while read chunk
do
  CHROM=$(echo $chunk | awk '{print $1}')
  START=$(echo $chunk | awk '{print $2}')
  END=$(echo $chunk | awk '{print $3}')
  GT_PGEN=$(echo $chunk | awk '{print $5}')

  echo "${TMP_DIR}/aggV2_chr${CHROM}_${START}_${END}_anc_alleles" >> ${TMP_DIR}/aggV2_anc_alleles_merge.list
  plink2 --pfile ${GT_PGEN} \
    --extract ${TMP_DIR}/gnomad_alleles.ids \
    --geno 0.05 \
    --make-pgen \
    --memory 10000 \
    --threads 1 \
    --out ${TMP_DIR}/aggV2_chr${CHROM}_${START}_${END}_anc_alleles
done

### Combine chunks and write to the aggV2_gnomad_sites data directory
plink2 --pmerge-list ${TMP_DIR}/aggV2_anc_alleles_merge.list \
  --make-bed \
  --memory 10000 \
  --threads 1 \
  --out ${DATA_DIR}/aggV2_all_chr_gnomad_anc