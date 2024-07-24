#!/bin/bash

##-------------------------------------------------------------------------
## Configs + packages

source ../config/config.sh

module purge
module load plink/2.0

###########################################################################
###########################################################################
###                                                                     ###
###                           DATA CREATE                               ###
###                                                                     ###
###########################################################################
###########################################################################

### Files + prefixes
DATA_DIR=${OUT_DIR}/covid_afs_data
POP_LIST=${OUT_DIR}/combined_data/aggCOVID_unrelated_pops.fam

### COVID
CHUNK_LIST=../config/aggCOVIDv5.file_list
TIERED_VARIANTS_TAB=${OUT_DIR}/labkey_data/tiering_sorted_variants.tsv

### Make key of tiered variants based on PLINK encoded variant IDs
sed 's/,/_/g' $TIERED_VARIANTS_TAB | \
awk '{print $1,":",$2,"_",$3}' | \
sed 's/ //g' > ${TMP_DIR}/tiering_sorted_variants_key.txt

## Run per file chunk
while read chunk
do
  CHROM=$(echo $chunk | awk '{print $1}')
  START=$(echo $chunk | awk '{print $2}')
  END=$(echo $chunk | awk '{print $3}')
  ANNOT_VCF=$(echo $chunk | awk '{print $4}')
  GT_PGEN=$(echo $chunk | awk '{print $5}')

  ## File chunk output prefix
  OUT_PREFIX=${CHROM}_${START}_${END}

  POP_ARRAY=( $( awk '{print $7}' $POP_LIST | sort | uniq ) )
  COUNTER=1
  for POP in "${POP_ARRAY[@]}"
  do
    grep $POP $POP_LIST | awk '{print $1,$2,$3,$4,$5,$6}' \
> ${TMP_DIR}/${OUT_PREFIX}_${POP}.fam

    ## Get frequencies for this specific group
    plink2 --pfile $GT_PGEN \
--freq counts \
--keep ${TMP_DIR}/${OUT_PREFIX}_${POP}.fam \
--out ${TMP_DIR}/${OUT_PREFIX}_${POP}

    if [[ "$COUNTER" == 1 ]]
    then
      awk '{print $1,$2,$3,$4}' ${TMP_DIR}/${OUT_PREFIX}_${POP}.acount \
> ${TMP_DIR}/${OUT_PREFIX}.acount
    fi

    sed "s/ALT_CTS/${POP}_ALT_CTS/g" ${TMP_DIR}/${OUT_PREFIX}_${POP}.acount | \
    sed "s/OBS_CT/${POP}_OBS_CTS/g" | \
    awk '{print $5,$6}' | paste ${TMP_DIR}/${OUT_PREFIX}.acount - > ${TMP_DIR}/${OUT_PREFIX}_tmp && \
mv ${TMP_DIR}/${OUT_PREFIX}_tmp ${TMP_DIR}/${OUT_PREFIX}.acount
  
    rm ${TMP_DIR}/${OUT_PREFIX}_${POP}*
    
    COUNTER=$(( $COUNTER + 1 ))
  done

  ## Get file header (replace each time)
  head -1 ${TMP_DIR}/${OUT_PREFIX}.acount | \
awk 'OFS="\t" {print}0' | \
sed 's/ /\t/g' > ${DATA_DIR}/aggCOVID_full.acount

done < $CHUNK_LIST

### Combine chunks
for acount_chunk in ${TMP_DIR}/${OUT_PREFIX}*acount
do
  tail -n+2 $acount_chunk | \
awk 'OFS="\t" {print}0' | \
sed 's/ /\t/g' >> ${DATA_DIR}/aggCOVID_full.acount
done

## Compress with gzip
gzip ${DATA_DIR}/aggCOVID_full.acount

## Get only tiered variants...
zcat ${DATA_DIR}/aggCOVID_full.acount.gz | \
awk -F"\t" 'FNR==NR {hash[$1]; next} $2 in hash' ${TMP_DIR}/tiering_sorted_variants_key.txt - | \
gzip -f > ${DATA_DIR}/aggCOVID_tiering.acount.gz

###########################################################################
###########################################################################
###                                                                     ###
###                         END OF SCRIPT                               ###
###                                                                     ###
###########################################################################
###########################################################################
