#!/bin/bash

##-------------------------------------------------------------------------
## Configs + packages

source ../config/config.sh

module purge
module load bcftools/1.16
module load tabix/1.18

###########################################################################
###########################################################################
###                                                                     ###
###                        DATA CREATE VEP                              ###
###                                                                     ###
###########################################################################
###########################################################################

### Variables
DATA_DIR=${OUT_DIR}/deleterious_predictions_data
TIERED_VARIANTS_TAB=${OUT_DIR}/labkey_data/tiering_sorted_variants.tsv

### Run per chunk
CHUNK_LIST=../config/aggV2.file_list

### Index tiered variants tab file.
bgzip $TIERED_VARIANTS_TAB && tabix -s1 -b2 -e2 ${TIERED_VARIANTS_TAB}.gz

while read chunk
do
  CHROM=$(echo $chunk | awk '{print $1}')
  START=$(echo $chunk | awk '{print $2}')
  END=$(echo $chunk | awk '{print $3}')
  ANNOT_VCF=$(echo $chunk | awk '{print $4}')

  ### Run
  bcftools +split-vep \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%CADD_RAW\t%CADD_PHRED\t%PrimateAI_pred\n' \
    -d \
    -s primary \
    -T ${TIERED_VARIANTS_TAB}.gz \
    -A tab ${ANNOT_VCF} > ${TMP_DIR}/${CHROM}_${START}_${END}.vep_tmp.tsv
  
  ## Get file header (replace each time)
  head -1 ${TMP_DIR}/${CHROM}_${START}_${END}.vep_tmp.tsv > ${DATA_DIR}/vep_preds_data.tsv
done < $CHUNK_LIST

### Combine chunks
for vep_chunk in ${TMP_DIR}/*vep_tmp.tsv
do
  tail -n+2 $vep_chunk >> ${DATA_DIR}/cadd_scores_tiered.tsv
done

## Compress with gzip
gzip -f ${DATA_DIR}/cadd_scores_tiered.tsv

###########################################################################
###########################################################################
###                                                                     ###
###                        DATA CREATE PAI3D                            ###
###                                                                     ###
###########################################################################
###########################################################################

### Make key of tiered variants based on PLINK encoded variant IDs
sed 's/,/_/g' $TIERED_VARIANTS_TAB | \
awk '{print $1,":",$2,"_",$3}' | \
sed 's/ //g' > ${TMP_DIR}/tiering_sorted_variants_key.txt

zcat ${AM_SCORES} | \
sed '/^#/d' | \
awk 'OFS="\t" { $(NF+1) = $1":"$2"_"$3"_"$4 }1' | \
awk -F"\t" 'FNR==NR {hash[$1]; next} $(NF) in hash' ${TMP_DIR}/tiering_sorted_variants_key.txt - | \
awk '{print $1,$2,$3,$4,$9,$10}' > ${TMP_DIR}/tmp_am_scores

printf "CHROM\tPOS\tREF\tALT\tam_pathogenicity\tam_class\n" > ${TMP_DIR}/tmp_am_header
cat tmp_am_header tmp_am_scores | gzip -f > ${DATA_DIR}/alpha_missense_scores_tiered.tsv.gz

rm tmp_am_header tmp_am_scores

###########################################################################
###########################################################################
###                                                                     ###
###                        DATA CREATE PAI3D                            ###
###                                                                     ###
###########################################################################
###########################################################################

#####################################################################
###### THESE FILES CANNOT BE MADE AVAILABLE TO ALL IN THE RE ########
#####################################################################

### DO YOU HAVE PRIMATEAI-3D ACCESS?
PAI3D_ACCESS=$1

if [[ "$PAI3D_ACCESS" == "PAI3D_ACCESS" ]]
then
  echo "Selected that you have access to PrimateAI-3D scores..."
  echo ${PAI3D_SCORES}
  zcat ${PAI3D_SCORES} | \
  awk -F"," 'OFS="\t" { $(NF+1) = $1":"$2"_"$3"_"$4 }1' | \
  awk -F"\t" 'FNR==NR {hash[$1]; next} $(NF) in hash' ${TMP_DIR}/tiering_sorted_variants_key.txt - | \
  awk '{print $1,$2,$3,$4,$8}' > ${TMP_DIR}/tmp_pai_scores

  printf "chromosome\tposition\treference\talternate\tPAI3D_score\n" > ${TMP_DIR}/tmp_pai_header
  cat tmp_pai_header tmp_pai_scores | gzip -f > ${DATA_DIR}/primate_ai_3d_scores_tiered.tsv.gz

  rm tmp_pai_header tmp_pai_scores

elif [[ "$PAI3D_ACCESS" == "NO_PAI3D_ACCESS" ]]
then
  echo "Not gathering PrimateAI-3D scores for analysis!"
  exit 1
else
  echo "Not gathering PrimateAI-3D scores for analysis!"
  exit 1
fi

####################################################################
###### ---------------------------------------------------- ########
####################################################################

###########################################################################
###########################################################################
###                                                                     ###
###                         END OF SCRIPT                               ###
###                                                                     ###
###########################################################################
###########################################################################
