#!/bin/bash

##-------------------------------------------------------------------------
## Configs + packages

source ../config/config.sh

module purge
module load bcftools/1.16
module load tabix/1.18
module load plink/2.0

##-------------------------------------------------------------------------
## Which file set to run for (aggV2 or aggCOVIDv5)

### COVID or AGGV2
FILE_SET=$1

if [[ "$FILE_SET" == "aggCOVIDv5" ]]
then
  CHUNK_LIST=../config/aggCOVIDv5.file_list
  OUT_FILE_PREFIX=aggCOVID_all
elif [[ "$FILE_SET" == "aggV2" ]]
then
  CHUNK_LIST=../config/aggV2.file_list
  OUT_FILE_PREFIX=aggV2_all
else
  echo "Need to specify which dataset to run (aggV2 or aggCOVIDv5) as string alongside script."
  exit 1
fi

###########################################################################
###########################################################################
###                                                                     ###
###                           DATA CREATE                               ###
###                                                                     ###
###########################################################################
###########################################################################

### Variables
DATA_DIR=${OUT_DIR}/pavs_gnomad_miss_data

## Run per file chunk
while read chunk
do
  CHROM=$(echo $chunk | awk '{print $1}')
  START=$(echo $chunk | awk '{print $2}')
  END=$(echo $chunk | awk '{print $3}')
  ANNOT_VCF=$(echo $chunk | awk '{print $4}')
  GT_PGEN=$(echo $chunk | awk '{print $5}')

  ### Output prefix
  TMP_PREFIX=${TMP_DIR}/${CHROM}_${START}_${END}

  ### Get PAVs targets file
  ### Use the annotations that define "Protein-altering" in the RD pipeline
  bcftools +split-vep \
-f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%CSQ\n' \
-i "%FILTER='PASS' | %FILTER='.'" \
-s worst \
-d \
-A tab ${ANNOT_VCF} | \
grep -P "splice_region_variant|MODERATE|HIGH|incomplete_terminal_codon_variant|initiator_codon_variant" | \
awk '{print $1 ":" $2 "_" $3 "_" $4}' \
> ${TMP_PREFIX}.pav.ids

  ### Filter out variants in the chunk that are present in gnomADv2 exomes or gnomADv3 genomes
  ### Evaluate only PASS variants in GEL aggregates, but include all variants in gnomAD aggregates.
  bcftools isec \
-i "%FILTER='PASS' | %FILTER='.'" \
-i- \
-i- \
-C \
-t ${CHROM}:${START}-${END} \
${ANNOT_VCF} \
${GNOMAD_V2_VCF} \
${GNOMAD_V3_VCF} | \
awk '{print $1 ":" $2 "_" $3 "_" $4}' \
> ${TMP_PREFIX}.gnomAD_missing.ids

  ### Get only shared files based on variant ID (modelled on those in the PLINK PGEN files)
  grep -Fxf ${TMP_PREFIX}.pav.ids ${TMP_PREFIX}.gnomAD_missing.ids > ${TMP_PREFIX}.gnomAD_missing_pav.ids

  ### Number of missing variants per individual
  plink2 --pfile ${GT_PGEN} \
--sample-counts \
'cols=homref,homalt,het,nonsnp,single,missing' \
--extract ${TMP_PREFIX}.gnomAD_missing_pav.ids \
--out ${TMP_PREFIX}_gnomAD_missing_pavs

  ## Get file header (replace each time)
  head -1 ${TMP_PREFIX}_gnomAD_missing_pavs.scount > ${DATA_DIR}/${OUT_FILE_PREFIX}_gnomAD_missing_pavs.scount

done < $CHUNK_LIST

### Combine chunks
for scount_chunk in ${TMP_PREFIX}*_gnomAD_missing_pavs.scount
do
  tail -n+2 $scount_chunk >> ${DATA_DIR}/${OUT_FILE_PREFIX}_gnomAD_missing_pavs.scount
done

## Compress with gzip
gzip -f ${DATA_DIR}/${OUT_FILE_PREFIX}_gnomAD_missing_pavs.scount

###########################################################################
###########################################################################
###                                                                     ###
###                         END OF SCRIPT                               ###
###                                                                     ###
###########################################################################
###########################################################################
