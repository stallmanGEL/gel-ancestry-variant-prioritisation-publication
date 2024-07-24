#!/usr/bin/env bash

##-------------------------------------------------------------------------
## Configs + packages

source ../config/config.sh

###########################################################################
###########################################################################
###                                                                     ###
###                        DATA CREATE gnomADv3                         ###
###                                                                     ###
###########################################################################
###########################################################################

### Variables
DATA_DIR=${OUT_DIR}/gnomad_afs_data
TIERED_VARIANTS_TAB=${OUT_DIR}/labkey_data/tiering_sorted_variants.tsv

### Make key of tiered variants based on PLINK encoded variant IDs
sed 's/,/_/g' $TIERED_VARIANTS_TAB | \
awk '{print $1,":",$2,"_",$3}' | \
sed 's/ //g' > ${TMP_DIR}/tiering_sorted_variants_key.txt

### For gnomADv3, we need to loop through each chromosome
for CHROM in `seq 1 22`
do

  ### Then add a header (with a specific title, v3)
  ### Replace each time
  title=v3
  tmp_header=${TMP_DIR}/tmp_header
  printf "chromosome\tposition\trsid\treference\talternate\tfilter\tAN_gnomad_${title}\tAF_gnomad_${title}\tAN_afr_gnomad_${title}\tAC_afr_gnomad_${title}\tAN_asj_gnomad_${title}\tAC_asj_gnomad_${title}\tAN_amr_gnomad_${title}\tAC_amr_gnomad_${title}\tAN_fin_gnomad_${title}\tAC_fin_gnomad_${title}\tAN_nfe_gnomad_${title}\tAC_nfe_gnomad_${title}\tAN_sas_gnomad_${title}\tAC_sas_gnomad_${title}\tAN_eas_gnomad_${title}\tAC_eas_gnomad_${title}\tAN_oth_gnomad_${title}\tAC_oth_gnomad_${title}\n" > $tmp_header

  ### Search for the relevant AN/AC terms 
  ### Only take the PASS variants...
  tmp_afs=${TMP_DIR}/tmp_${CHROM}_afs
  zcat < $GNOMAD_V3_VCF | sed '/^#/d' | \
sed 's/;/ /g' | \
grep "PASS" | \
awk -v search="AN,AF,AN_afr,AC_afr,AN_asj,AC_asj,AN_amr,AC_amr,AN_fin,AC_fin,AN_nfe,AC_nfe,AN_sas,AC_sas,AN_eas,AC_eas,AN_oth,AC_oth" \
'BEGIN {OFS="\t"} { n=split(search, arr, /,/) } { for(i=1; i in arr; i++) printf("%s%s", (match($0,"(^| )"arr[i]"=[^ ]*") ? substr($0,(RSTART>1?RSTART+1:RSTART),(RSTART>1?RLENGTH-1:RLENGTH)) : ""), i==n ? ORS : OFS) }' > $tmp_afs

  ### Then add the relevant bits of the VCF file
  tmp_sites=${TMP_DIR}/tmp_${CHROM}_sites
zcat < $GNOMAD_V3_VCF | sed '/^#/d' | \
grep "PASS" | \
awk 'OFS="\t" {print $1,$2,$3,$4,$5,$7}' > $tmp_sites

  ### Paste them together
  paste $tmp_sites $tmp_afs | \
sed 's/AN=//g' | \
sed 's/AF=//g' | \
sed 's/AC_afr=//g' | \
sed 's/AC_fin=//g' | \
sed 's/AC_amr=//g' | \
sed 's/AC_asj=//g' | \
sed 's/AC_sas=//g' | \
sed 's/AC_eas=//g' | \
sed 's/AC_nfe=//g' | \
sed 's/AC_oth=//g' | \
sed 's/AN_afr=//g' | \
sed 's/AN_fin=//g' | \
sed 's/AN_amr=//g' | \
sed 's/AN_asj=//g' | \
sed 's/AN_sas=//g' | \
sed 's/AN_eas=//g' | \
sed 's/AN_nfe=//g' | \
sed 's/AN_oth=//g' | \
sed 's/ /\t/g' > ${TMP_DIR}/gnomad_v3_genomes_chr${CHROM}.afs

  ### Remove the temporary files.
  rm $tmp_afs $tmp_sites
done

### Combine chromosomes
for gnomad_chunk in ${TMP_DIR}/gnomad_v3_genomes_chr*.afs
do
  tail -n+2 $gnomad_chunk >> ${DATA_DIR}/gnomad_v3_genomes_all.afs
done

## Compress with gzip
gzip -f ${DATA_DIR}/gnomad_v3_genomes_all.afs

### Extract only tiered variants...
zcat ${DATA_DIR}/gnomad_v3_genomes_all.afs.gz | \
awk 'OFS="\t" { $(NF+1) = $1":"$2"_"$4"_"$5 }1' | \
awk -F"\t" 'FNR==NR {hash[$1]; next} $(NF) in hash' ${TMP_DIR}/tiering_sorted_variants_key.txt - | \
gzip -f > ${DATA_DIR}/gnomad_v3_genomes_tiering.afs.gz 

###########################################################################
###########################################################################
###                                                                     ###
###                        DATA CREATE gnomADv2                         ###
###                                                                     ###
###########################################################################
###########################################################################

random_string=$(echo $RANDOM | md5sum | head -c 20; echo)

tmp_afs=${TMP_DIR}/tmp_${random_string}_afs
zcat < /public_data_resources/gnomad/2.1.1/grch38_lift_over/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz | sed '/^#/d' | \
sed 's/;/ /g' | \
grep "PASS" | \
awk -v search="AN,AF,AN_afr,AC_afr,AN_asj,AC_asj,AN_amr,AC_amr,AN_fin,AC_fin,AN_nfe,AC_nfe,AN_sas,AC_sas,AN_eas,AC_eas,AN_oth,AC_oth" \
'BEGIN {OFS="\t"} { n=split(search, arr, /,/) } { for(i=1; i in arr; i++) printf("%s%s", (match($0,"(^| )"arr[i]"=[^ ]*") ? substr($0,(RSTART>1?RSTART+1:RSTART),(RSTART>1?RLENGTH-1:RLENGTH)) : ""), i==n ? ORS : OFS) }' > $tmp_afs

tmp_sites=${TMP_DIR}/tmp_${random_string}_sites
zcat < /public_data_resources/gnomad/2.1.1/grch38_lift_over/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz | sed '/^#/d' | \
grep "PASS" | \
awk 'OFS="\t" {print $1,$2,$3,$4,$5,$7}' > $tmp_sites

title=v2
tmp_header=${TMP_DIR}/tmp_${random_string}_header
printf "chromosome\tposition\trsid\treference\talternate\tfilter\tAN_gnomad_${title}\tAF_gnomad_${title}\tAN_afr_gnomad_${title}\tAC_afr_gnomad_${title}\tAN_asj_gnomad_${title}\tAC_asj_gnomad_${title}\tAN_amr_gnomad_${title}\tAC_amr_gnomad_${title}\tAN_fin_gnomad_${title}\tAC_fin_gnomad_${title}\tAN_nfe_gnomad_${title}\tAC_nfe_gnomad_${title}\tAN_sas_gnomad_${title}\tAC_sas_gnomad_${title}\tAN_eas_gnomad_${title}\tAC_eas_gnomad_${title}\tAN_oth_gnomad_${title}\tAC_oth_gnomad_${title}\n" > $tmp_header

paste $tmp_sites $tmp_afs | \
sed 's/AN=//g' | \
sed 's/AF=//g' | \
sed 's/AC_afr=//g' | \
sed 's/AC_fin=//g' | \
sed 's/AC_amr=//g' | \
sed 's/AC_asj=//g' | \
sed 's/AC_sas=//g' | \
sed 's/AC_eas=//g' | \
sed 's/AC_nfe=//g' | \
sed 's/AC_oth=//g' | \
sed 's/AN_afr=//g' | \
sed 's/AN_fin=//g' | \
sed 's/AN_amr=//g' | \
sed 's/AN_asj=//g' | \
sed 's/AN_sas=//g' | \
sed 's/AN_eas=//g' | \
sed 's/AN_nfe=//g' | \
sed 's/AN_oth=//g' | \
sed 's/ /\t/g' | \
cat $tmp_header - | \
gzip -f > ${DATA_DIR}/gnomad_v2_exomes_all.afs.gz

rm $tmp_afs $tmp_header $tmp_sites

### Extract only tiered variants...
zcat ${DATA_DIR}/gnomad_v2_exomes_all.afs.gz | \
awk 'OFS="\t" { $(NF+1) = $1":"$2"_"$4"_"$5 }1' | \
awk -F"\t" 'FNR==NR {hash[$1]; next} $(NF) in hash' ${TMP_DIR}/tiering_sorted_variants_key.txt - | \
gzip -f > ${DATA_DIR}/gnomad_v2_exomes_tiering.afs.gz

###########################################################################
###########################################################################
###                                                                     ###
###                         END OF SCRIPT                               ###
###                                                                     ###
###########################################################################
###########################################################################