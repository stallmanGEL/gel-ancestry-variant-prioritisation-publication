##-------------------------------------------------------------------------
## Configs + packages

source ../config/config.cfg

module purge
module load bcftools/1.16
module load tabix/1.18

### Variables
DATA_DIR=${OUT_DIR}/relate_data

## Some library path stuff to get Relate running
LD_LIBRARY_PATH="/resources/tools/apps/software/lang/Python/3.7.4-GCCcore-8.3.0/lib:/resources/tools/apps/software/R/4.0.2-foss-2019b/lib64/R/lib:/resources/tools/apps/software/R/4.0.2-foss-2019b/lib64:/resources/tools/apps/software/vis/ImageMagick/7.0.9-5-GCCcore-8.3.0/lib:/resources/tools/apps/software/vis/LittleCMS/2.9-GCCcore-8.3.0/lib:/resources/tools/apps/software/vis/JasPer/2.0.14-GCCcore-8.3.0/lib64:/resources/tools/apps/software/tools/Ghostscript/9.50-GCCcore-8.3.0/lib:/resources/tools/apps/software/numlib/GSL/2.6-GCC-8.3.0/lib:/resources/tools/apps/software/phys/UDUNITS/2.2.26-GCCcore-8.3.0/lib:/resources/tools/apps/software/data/HDF5/1.10.5-gompi-2019b/lib:/resources/tools/apps/software/tools/Szip/2.1.1-GCCcore-8.3.0/lib:/resources/tools/apps/software/ICU/65.1-GCCcore-8.3.0/lib:/resources/tools/apps/software/lib/libsndfile/1.0.28-GCCcore-8.3.0/lib:/resources/tools/apps/software/numlib/NLopt/2.6.1-GCCcore-8.3.0/lib64:/resources/tools/apps/software/tools/cURL/7.66.0-GCCcore-8.3.0/lib:/resources/tools/apps/software/vis/Tk/8.6.9-GCCcore-8.3.0/lib:/resources/tools/apps/software/lang/Java/11.0.2/lib:/resources/tools/apps/software/lib/LibTIFF/4.0.10-GCCcore-8.3.0/lib:/resources/tools/apps/software/lib/libjpeg-turbo/2.0.3-GCCcore-8.3.0/lib64:/resources/tools/apps/software/PCRE2/10.33-GCCcore-8.3.0/lib:/resources/tools/apps/software/devel/SQLite/3.29.0-GCCcore-8.3.0/lib:/resources/tools/apps/software/lang/Tcl/8.6.9-GCCcore-8.3.0/lib:/resources/tools/apps/software/lib/libreadline/8.0-GCCcore-8.3.0/lib:/resources/tools/apps/software/vis/cairo/1.16.0-GCCcore-8.3.0/lib:/resources/tools/apps/software/vis/GLib/2.62.0-GCCcore-8.3.0/lib:/resources/tools/apps/software/devel/PCRE/8.43-GCCcore-8.3.0/lib:/resources/tools/apps/software/tools/gettext/0.20.1-GCCcore-8.3.0/lib:/resources/tools/apps/software/lib/libffi/3.2.1-GCCcore-8.3.0/lib64:/resources/tools/apps/software/lib/libffi/3.2.1-GCCcore-8.3.0/lib:/resources/tools/apps/software/vis/pixman/0.38.4-GCCcore-8.3.0/lib:/resources/tools/apps/software/vis/libGLU/9.0.1-GCCcore-8.3.0/lib:/resources/tools/apps/software/vis/Mesa/19.1.7-GCCcore-8.3.0/lib:/resources/tools/apps/software/lib/libunwind/1.3.1-GCCcore-8.3.0/lib:/resources/tools/apps/software/compiler/LLVM/9.0.0-GCCcore-8.3.0/lib:/resources/tools/apps/software/lib/libdrm/2.4.99-GCCcore-8.3.0/lib:/resources/tools/apps/software/lib/nettle/3.5.1-GCCcore-8.3.0/lib64:/resources/tools/apps/software/math/GMP/6.1.2-GCCcore-8.3.0/lib:/resources/tools/apps/software/vis/X11/20190717-GCCcore-8.3.0/lib:/resources/tools/apps/software/vis/fontconfig/2.13.1-GCCcore-8.3.0/lib:/resources/tools/apps/software/tools/util-linux/2.34-GCCcore-8.3.0/lib:/resources/tools/apps/software/devel/ncurses/6.1-GCCcore-8.3.0/lib:/resources/tools/apps/software/vis/freetype/2.10.1-GCCcore-8.3.0/lib:/resources/tools/apps/software/lib/libpng/1.6.37-GCCcore-8.3.0/lib:/resources/tools/apps/software/tools/expat/2.2.7-GCCcore-8.3.0/lib:/resources/tools/apps/software/tools/bzip2/1.0.8-GCCcore-8.3.0/lib:/resources/tools/apps/software/numlib/ScaLAPACK/2.0.2-gompi-2019b/lib:/resources/tools/apps/software/numlib/FFTW/3.3.8-gompi-2019b/lib:/resources/tools/apps/software/numlib/OpenBLAS/0.3.7-GCC-8.3.0/lib:/resources/tools/apps/software/mpi/OpenMPI/3.1.4-GCC-8.3.0/lib:/resources/tools/apps/software/system/hwloc/1.11.12-GCCcore-8.3.0/lib:/resources/tools/apps/software/system/libpciaccess/0.14-GCCcore-8.3.0/lib:/resources/tools/apps/software/lib/libxml2/2.9.9-GCCcore-8.3.0/lib:/resources/tools/apps/software/tools/XZ/5.2.4-GCCcore-8.3.0/lib:/resources/tools/apps/software/tools/numactl/2.0.12-GCCcore-8.3.0/lib:/resources/tools/apps/software/tools/binutils/2.32-GCCcore-8.3.0/lib:/resources/tools/apps/software/lib/zlib/1.2.11-GCCcore-8.3.0/lib:/resources/tools/apps/software/compiler/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0:/resources/tools/apps/software/compiler/GCCcore/8.3.0/lib64:/resources/tools/apps/software/compiler/GCCcore/8.3.0/lib:/usr/share/lsf/10.1/linux3.10-glibc2.17-x86_64/lib"

## Go to the relate analysis directory
cd $TMP_DIR
if [ ! -d ./phased_vcfs ]
then
  mkdir -p ./phased_vcfs
fi

###########################################################################
###########################################################################
###                                                                     ###
###                          ANALYSIS RUN                               ###
###                                                                     ###
###########################################################################
###########################################################################

### Population labels for Relate (10 perfectly phased probands per pop)
tail -n+2 ${OUT_DIR}/combined_data/relate_ss_poplabels.tsv | awk '{print $1}' > relate_ss.inds

### Run script per CHROM
for CHROM in `seq 1 22`
do

  ## Get SNPs from phased aggV2 VCF file...
  bcftools view --min-ac 1 \
--threads 8 \
--phased \
--types snps \
--samples-file relate_ss.inds \
--output-type z \
--output ./phased_vcfs/phased_panel_relate_subset_chr${CHROM}.vcf.gz \
${PHASED_AGGV2_VCF}

  ## Index
  tabix -f ./phased_vcfs/phased_panel_relate_subset_chr${CHROM}.vcf.gz

  ## Convert to HAPS/SAMPLE
  $PATH_TO_RELATE/bin/RelateFileFormats \
--mode ConvertFromVcf \
--haps ./phased_vcfs/phased_panel_relate_subset_chr${CHROM}.haps \
--sample ./phased_vcfs/phased_panel_relate_subset_chr${CHROM}.sample \
-i ./phased_vcfs/phased_panel_relate_subset_chr${CHROM} \
--chr ${CHROM}

  ## Remove the chr prefix to align with other files...
  sed -i 's/chr//g' ./phased_vcfs/phased_panel_relate_subset_chr${CHROM}.haps

  $PATH_TO_RELATE/scripts/PrepareInputFiles/PrepareInputFiles.sh \
--haps ./phased_vcfs/phased_panel_relate_subset_chr${CHROM}.haps \
--sample ./phased_vcfs/phased_panel_relate_subset_chr${CHROM}.sample \
--ancestor $ANCESTRAL_GENOME \
--poplabels ${OUT_DIR}/combined_data/relate_ss_poplabels.tsv \
--mask $TGP_MASK \
-o ./phased_vcfs/phased_panel_relate_subset_anc_chr${CHROM}

  ## Make a genetics map directory
  if [ ! -d ./genetic_maps ]
  then
    mkdir -p ./genetic_maps
  fi

  ## Change genetic map accordingly
  awk '{print $2,$3,$4}' $GENETIC_MAP > ./genetic_maps/chr${CHROM}.b38.map
  cd phased_vcfs/

  ### Create trees in parallel
  $PATH_TO_RELATE/scripts/RelateParallel/RelateParallel.sh \
-m 1.25e-8 \
-N 30000 \
--haps phased_panel_relate_subset_anc_chr${CHROM}.haps.gz \
--sample phased_panel_relate_subset_anc_chr${CHROM}.sample.gz \
--map ../genetic_maps/chr${CHROM}.b38.map \
--dist phased_panel_relate_subset_anc_chr${CHROM}.dist.gz \
--seed 1 \
-o phased_panel_relate_subset_anc_chr${CHROM} \
--threads 8

done

###### Run popsizes script
## Set seed (1)
$PATH_TO_RELATE/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
-i phased_vcfs/phased_panel_relate_subset_anc \
-m 1.25e-8 \
--poplabels ${OUT_DIR}/combined_data/relate_ss_poplabels.tsv  \
--seed 1 \
--first_chr 1 \
--last_chr 22 \
-o ${DATA_DIR}/phased_panel_relate_subset_anc \
--threads 8

###########################################################################
###########################################################################
###                                                                     ###
###                         END OF SCRIPT                               ###
###                                                                     ###
###########################################################################
###########################################################################