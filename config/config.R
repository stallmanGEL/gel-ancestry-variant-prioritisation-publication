##########################
## WILDCARD WEKA CONFIG ##
##########################

analysis_dir <- "/pgen_int_work/BRS/dd/sam/data/gel-ancestry-variant-prioritisation-publication/analysis/"
plots_dir <- "/pgen_int_work/BRS/dd/sam/data/gel-ancestry-variant-prioritisation-publication/plots/"
tmpDir <- "/re_scratch/dd/sam/temp/"

#### ONLY CHANGE THIS TO YOUR PERSONAL DIRECTORY IF YOU WISH TO REPLICATE ALL DATA GENERATION (WILL TAKE SOME TIME!!)
out_dir <- "/pgen_int_work/BRS/dd/sam/data/gel-ancestry-variant-prioritisation-publication/data/"
####################################################################################################################

#####################################
##### FILES FOR DATA GENERATION #####
#####################################

hq_snps_plink <- "/pgen_int_work/BRS/covid-19/aggCOVID_v5_aggV2/HQ_SNPs/aggCOVID_v5_aggV2_HQSNPs"
kinship_matrix <- "/pgen_int_work/BRS/covid-19/aggCOVID_v5_aggV2/relatedness/aggCOVID_v5_aggV2_HQSNPs.kin0"

#### Uncomment the below if you are unable to access aggCOVIDv5 in the RE BUT WISH TO RUN DATA GENERATION.
#### NOTE: OUTPUTS MAY NOT BE IDENTICAL WITH THE PUBLISHED ANALYSIS.

# hq_snps_plink <- "/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/additional_data/HQ_SNPs/GELautosomes_LD_pruned_1kgp3Intersect_maf0.05_mpv10"
# kinship_matrix <- "/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/additional_data/PCs_relatedness/relatedness/GEL_aggV2_MAF5_mp10_0.0442.kin0"

prive_ukbb_build38_afs <- "/public_data_resources/ukbb_prive/uk_biobank_1kg_grch38_freqs.csv.gz"
prive_ukbb_build38_loadings <- "/public_data_resources/ukbb_prive/uk_biobank_1kg_grch38_projection.csv.gz"

###################
## LABKEY CONFIG ##
###################

labkey_v17 <- "/main_programme/main-programme_v17_2023-03-30"
labkey_v18 <- "/main_programme/main-programme_v18_2023-12-21"
baseUrl <- "https://labkey-embassy.gel.zone/labkey/"

### Set seed for replication of randomness.
set.seed(1)

###############
##### END #####
###############
