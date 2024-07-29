#!/usr/bin/env Rscript

##-------------------------------------------------------------------------
## Configs + packages

library(tidyverse)
library(data.table)
library(magrittr)

source("../config/config.R")

###########################################################################
###########################################################################
###                                                                     ###
###                       SECTION 1 FUNCTIONS                           ###
###                                                                     ###
###########################################################################
###########################################################################

add_varkey <- function(data, varid = TRUE) {
  data <- data %>%
    mutate(varkey = paste0("chr", chromosome, ":",
                           position, "_", reference, "_", alternate))
  if (varid) {
    data <- data %>%
      mutate(varid = paste0(participant_id, "_", varkey))
  }
  return(data)
}

acmg_tiers_adder <- function(exit, tiers, dd = FALSE,
                             exclude_variants = NULL,
                             diagnosis_definition = c("yes", "partially")) {

  tiers %<>% dplyr::select(-c(participant_id, varkey))

  ## Diagnostic discovery?
  if(dd) {
    exit %<>%
      filter(genome_build == "GRCh38") %>%
      mutate(acmg_classification = "likely_pathogenic_variant",
             case_solved_family = "diagnostic_discovery")
  }

  exit %<>%
    left_join(tiers, by = "varid") %>%
    mutate(acmg_broad = case_when(acmg_classification %in% c("likely_pathogenic_variant",
                                                             "pathogenic_variant") ~ "Pathogenic",
                                  acmg_classification == c("variant_of_unknown_clinical_significance") ~ "VUS",
                                  acmg_classification %in% c("likely_benign_variant",
                                                             "benign_variant") ~ "Benign"),
           tier = ifelse(is.na(tier), "UNTIERED", tier),
           tier = paste(tier, "ACMG", acmg_broad, sep = "_")) %>%
    dplyr::select(varkey, participant_id, varid, tier,
                  acmg_broad, case_solved_family) %>%
    distinct()

  if (!is.null(exclude_variants)) {
    exit %<>%
      filter(!varid %in% exclude_variants)
  }
  return(exit)
}

tiering_frequency_filter <- function(tiers, freq_data, rarity = "all") {
  ### Filtering thresholds are defined by the Rare Disease Guide
  variants_ultra_rare <- freq_data %>%
    filter_at(vars(c("AF_nfe_gnomad_v2",
                     "AF_afr_gnomad_v2",
                     "AF_sas_gnomad_v2",
                     "AF_amr_gnomad_v2",
                     "AF_eas_gnomad_v2",
                     "AF_asj_gnomad_v2",
                     "AF_oth_gnomad_v2",
                     "AF_fin_gnomad_v2",
                     "AF_afr_gnomad_v3",
                     "AF_amr_gnomad_v3",
                     "AF_eas_gnomad_v3",
                     "AF_sas_gnomad_v3",
                     "AF_fin_gnomad_v3",
                     "AF_nfe_gnomad_v3",
                     "uk10k_alspac",
                     "uk10k_twinsuk",
                     "gel_gl_6628_af")), all_vars(. < 0.001)) %>%
    filter_at(vars(c("AF_oth_gnomad_v3",
                     "AF_asj_gnomad_v3")), all_vars(. < 0.002)) %>%
    pull(varkey)

  variants_missing <- freq_data %>%
    filter_at(vars(contains("_")), all_vars(. == 0)) %>%
    pull(varkey)

  variants_missing_gnomad <- freq_data %>%
    filter_at(vars(contains("_gnomad")), all_vars(. == 0)) %>%
    pull(varkey)

  if (rarity == "all") {
    variants <- freq_data$varkey

  } else if (rarity == "rare") {
    variants <- freq_data %>%
      filter(!varkey %in% variants_ultra_rare) %>%
      pull(varkey)

  } else if (rarity == "ultra-rare") {
    variants <- variants_ultra_rare

  } else if (rarity == "missing") {
    variants <- variants_missing
  } else if (rarity == "missing-gnomad") {
    variants <- variants_missing_gnomad
  }
  ## Output
  tiers %>%
    filter(varkey %in% variants)
}

# forward method by James Ware:
# given a desired AF to filter on for Mendelian analyses (e.g. AF = 1%), what is the max AC that would be expected in a dataset of a given AN?
find_max_ac = function(af, an, ci = .95) {
  if (af == 0) {
    return (0)
  } else {
    quantile_limit = ci # ci for one-sided, 1-(1-ci)/2 for two-sided
    max_ac = qpois(quantile_limit, an*af)
    return (max_ac)
  }
}

find_af_filter = function(ac, an, ci=.95, lower=(.1/(2*50000)), upper=2, tol=1e-7, precision=1e-6) {
  ### Avoid singletons
  if (is.na(ac) | is.na(an) | ac == 0 | an == 0 | ac == 1) {
    return (0.0)
  } else {
    quantile_limit = ci # ci for one-sided, 1-(1-ci)/2 for two-sided
    attempt_uniroot = tryCatch({
      uniroot_result = uniroot(f = function(af,ac,an) { return (ac - 1 - qpois(p=quantile_limit,lambda=an*af)) },lower=lower,upper=upper,ac=ac,an=an,tol=tol)
    }, warning = function(w) {
      print(paste("ac= ",as.character(ac),", an= ",as.character(an)," warning = ",as.character(w),sep=''))
      return (0.0)
    }, error = function(e) {
      print(paste("ac= ",as.character(ac),", an= ",as.character(an)," error = ",as.character(e),sep=''))
      return (0.0)
    }, finally = {
    })
    max_af = round(uniroot_result$root,-log10(precision)) # round to nearest millionth
    while(find_max_ac(af = max_af, an = an) < ac) {
      max_af = max_af + precision # increase by millionths until you find the upper bound - the highest AF for which 95%CI AC is still less than observed AC
    }
    max_af = max_af - precision # back off one unit from the AF that violated the 95%CI AC < obs AC condition
    return (max_af)
  }
}

# This finds the af_filter for, say, a vector of AC at 10,000 differnet alleles and a corresponding vector of 10,000 AN values
# it's basically just a wrapper for mapply, which is trivial but makes more readable the code below in find_highest_af_filter.
# For applying to pop-specific COVID-19 AN/ACs...
find_af_filter_vectorized = function(ac_vector, an_vector, ci=.95, lower=(.1/(2*60706)), upper=2, tol=1e-7) { 
  return (mapply(find_af_filter, ac_vector, an_vector, ci=ci, lower=lower, upper=upper, tol=tol))
}

count_common_pops <- function(data, AF = 0.01, pop_id = "_") {
  data_afs <- data %>%
    dplyr::select(contains(pop_id))
  af_filter_pops <- lapply(apply(data_afs, 1,
                         function(x) which(x >= AF)),
                         names)
  pops_unique <- colnames(data_afs)
  pop_common_counts <- data.frame(ref_pop = NA_character_, num_instances = NA_integer_, num_instances_unique = NA_integer_)
  for (index in 1:length(pops_unique)) {
    pop_common_counts[index, 1] <- pops_unique[index]
    pop_common_counts[index, 2] <- length(grep(pops_unique[index], af_filter_pops))
    pop_common_counts[index, 3] <- length(grep(pops_unique[index], af_filter_pops[sapply(af_filter_pops, function(i) length(i) == 1)]))
  }
  return(pop_common_counts)
}

###########################################################################
###########################################################################
###                                                                     ###
###                       SECTION 2 DATA INGEST                         ###
###                      AF / FAF95 CALCULATIONS                        ###
###                                                                     ###
###########################################################################
###########################################################################

##### Meta data combined from various sources...
meta_combined_data <- read_csv(paste0(out_dir, "combined_data/combined_meta_data.csv.gz"))
nrow(meta_combined_data)
### 61,512 individuals

probands_meta_combined_data <- meta_combined_data %>%
  filter(participant_type == "Proband")
nrow(probands_meta_combined_data)
### 29,425 families / probands

#### Tiering frequency data
tiering_frequency_data <- read_csv(paste0(out_dir, "combined_data/combined_tiering_frequency_data.csv.gz"))

#### Tiering data for probands (filtered)
tiering_probands_filtered_data <- read_csv(paste0(out_dir, "combined_data/filtered_retiering_probands_data.csv.gz")) %>%
  filter(participant_id %in% probands_meta_combined_data$participant_id)

#### Only ultra-rare variants (exclude those >0.01% frequency in gnomAD)
#### With some sample size adjustments according to the Rare Disease Analysis Guide
tiering_probands_ur <- tiering_probands_filtered_data %>%
  tiering_frequency_filter(tiering_frequency_data, "ultra-rare")

#### Only ultra-rare variants (exclude those >0.01% frequency in gnomAD)
#### With some sample size adjustments according to the Rare Disease Analysis Guide
tiering_probands_miss <- tiering_probands_filtered_data %>%
  tiering_frequency_filter(tiering_frequency_data, "missing")

### COVID AFs for Tiered variants
covid_afs <- read_tsv(paste0(out_dir, "covid_afs_data/aggCOVID_tiering.acount.gz")) %>%
  rename(varkey = ID) %>%
  mutate(total_OBS_CTS = dplyr::select(., ends_with("OBS_CTS")) %>% rowSums(na.rm = TRUE),
         total_ALT_CTS = dplyr::select(., ends_with("ALT_CTS")) %>% rowSums(na.rm = TRUE),
         total_AF = total_ALT_CTS / total_OBS_CTS,
         missingness = 1 - (total_OBS_CTS / max(total_OBS_CTS))) %>%
  filter(missingness < 0.1 & total_AF < 0.01) %>%
  dplyr::select(-contains(c("Japan", "Africa_North", "South_America", "Finland"))) %>%
  filter(total_ALT_CTS > 0) %>%
  filter(varkey %in% tiering_probands_filtered_data$varkey)

##### Annotate exit data with tiering info...
exit_data <- read_csv(paste0(out_dir, "labkey_data/exit_data.csv.gz")) %>%
  add_varkey()

### Ultra-rare
exit_tiers_ur_data <- exit_data %>%
  acmg_tiers_adder(tiers = tiering_probands_ur)

###########################################################################
###########################################################################
###                                                                     ###
###                       SECTION 3 ANALYSIS                            ###
###                  FILTERING cPAVs + P/LP variants                    ###
###                                                                     ###
###########################################################################
###########################################################################

#### FAF95 calculation
cat(paste0("Performing FAF95 calculation for ", ncol(covid_afs %>% dplyr::select(contains("ALT_CTS"))) - 1, " groups in the cohort.\n"))
cat(paste0("Running for ", covid_afs %>% pull(varkey) %>% unique() %>% length()), "variants in the cohort triaged as cPAVs in the 100kGP.\n")
ac_list <- as.list(covid_afs %>% dplyr::select(contains("ALT_CTS")))
an_list <- as.list(covid_afs %>% dplyr::select(contains("OBS_CTS")))
af_filter_list = mapply(find_af_filter_vectorized, ac_list, an_list, SIMPLIFY=FALSE)
af_filter_data <- data.frame(af_filter_list) %>%
  rename_with(~gsub("ALT_CTS", "FAF95", .x, fixed = TRUE))

#### Write FAF95 CSV file
af_filter_fafs <- cbind(covid_afs, af_filter_data) %>%
  mutate(Africa_East_AF = Africa_East_ALT_CTS / Africa_East_OBS_CTS,
       Africa_South_AF = Africa_South_ALT_CTS / Africa_South_OBS_CTS,
       Africa_West_AF = Africa_West_ALT_CTS / Africa_West_OBS_CTS,
       Middle_East_AF = Middle_East_ALT_CTS / Middle_East_OBS_CTS,
       Mediterranean_AF = Mediterranean_ALT_CTS / Mediterranean_OBS_CTS,
       Ashkenazi_AF = Ashkenazi_ALT_CTS / Ashkenazi_OBS_CTS,
       Europe_North_West_AF = Europe_North_West_ALT_CTS / Europe_North_West_OBS_CTS,
       Europe_South_West_AF = Europe_South_West_ALT_CTS / Europe_South_West_OBS_CTS,
       Europe_East_AF = Europe_East_ALT_CTS / Europe_East_OBS_CTS,
       Asia_East_AF = Asia_East_ALT_CTS / Asia_East_OBS_CTS,
       Philippines_AF = Philippines_ALT_CTS / Philippines_OBS_CTS,
       Pakistan_AF = Pakistan_ALT_CTS / Pakistan_OBS_CTS,
       Bangladesh_AF = Bangladesh_ALT_CTS / Bangladesh_OBS_CTS,
       Sri_Lanka_AF = Sri_Lanka_ALT_CTS / Sri_Lanka_OBS_CTS
)

## Variants >1% in at least one GIA group in COVID-19 cohort
common_varkeys <- af_filter_fafs %>%
  dplyr::select(-total_FAF95) %>%
  filter_at(vars(ends_with("FAF95")), any_vars(. >= 0.01)) %>%
  pull(varkey)

## Common cPAVs IDd as ultra-rare
tiering_ur_common <- tiering_probands_ur %>%
  filter(varkey %in% common_varkeys)
nrow(tiering_ur_common)
# 13,208 ultra-rare cPAVs common in COVID-19 cohort

length(unique(tiering_ur_common$varkey))
# 2,046 distinct ultra-rare cPAVs common in COVID-19 cohort

#### Common GPcPAVs IDd as ultra-rare
gpcpav_ur_common <- tiering_probands_ur %>%
  filter(tier %in% c("TIER2", "TIER1")) %>%
  filter(varkey %in% common_varkeys)

nrow(gpcpav_ur_common)
# 173
length(unique(gpcpav_ur_common$varkey))
# 111

exit_tiers_ur_data %>%
  filter(varkey %in% gpcpav_ur_common$varkey) %>%
  drop_na(acmg_broad)
# No P/LP or B/LB variants or submitted to diagnostic discovery

## In which GIA groups in GenOMICC are these occurring?
####### Pop-specific variation
pop_common_tiering <- af_filter_fafs %>%
  filter(varkey %in% tiering_ur_common$varkey) %>%
  dplyr::select(-total_FAF95) %>%
  count_common_pops(AF = 0.01, pop_id = "_FAF95") %>%
  mutate(participant_genetic_category = gsub("_", " ", gsub("_FAF95", "", ref_pop)),
         test = "rare_cPAV_common_FAF95")

pop_unique_counts <- af_filter_fafs %>%
  filter(varkey %in% tiering_ur_common$varkey) %>%
  dplyr::select(-total_ALT_CTS) %>%
  count_common_pops(AF = 1, pop_id = "ALT_CTS")

## Variants >0.1% in at least one GIA group in COVID-19 cohort
rare_varkeys <- af_filter_fafs %>%
  dplyr::select(-total_FAF95) %>%
  filter_at(vars(ends_with("FAF95")), any_vars(. >= 0.001)) %>%
  pull(varkey)

## Common cPAVs IDd as ultra-rare
tiering_miss_rare <- tiering_probands_miss %>%
  filter(varkey %in% rare_varkeys)
nrow(tiering_miss_rare)
# 4,961 ultra-rare common cPAVs

length(unique(tiering_miss_rare$varkey))
# 1,941 distinct ultra-rare common cPAVs

## Common cPAVs IDd as ultra-rare
tiering_ur_rare <- tiering_probands_ur %>%
  filter(varkey %in% rare_varkeys)
nrow(tiering_ur_rare)
# 4,961 ultra-rare common cPAVs

length(unique(tiering_ur_rare$varkey))
# 1,941 distinct ultra-rare common cPAVs

#### Common GPcPAVs IDd as ultra-rare
gpcpav_miss_rare <- tiering_probands_miss %>%
  filter(tier %in% c("TIER2", "TIER1")) %>%
  filter(varkey %in% rare_varkeys)

pop_miss_rare_tiering <- af_filter_fafs %>%
  filter(varkey %in% tiering_miss_rare$varkey) %>%
  dplyr::select(-total_FAF95) %>%
  count_common_pops(AF = 0.001, pop_id = "_FAF95") %>%
  mutate(participant_genetic_category = gsub("_", " ", gsub("_FAF95", "", ref_pop)),
         test = "missing_cPAV_uncommon_FAF95")

#########################################################################
#########################################################################
###                                                                   ###
###                          WRITE DATA                               ###
###                                                                   ###
#########################################################################
#########################################################################

## Analysis data
data_dir <- paste0(analysis_dir, "covid_filter_analysis/")

## Create directory.
dir.create(file.path(data_dir))
cat(paste0("Writing analyses to: ", data_dir, "\n"))

## Number of filtered cPAVs / GPcPAVs data table
rbind(pop_common_tiering, pop_miss_rare_tiering) %>%
  write.csv(paste0(data_dir, "aggCOVID_faf95_filter_tiering.csv"),
            row.names = FALSE)

###########################################################################
###########################################################################
###                                                                     ###
###                         END OF SCRIPT                               ###
###                                                                     ###
###########################################################################
###########################################################################
