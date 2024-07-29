#!/usr/bin/env Rscript

##-------------------------------------------------------------------------
## Configs + packages

source("../config/config.R")

library(tidyverse)
library(data.table)
library(magrittr)
library(zoo)
library(broom)
library(MASS)

###########################################################################
###########################################################################
###                                                                     ###
###                           FUNCTIONS                                 ###
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

  if (rarity == "all") {
    variants <- freq_data$varkey

  } else if (rarity == "rare") {
    variants <- freq_data %>%
      filter(!varkey %in% variants_ultra_rare) %>%
      pull(varkey)

  } else if (rarity == "ultra-rare-or-missing") {
    variants <- variants_ultra_rare
  
  } else if (rarity == "ultra-rare") {
    variants <- variants_ultra_rare[which(!variants_ultra_rare %in% variants_missing)]

  } else if (rarity == "missing") {
    variants <- variants_missing
  }
  ## Output
  tiers %>%
    filter(varkey %in% variants)
}

tiers_multi_filter <- function(tiers, zygosity = NULL,
                               varkeys = NULL,
                               exclude_varkeys = NULL) {

  if (!is.null(varkeys)) {
    tiers %<>% filter(varkey %in% varkeys)
  }

  if (!is.null(exclude_varkeys)) {
    tiers %<>% filter(!varkey %in% exclude_varkeys)
  }

  if (!is.null(zygosity)) {
    if (zygosity %in% tiers$genotype) {
      tiers %<>% filter(genotype %in% zygosity)
      cat(paste0(
        "Summarising tiered variants that passed the ",
        zygosity,
        " filter only!\n"))
    } else {
      stop(paste0(
        "No individuals with genotype: ",
        zygosity, " in dataset.\n"))
    }
  }
  return(tiers)
}

pred_annot_adder <- function(tiers = NULL, exit = NULL,
                             annots, algo) {
  annot_pred <- rlang::sym(paste0(algo, "_pred"))
  annots_comb <- annots %>%
    dplyr::select(varkey, {{annot_pred}}) %>%
    filter({{annot_pred}} != "")
  if (!is.null(tiers)) {
    out <- tiers %>%
      left_join(annots_comb,
                by = "varkey") %>%
      drop_na({{annot_pred}}) %>%
      mutate(tier = paste(tier, algo, {{annot_pred}}, sep = "_"))
  } else if (!is.null(exit)) {
    out <- exit %>%
      left_join(annots, by = "varkey")
  }
  return(out)
}

acmg_tiers_adder <- function(exit, tiers, dd = FALSE,
                             exclude_variants = NULL,
                             diagnosis_definition = c("yes", "partially", "diagnostic_discovery")) {

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

### Variable
configure_tiering_counts <- function(meta, tiers, merge = c("TIER1", "TIER2"),
                                     merge_to = "GPcPAV", source = NULL) {
  merge_str <- paste(merge, collapse = "|")

  tiers_counts <- tiers %>%
    dplyr::select(participant_id, varkey, tier) %>%
    distinct() %>%
    filter(grepl(merge_str, tier)) %>%
    mutate(tier = gsub(merge_str, merge_to, tier)) %>%
    group_by(participant_id, tier) %>%
    tally() %>%
    pivot_wider(names_from = tier,
      values_from = n,
      values_fill = list(n = 0)) %>%
    ungroup()

  if (!is.null(source)) {
    merge_source_sym <- rlang::sym(paste(merge_to, source, sep = "_"))
    tiers_counts %<>%
      mutate({{merge_source_sym}} := rowSums(.[,-1]))
  }

  meta %>%
    left_join(tiers_counts, by = "participant_id") %>%
    mutate_at(vars(contains(merge_to)), replace_na, 0)
}

configure_diagnosed <- function(meta, patho_variants,
                                diagnosis_definition = c("yes", "partially"),
                                with = NULL,
                                diagnosed_label = "DIAGNOSED") {

  diagnosed_label_sym <- rlang::sym(diagnosed_label)
  diagnosed_inds <- patho_variants %>%
    filter(case_solved_family %in% diagnosis_definition) %>%
    pull(participant_id)

  meta %<>%
    mutate({{diagnosed_label_sym}} := ifelse(participant_id %in% diagnosed_inds,
                                             "Yes",
                                             "No"))

  if (!is.null(with)) {
    with_sym <- rlang::sym(with)
    diagnosed_with_sym <- rlang::sym(paste(diagnosed_label, with, sep = "_"))
    meta %<>%
      mutate({{diagnosed_with_sym}} := ifelse({{diagnosed_label_sym}} == "Yes" & {{with_sym}} > 0,
                                              "Yes", "No"))
  }
  return(meta)
}

relevel_covariates <- function(meta) {
  meta %<>%
    mutate(family_group_type = relevel(factor(family_group_type),
                                                  ref = "Singleton"),
           penetrance = relevel(factor(penetrance),
                                ref = "incomplete"),
           only_proband_affected = relevel(factor(only_proband_affected),
                                             ref = "No"))
}

check_covariate_variation <- function(df, covariates) {
  ## Remove covariate from list if it is...
  ## Not in the dataframe or not variable
  covar_vec_mult <- str_trim(unlist(strsplit(covariates, split = "[+]")))
  covar_vec <- str_trim(unlist(strsplit(covar_vec_mult, split = "[:]")))

  for (covar in covar_vec) {
    if (!covar %in% colnames(df)) {
      cat(paste0("Covariate ", covar,
                " does not exist. Removing...\n"))
      covar_vec <- covar_vec[!covar_vec %in% covar]
      covar_vec_mult <- covar_vec_mult[!grepl(covar, covar_vec_mult)]
    } else {
      covar_sym <- rlang::sym(covar)
      unique_covar_vals <- unique(pull(df, covar_sym))
      if (length(unique_covar_vals) <= 1) {
        cat(paste0("Covariate ", covar,
                   " is not variable in dataset. Removing...\n"))
        covar_vec <- covar_vec[!covar_vec %in% covar]
        covar_vec_mult <- covar_vec_mult[!grepl(covar, covar_vec_mult)]
      }
    }
  }
  covar_vec_mult2 <- str_trim(unlist(strsplit(covar_vec_mult, split = "[:]")))
  covariates_char <- paste(c(covar_vec_mult,
    setdiff(covar_vec, covar_vec_mult2)), collapse = " + ")
  return(covariates_char)
}

glm_run <- function(df, group_label = NULL, response,
                    reference = "Europe North West",
                    covariates = "", covariates_zero = "",
                    type = "overdispersed-count", total_var = NULL) {

  if (!is.null(group_label)) {

    group_label_sym <- rlang::sym(group_label)
    ## Reorder reference as factor level 1
    group_labels <- pull(df, {{group_label_sym}}) %>%
      unique()

    ## Mean center predictors?
    ## P value correction for number of covariates...
    num_tests <- length(group_labels) - 1

    df %<>%
      mutate({{group_label_sym}} := relevel(factor({{group_label_sym}}),
                                            ref = reference))

    ## Create + run the model
    if (covariates == "") {
      explan_char <- group_label
    } else {
      explan_char <- paste(group_label, covariates, sep = " + ")
    }

  } else {
    explan_char <- covariates
    num_tests <- 1
  }

  if (type == "overdispersed-count") {
    cat(paste("Running negative-binomial regression model for variable:",
      response, "\n"))
    cat(paste("Model: ", response, "~", explan_char, "\n"), sep = "")
    model <- glm.nb(as.formula(paste(response, "~", explan_char)),
                    data = df)

  } else if (type == "count") {
    cat(paste("Running poisson regression model for variable:",
      response, "\n"))
    cat(paste("Model: ", response, "~", explan_char, "\n"), sep = "")
    model <- glm(as.formula(paste(response, "~", explan_char)),
                 data = df, family = "poisson")

  } else if (type == "binary") {
    cat(paste("Running logistic regression model for variable:",
      response, "\n"))
    cat(paste("Model: ", response, "~", explan_char, "\n"), sep = "")
    response_sym <- rlang::sym(response)
    df %>% mutate({{response_sym}} := as.factor({{response_sym}})) %$%
      glm(as.formula(paste(response, "~", explan_char)),
          family = "binomial") -> model

  } else if (type == "proportion") {
    if (is.null(total_var)) {
      paste("Total variable needed for logistic regression (as proportion)!")
    } else {
      cat(paste("Running logistic regression model (proportion) for variable:",
      response, "/", total_var, "\n"))
      cat(paste("Model: ", "cbind(", response, ",", total_var, "-", response, ")", "~", explan_char, "\n"), sep = "")
      model <- glm(as.formula(paste0("cbind(", response, ",", total_var, "-", response, ") ~", explan_char)),
                   family = binomial(logit), data = df)
    }

  } else if (type == "zeroinfl") {
    cat(paste("Running zero-inflated regression (negbin for counts) model for variable:",
      response, "\n"))
    cat(paste("Model: ", response, "~", explan_char, "\n"), sep = "")
    model <- pscl::zeroinfl(as.formula(paste(response, "~", explan_char, "|", covariates_zero)),
                      data = df, dist = "negbin")
  }

  ## Model output
  if (type == "zeroinfl") {
    ci_model  <- data.frame(confint.default(model))
    colnames(ci_model) <- c("lower", "upper")
    rownames(ci_model) <- gsub("count_", "", rownames(ci_model))
    model_tidy <- summary(model)[[1]]$count %>%
      as.data.frame() %>%
      mutate(term = rownames(.)) %>%
      rename(p.value = `Pr(>|z|)`, estimate = Estimate) %>%
      filter(term != "Log(theta)")

  } else {
    ci_model  <- data.frame(confint.default(model))
    colnames(ci_model) <- c("lower", "upper")
    model_tidy <- model %>%
    tidy()
  }
  upper_vec <- pull(ci_model, upper)[rownames(ci_model) %in% model_tidy$term]
  lower_vec <- pull(ci_model, lower)[rownames(ci_model) %in% model_tidy$term]
  model_tidy %<>%
    mutate(upper = upper_vec, lower = lower_vec,
           p.adjust = p.value * num_tests,
           p.adjust = ifelse(p.adjust > 1, 1, p.adjust))

  if (!is.null(group_label)) {
    model_tidy_rm <- model_tidy %>%
      rename({{group_label_sym}} := term)

    model_out <- model_tidy_rm %>%
      mutate({{group_label_sym}} := gsub(group_label, "",
                                         pull(model_tidy_rm, {{group_label_sym}}))) %>%
      filter({{group_label_sym}} %in%  group_labels) %>%
      add_row(tibble_row({{group_label_sym}} := reference, estimate = 0))
  } else {
    model_out <- model_tidy
  }
  returnList <- list("model" = model,
                     "model_tidy" = model_out %>%
                        mutate(response = response))
  return(returnList)
}

glm_looper <- function(df, cat_var = NULL, group_label = NULL,
                       response_vec = "TIERED", reference = "Europe North West",
                       covariates = "", covariates_zero = 1,
                       type = "overdispersed-count", total_var_vec = NULL) {

  glm_list <- vector(mode = "list", length = length(response_vec))
  model_list <- list()
  counter <- 1
  for (response_index in 1:length(response_vec)) {
    if (!is.null(cat_var)) {
      cat_var_sym <- rlang::sym(cat_var)
    } else {
      df$dummy <- 0
      cat_var_sym <- rlang::sym("dummy")
    }
    cat_var_vec <- df %>% pull({{cat_var_sym}}) %>% unique()
    cat_glm_list <- vector(mode = "list", length = length(cat_var_vec))
    #### Loop through covariate (stratified if specified)
    for (cat_index in 1:length(cat_var_vec)) {
      #### Filter for specific disease
      df_filt <- df %>%
        filter({{cat_var_sym}} == cat_var_vec[cat_index])

      #### Per disease GLM
      covars_filt <- check_covariate_variation(df_filt,
                                               covariates)
      glm <- glm_run(df = df_filt,
                     group_label = group_label,
                     response = response_vec[response_index],
                     reference = reference,
                     covariates = covars_filt,
                     type = type,
                     covariates_zero = covariates_zero,
                     total_var = total_var_vec[response_index])
      cat_glm_list[[cat_index]] <- glm[[2]] %>%
        mutate(variable = cat_var_vec[cat_index])
      model_list[[counter]] <- glm[[1]]
      counter = counter + 1
    }
    glm_list[[response_index]] <- do.call(rbind, cat_glm_list)
  }
  glm_out <- do.call(rbind, glm_list) %>%
    mutate(p.adjust.groups = p.value * length(cat_var_vec),
           p.adjust.groups = ifelse(p.adjust.groups > 1, 1, p.adjust.groups))
  return(c(model_list, list(glm_out)))
}

###########
compare_coefs <- function(model1, model2) {
  model1 %<>% drop_na()
  model2 %<>% drop_na()
  # imports map_dbl() and pluck() functions from purrr library
  b1 <- model1 %>% pull(estimate)
  se1 = model1 %>% pull(std.error)
  upper1 = model1 %>% pull(upper)
  lower1 = model1 %>% pull(lower)

  b2 <- model2 %>% pull(estimate)
  se2 = model2 %>% pull(std.error)
  upper2 = model2 %>% pull(upper)
  lower2 = model2 %>% pull(lower)

  # Clogg et al. (1995) formula as cited by Ray Paternoster et al. (1998)
  b = b1 - b2
  s1 = se1^2
  s2 = se2^2
  sc = s1 + s2
  v = b / sqrt(sc)

  correlation = cor.test(b1, b2)
  out = data.frame(model1 = b1, upper951 = upper1, lower951 = lower1,
             model2 = b2, upper952 = upper2, lower952 = lower2,
             diff=b, zdiff=v, `p-value`= format(2*pnorm(-abs(v))))

  out$participant_genetic_category <- model1$participant_genetic_category
  return(list(correlation, out))
}

group_pheno_multi <- function(df, min_size = 10) {
  df %>%
    group_by(phenotype) %>%
    mutate(phenotype_narrow = ifelse(n() < min_size,
                                   "Other phenotypes",
                                   phenotype)) %>%
    ungroup()
}

P_disp <- function(x) {
  pr <- sum(residuals(x, type="pearson")^2)
  dispersion <- pr/x$df.residual
  cat("\n Pearson Chi2 = ", pr , "\n Dispersion = ", dispersion, "\n")
}

disperse_test <- function(x, data, response) {
  cat(paste0("Testing if response: ", response, " is overdispersed."))
  P_disp(glm_tier_poiss[[1]])
  response_sym <- rlang::sym(response)
  mu <- predict(x, type = "response")
  response_vec <- pull(data, {{response_sym}})
  z <- ((response_vec - mu)^2 - response_vec)/ (mu * sqrt(2))
  summary(zscore <- lm(z ~ 1))
}

###########################################################################
###########################################################################
###                                                                     ###
###                           DATA INGEST                               ###
###                                                                     ###
###########################################################################
###########################################################################

##### Meta data combined from various sources...
meta_combined_data <- read_csv(paste0(out_dir, "combined_data/combined_meta_data.csv.gz"))
nrow(meta_combined_data)
### 61,512 individuals

### Probands only
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
  tiering_frequency_filter(tiering_frequency_data, "ultra-rare-or-missing")

##### A small number of cPAVs are not ultra-rare...
tiering_probands_rare <- tiering_probands_filtered_data %>%
  tiering_frequency_filter(tiering_frequency_data, "rare")

#### Annotate tiering with Pathogenicity predictions...
predictions_data <- read_csv(paste0(out_dir, "combined_data/combined_predictions_data.csv.gz")) %>%
  rename(CADD_score = CADD_PHRED)

### Annotate with AlphaMissense
tiering_probands_ur_am <- tiering_probands_ur %>%
  pred_annot_adder(exit = NULL, annots = predictions_data, algo = "AM")

### Annotate with PrimateAI-3D (or PrimateAI if)
tiering_probands_ur_pai <- tiering_probands_ur %>%
  pred_annot_adder(exit = NULL, annots = predictions_data, algo = "PAI3D")

### Annotate with CADD (>30)
tiering_probands_ur_cadd <- tiering_probands_ur %>%
  pred_annot_adder(exit = NULL, annots = predictions_data, algo = "CADD")

##### Annotate exit data with tiering info...
exit_data <- read_csv(paste0(out_dir, "labkey_data/exit_data.csv.gz")) %>%
  add_varkey()
exit_tiers_filtered_data <- exit_data %>%
  acmg_tiers_adder(tiers = tiering_probands_filtered_data)

###########################################################################
###########################################################################
###                                                                     ###
###                           ANALYSIS RUN                              ###
###                                                                     ###
###########################################################################
###########################################################################

#### Configure response varianbles for running GLMs
analysis_data <- probands_meta_combined_data %>%
  relevel_covariates() %>%
  group_pheno_multi(min_size = 10) %>%
  configure_tiering_counts(tiering_probands_filtered_data,
                           merge = c("TIER1", "TIER2", "TIER3", "TIER0"),
                           merge_to = "cPAV") %>%
  configure_tiering_counts(tiering_probands_filtered_data,
                           merge = c("TIER1", "TIER2"),
                           merge_to = "GPcPAV") %>%
  configure_tiering_counts(tiering_probands_ur,
                           merge = c("TIER1", "TIER2", "TIER3", "TIER0"),
                           merge_to = "rare_cPAV") %>%
  configure_tiering_counts(tiering_probands_ur,
                           merge = c("TIER1", "TIER2"),
                           merge_to = "rare_GPcPAV")

##-------------------------------------------------------------------------
## cPAVs / GPcPAV - Ultra Rare
##-------------------------------------------------------------------------

### Total counts stratified by rarity classification.
### Uncommon cPAVs
counts_uncommon_cpavs <- tiering_probands_rare %>%
  group_by(varkey) %>%
  tally() %>%
  mutate(freq = "cPAVs popmax 1% > AF > 0.1%\nn = 86,570")

### Rare cPAVs
counts_ur_cpavs <- tiering_probands_ur %>%
  group_by(varkey) %>%
  tally() %>%
  mutate(freq = "cPAVs popmax AF < 0.1%\nn = 1,865,089")

#### Basic statistics...
analysis_data %>%
  cont_var_stats("participant_type", "cPAV")

# participant_type  mean median    sd upper_95 lower_95 variance   min   max
# Proband           172.    204  120.      392        8   14302.     0  1013

analysis_data %>%
  filter(penetrance == "complete") %>%
  cont_var_stats("family_group_type", "cPAV")

#  family_group_type     mean median    sd upper_95 lower_95 variance   min   max
#  <chr>                <dbl>  <dbl> <dbl>    <dbl>    <dbl>    <dbl> <dbl> <dbl>
#1 Duo with Mother or … 122.   113    39.7     217.    77       1574.     0   388
#2 Duo with other Biol… 167.   139    80.1     361.    63       6411.     0   605
#3 Families with more …  25.2   12    40.9     118.     1       1672.     0   404
#4 Singleton            273.   257    68.5     431    192       4688.     0  1013
#6 Trio with Mother or…  73.0   69.5  42.3     172.     9.05    1788.     0   424
#5 Trio with Mother an…  23.3   18    21.6      92      6        466.     0   422
#7 Trio with other Bio… 123.    83    90.1     342.    27.5     8119.     0   485

analysis_data %>%
  filter(penetrance == "incomplete") %>%
  cont_var_stats("participant_type", "cPAV")

# 1 Proband           216.    212  74.2      370       75    5506.    16   760

#### cPAV / GPcPAV counts GLM (poisson)
glm_tier_poiss <- glm_looper(
  df = analysis_data,
  group_label = "participant_genetic_category",
  response_vec = c("cPAV", "GPcPAV"),
  type = "count",
  reference = "Europe North West",
  covariates = ""
)

disperse_test(glm_tier_poiss[[1]], analysis_data, "cPAV")
## p < 2-16

disperse_test(glm_tier_poiss[[2]], analysis_data, "GPcPAV")
## p < 2-16

#### Penetrance is only relevant if family_group_type != Singleton (singletons encoded as dummy variable)
#### Penetrance/fgt as an interaction effect only.
covariates_vec <- "family_group_type + penetrance + family_group_type:penetrance + age + c_roh_mb + n_roh + only_proband_affected + karyotype + phenotype_narrow + date_months + unique_genes_panels_applied_mb + panels_applied_count"

#### Quality covariates
qc_covariates <- "samtools_error_rate + samtools_reads_mapped_percentage + illumina_autosome_mean_coverage"

### Add handling GMC trust as a covariate moving forward.
covariates_yield <- "handling_gmc_trust"
covariates_no_phenotype <- "family_group_type + penetrance + age + c_roh_mb + n_roh + only_proband_affected + karyotype + date_months + unique_genes_panels_applied_mb + panels_applied_count"

#### Run full model with covariates...
glm_tier_nb_cov <- glm_looper(
  df = analysis_data,
  group_label = "participant_genetic_category",
  response_vec = c("rare_cPAV", "rare_GPcPAV", "cPAV", "GPcPAV"),
  type = "overdispersed-count",
  reference = "Europe North West",
  covariates = paste(covariates_vec, qc_covariates, sep = "+")
)

#### Are the effects of GIA group on cPAVs + GPcPAVs essentially identical?
cpav_vs_gpcpav <- compare_coefs(glm_tier_nb_cov[[length(glm_tier_nb_cov)]] %>%
                                  filter(response == "cPAV"),
                                glm_tier_nb_cov[[length(glm_tier_nb_cov)]] %>%
                                  filter(response == "GPcPAV"))
min(as.numeric(cpav_vs_gpcpav[[2]]$p.value))
cpav_vs_gpcpav[[1]]
### min p = 0.03595557

glm_tier_nb_cov_miss <- glm_looper(
  df = analysis_data %>%
    mutate(missing_pavs_count = round(missing_pavs_count / 10)),
  response_vec = c("rare_cPAV", "rare_GPcPAV", "GPcPAV", "cPAV"),
  type = "overdispersed-count",
  covariates = paste("missing_pavs_count", covariates_vec, qc_covariates, sep = "+")
)

##-------------------------------------------------------------------------
## PAVs missing from gnomAD
##-------------------------------------------------------------------------

#### Missing PAVs from gnomAD
analysis_data_plots %>%
  cont_var_stats("participant_genetic_category", "missing_pavs_count") %>%
  as.data.frame()

glm_pavs_miss <- glm_looper(
  df = analysis_data,
  group_label = "participant_genetic_category",
  response_vec = "missing_pavs_count",
  type = "overdispersed-count",
  reference = "Europe North West",
  covariates = qc_covariates
)

##-------------------------------------------------------------------------
## Ultra-rare cPAVs / cPAVs applied panels - Pathogenicity Predictions
##-------------------------------------------------------------------------

analysis_data %<>%
  configure_tiering_counts(tiering_probands_ur_am,
                           merge = c("TIER3", "TIER2", "TIER1", "TIER0"),
                           merge_to = "rare_cPAV", source = "AM") %>%
  mutate(rare_cPAV_AM_nDeleterious = rare_cPAV_AM - rare_cPAV_AM_Deleterious) %>%
  configure_tiering_counts(tiering_probands_ur_am,
                           merge = c("TIER2", "TIER1"),
                           merge_to = "rare_GPcPAV", source = "AM") %>%
  mutate(rare_GPcPAV_AM_nDeleterious = rare_GPcPAV_AM - rare_GPcPAV_AM_Deleterious) %>%
  configure_tiering_counts(tiering_probands_ur_pai,
                           merge = c("TIER3", "TIER2", "TIER1", "TIER0"),
                           merge_to = "rare_cPAV", source = "PAI3D") %>%
  mutate(rare_cPAV_PAI3D_nDeleterious = rare_cPAV_PAI3D - rare_cPAV_PAI3D_Deleterious) %>%
  configure_tiering_counts(tiering_probands_ur_pai,
                           merge = c("TIER2", "TIER1"),
                           merge_to = "rare_GPcPAV", source = "PAI3D") %>%
  mutate(rare_GPcPAV_PAI3D_nDeleterious = rare_GPcPAV_PAI3D - rare_GPcPAV_PAI3D_Deleterious) %>%
  configure_tiering_counts(tiering_probands_ur_cadd,
                           merge = c("TIER3", "TIER2", "TIER1", "TIER0"),
                           merge_to = "rare_cPAV", source = "CADD") %>%
  mutate(rare_cPAV_CADD_nDeleterious = rare_cPAV_CADD - rare_cPAV_CADD_Deleterious) %>%
  configure_tiering_counts(tiering_probands_ur_cadd,
                           merge = c("TIER2", "TIER1"),
                           merge_to = "rare_GPcPAV", source = "CADD") %>%
  mutate(rare_GPcPAV_CADD_nDeleterious = rare_GPcPAV_CADD - rare_GPcPAV_CADD_Deleterious)

##### We want to know if missing variants from gnomAD is correlated with numbers of cPAVs
##### of each pathogenicity class?
##### Across multiple distinct methods of predicting pathogenicity.

#### Run full model with covariates...
glm_tier_nb_cov_delet <- glm_looper(
  df = analysis_data,
  group_label = "participant_genetic_category",
  response_vec = c(
    "rare_cPAV_CADD_Deleterious",
    "rare_cPAV_PAI3D_Deleterious",
    "rare_cPAV_AM_Deleterious",
    "rare_GPcPAV_CADD_Deleterious",
    "rare_GPcPAV_PAI3D_Deleterious",
    "rare_GPcPAV_AM_Deleterious",
    "rare_cPAV_CADD_nDeleterious",
    "rare_cPAV_PAI3D_nDeleterious",
    "rare_cPAV_AM_nDeleterious",
    "rare_GPcPAV_CADD_nDeleterious",
    "rare_GPcPAV_PAI3D_nDeleterious",
    "rare_GPcPAV_AM_nDeleterious"
  ),
  type = "overdispersed-count",
  reference = "Europe North West",
  covariates = paste(covariates_vec, qc_covariates, sep = "+")
)

##### Is there variation in the proportion (%) of deleterious cPAVs
##### across ancestry groups?
glm_tier_patho_ukbb <- glm_looper(
  df = analysis_data,
  group_label = "participant_genetic_category",
  response_vec = c(
    "rare_cPAV_AM_Deleterious",
    "rare_GPcPAV_AM_Deleterious",
    "rare_cPAV_PAI3D_Deleterious",
    "rare_GPcPAV_PAI3D_Deleterious",
    "rare_cPAV_CADD_Deleterious",
    "rare_GPcPAV_CADD_Deleterious"
  ),
  type = "proportion",
  reference = "Europe North West",
  total_var_vec = c(
    "rare_cPAV_AM",
    "rare_GPcPAV_AM",
    "rare_cPAV_PAI3D",
    "rare_GPcPAV_PAI3D",
    "rare_cPAV_CADD",
    "rare_GPcPAV_CADD"
  ),
  covariates = paste(covariates_vec, qc_covariates, sep = "+")
)

##-------------------------------------------------
## cPAVs / cPAVs applied panels - Diagnostics
##-------------------------------------------------

### Update the analysis dataframe
analysis_data %<>%
  configure_tiering_counts(exit_tiers_filtered_data,
                           merge = c("TIER2", "TIER1"),
                           merge_to = "GPcPAV",
                           source = "ACMG") %>%
  configure_diagnosed(patho_variants = exit_tiers_filtered_data,
                      with = "GPcPAV_ACMG_Pathogenic",
                      diagnosis_definition = c(
                        "yes",
                        "partially"
                      )) %>%
  configure_tiering_counts(exit_tiers_filtered_data,
                           merge = c("TIER2", "TIER1", "TIER3", "TIER0"),
                           merge_to = "cPAV",
                           source = "ACMG") %>%
  configure_diagnosed(patho_variants = exit_tiers_filtered_data,
                      with = "cPAV_ACMG_Pathogenic",
                      diagnosis_definition = c(
                        "yes",
                        "partially"
                      )) %>%
  mutate(GPcPAV_VUS_or_noclass = GPcPAV - (GPcPAV_ACMG_Benign + GPcPAV_ACMG_Pathogenic))

##--------------------------------------------------
## cPAVs applied panels - Positive Predictive Value
##--------------------------------------------------

### Remove the interaction effect between family_group_type + penetrance.
covariates_vec <- "family_group_type + penetrance + age + c_roh_mb + n_roh + only_proband_affected + karyotype + phenotype_narrow + date_months + unique_genes_panels_applied_mb + panels_applied_count"

### Positive predictive value (PPV) main model
glm_tiered_ppv <- glm_looper(
  df = analysis_data,
  group_label = "participant_genetic_category",
  response_vec = c(
    "GPcPAV_ACMG_Pathogenic",
    "cPAV_ACMG_Pathogenic"
  ),
  reference = "Europe North West",
  type = "proportion",
  total_var_vec = c(
    "GPcPAV",
    "cPAV"
  ),
  covariates = paste(covariates_vec, covariates_yield,
                     qc_covariates, sep = " + ")
)

##-------------------------------------------------
## cPAVs / cPAVs applied panels - Yield
##-------------------------------------------------

##### Check whether there is variation in...
##### Diagnostic yield (overall)
##### Diagnostic yield from GPcPAVs

glm_tiered_yield <- glm_looper(
  df = analysis_data,
  group_label = "participant_genetic_category",
  response_vec = c(
    "DIAGNOSED_GPcPAV_ACMG_Pathogenic",
    "DIAGNOSED_cPAV_ACMG_Pathogenic"
  ),
  reference = "Europe North West",
  type = "binary",
  covariates = paste(covariates_yield, qc_covariates, covariates_vec, sep = " + ")
)

DescTools::PseudoR2(glm_tiered_yield[[1]], which = "Nagelkerke")
### Model explains 0.1874801 of the variation in diagnostic yield.

##### Analysis of deviance (Type II)
anova_tiered_yield <- car::Anova(glm_tiered_yield[[1]], type = "II")

### Full model (for plotting)
full_model_tiered_yield <- glm_tiered_yield[[1]] %>%
  tidy()

ci_model <- as.data.frame(confint.default(glm_tiered_yield[[1]]))
lower_vec <- pull(ci_model, `2.5 %`)[rownames(ci_model) %in% full_model_tiered_yield$term]
upper_vec <- pull(ci_model, `97.5 %`)[rownames(ci_model) %in% full_model_tiered_yield$term]

### Without phenotype as a covariate...
glm_tiered_yield_nopheno <- glm_looper(
  df = analysis_data,
  group_label = "participant_genetic_category",
  response_vec = c(
    "DIAGNOSED_GPcPAV_ACMG_Pathogenic"
  ),
  reference = "Europe North West",
  type = "binary",
  covariates = paste(covariates_yield, qc_covariates, covariates_no_phenotype, sep = " + ")
)

DescTools::PseudoR2(glm_tiered_yield_nopheno[[1]], which = "Nagelkerke")
### Model (without disease phenotype) explains 0.04278795 of the variation in diagnostic yield.

##### Instead, utilise ancestry groups provided by gnomAD
##### Only groups with n > 1000 probands.
glm_tiered_yield_gnomad <- glm_looper(
  df = analysis_data %>%
    ## Only keep largest ancestry groups
    filter(participant_gnomad_category %in% c("nfe", "afr", "sas", "eas")) %>%
    group_pheno_multi(min_size = 10),
  group_label = "participant_gnomad_category",
  response_vec = c(
    "DIAGNOSED_cPAV_ACMG_Pathogenic",
    "DIAGNOSED_GPcPAV_ACMG_Pathogenic"
  ),
  reference = "nfe",
  type = "binary",
  covariates = paste(covariates_vec, qc_covariates,
                     covariates_yield, sep = " + ")
)

##--------------------------------------------------
## cPAVs applied panels - VUS
##--------------------------------------------------

#### Run full model with covariates...
glm_vus_nb_cov <- glm_looper(
  df = analysis_data,
  cat_var = "DIAGNOSED",
  group_label = "participant_genetic_category",
  response_vec = "GPcPAV_VUS_or_noclass",
  type = "overdispersed-count",
  reference = "Europe North West",
  covariates = paste(covariates_vec, qc_covariates, covariates_yield, sep = "+")
)

#########################################################################
#########################################################################
###                                                                   ###
###                          WRITE DATA                               ###
###                                                                   ###
#########################################################################
#########################################################################

## Analysis data
data_dir <- paste0(analysis_dir, "glm_analysis/")

## Create directory.
dir.create(file.path(data_dir))
cat(paste0("Writing data to: ", data_dir, "\n"))

### Write final analysis dataframe
fwrite(analysis_data, paste0(data_dir, "analysis_data.csv.gz"))

### Write this to tables.
rbind(counts_ur_cpavs, counts_uncommon_cpavs) %>%
write.csv(paste0(data_dir, "cpavs_counts_rd.csv"),
          row.names = FALSE)

### Write output files...
glm_tier_nb_cov[[length(glm_tier_nb_cov)]] %>%
  mutate(risk_ratio = exp(estimate)) %>%
  write.csv(paste0(data_dir, "ukbb_cpavs_glm.csv"),
            row.names = FALSE)

### Predict n cPAVs with n PAVs missing from gnomAD
glm_tier_nb_cov_miss[[length(glm_tier_nb_cov_miss)]] %>%
  mutate(risk_ratio = exp(estimate)) %>%
  filter(term == "missing_pavs_count") %>%
  write.csv(paste0(data_dir, "missing_gnomad_cpavs_glm.csv"),
            row.names = FALSE)

## GLM PAVs missing from gnomAD GLM
glm_pavs_miss[[length(glm_pavs_miss)]] %>%
  mutate(risk_ratio = exp(estimate)) %>%
  write.csv(paste0(data_dir, "ukbb_missing_glm.csv"),
            row.names = FALSE)

## Pathogenicity counts GLM (deleterious / non-deleterious)
glm_tier_nb_cov_delet[[length(glm_tier_nb_cov_delet)]] %>%
  mutate(risk_ratio = exp(estimate)) %>%
  write.csv(paste0(data_dir, "ukbb_patho_count_glm.csv"),
            row.names = FALSE)

## Pathogenicity prediction GLM.
glm_tier_patho_ukbb[[length(glm_tier_patho_ukbb)]] %>%
  mutate(risk_ratio = exp(estimate)) %>%
  write.csv(paste0(data_dir, "ukbb_patho_prop_glm.csv"),
            row.names = FALSE)

## Diagnostic yield GLM.
glm_tiered_yield[[length(glm_tiered_yield)]] %>%
  mutate(risk_ratio = exp(estimate)) %>%
  write.csv(paste0(data_dir, "ukbb_diagnostic_yield_glm.csv"),
            row.names = FALSE)

## Diagnostic yield ANOVA
write.csv(anova_tiered_yield %>% as.data.frame() %>% mutate(covar = rownames(.)),
          paste0(data_dir, "ukbb_diagnostic_yield_anova.csv"),
          row.names = FALSE)

## Full model tiered yield
write.csv(full_model_tiered_yield,
          paste0(data_dir, "ukbb_diagnostic_yield_full_model_glm.csv"),
          row.names = FALSE)

## Diagnostic yield GLM gnomAD
glm_tiered_yield_gnomad[[length(glm_tiered_yield_gnomad)]] %>%
  mutate(risk_ratio = exp(estimate)) %>%
  write.csv(paste0(data_dir, "gnomad_diagnostic_yield_glm.csv"),
            row.names = FALSE)

## PPV GLM.
glm_tiered_ppv[[length(glm_tiered_ppv)]] %>%
  mutate(risk_ratio = exp(estimate)) %>%
  write.csv(paste0(data_dir, "ukbb_ppv_glm.csv"),
            row.names = FALSE)

### VUS / unclassified GLM
glm_vus_nb_cov[[length(glm_vus_nb_cov)]] %>%
  mutate(risk_ratio = exp(estimate)) %>%
  write.csv(paste0(data_dir, "ukbb_vus_glm.csv"),
            row.names = FALSE)

###########################################################################
###########################################################################
###                                                                     ###
###                         END OF SCRIPT                               ###
###                                                                     ###
###########################################################################
###########################################################################
