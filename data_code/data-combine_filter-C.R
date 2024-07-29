#!/usr/bin/env Rscript

##-------------------------------------------------------------------------
## Configs + packages

source("../config/config.R")

library(tidyverse)
library(data.table)
library(magrittr)

### PrimateAI-3D access argument
args <- commandArgs(trailingOnly = TRUE)

### Set no as the default if no arguments are selected.
if (length(args) != 2) {
  args[1] <- "NO_PAI3D_ACCESS"
  args[2] <- "NO_COVID_ACCESS"
}

###########################################################################
###########################################################################
###                                                                     ###
###                  JOINING AND FILTERING FUNCTIONS                    ###
###                              ALL                                    ###
###                                                                     ###
###########################################################################
###########################################################################

### Static
tiers_new_freq_join <- function(freq_data, new_freq_data) {
  freq_data %>%
    left_join(new_freq_data, by = "varkey") %>%
    mutate_at(vars(contains("_")), ~replace_na(as.numeric(.), 0))
}

count_common_pops <- function(data, AF = 0.01, pop_id = "_") {
  af_filter_pops <- lapply(apply(data %>%
                         dplyr::select(contains(pop_id)), 1,
                         function(x) which(x >= AF)),
                         names)

  pops_unique <- colnames(data %>%
                        dplyr::select(contains("_")))
  
  pop_common_counts <- data.frame(participant_genetic_category = NA_character_, num_instances = NA_integer_)
  for (index in 1:length(pops_unique)) {
    pop_common_counts[index, 1] <- pops_unique[index]
    pop_common_counts[index, 2] <- length(grep(pops_unique[index], af_filter_pops))
    pop_common_counts[index, 3] <- length(grep(pops_unique[index], af_filter_pops[sapply(af_filter_pops, function(i) length(i) == 1)]))
  }
  return(pop_common_counts)
}

tiers_harmonize_afs_filter <- function(tiers, freq_data) {
  ## Filter freq data
  freq_data <- freq_data %>% filter(varkey %in% tiers$varkey)

  ### Common tiered variants
  tiered_variants_common_big <- freq_data %>%
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
                       "gel_gl_6628_af")), any_vars(. >= 0.01))

  tiered_variants_common_small <- freq_data %>%
      filter_at(vars(c("AF_oth_gnomad_v3",
                       "AF_asj_gnomad_v3")), any_vars(. >= 0.02))

  tiered_variants_common <- rbind(tiered_variants_common_big,
                                  tiered_variants_common_small)
  cat(paste0(length(unique(tiered_variants_common$varkey)), " common (AF > 1%) variants\n"))
  
  ### Rare tiered variants
  tiered_variants_rare_big <- freq_data %>%
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
                       "gel_gl_6628_af")), any_vars(. >= 0.001))

  tiered_variants_rare_med <- freq_data %>%
      filter_at(vars(c("AF_oth_gnomad_v3",
                       "AF_asj_gnomad_v3")), any_vars(. >= 0.002))

  tiered_variants_rare <- rbind(tiered_variants_rare_big,
                                tiered_variants_rare_med)
  cat(paste0(length(unique(tiered_variants_rare$varkey)), " uncommon (AF > 0.1%) variants\n"))

  nvar <- length(unique(tiers$varkey))
  tiers_ur_mono_filt <- tiers %>%
    filter(!(mode_of_inheritance %in% c("monoallelic", "monoallelic_not_imprinted",
                                        "monoallelic_paternally_imprinted",
                                        "monoallelic_maternally_imprinted",
                                        "xlinked_monoallelic", "mitochondrial") &
            segregation_pattern != "deNovo" &
            varkey %in% tiered_variants_rare$varkey))

  nvar_ur_mono <- length(unique(tiers_ur_mono_filt$varkey))
  cat(paste0(nvar - nvar_ur_mono, " variants removed after not passing dominant AF filter!\n"))
  # Which pops were these rare variants found in?
  tiered_variants_rare_filt <- tiered_variants_rare %>%
    filter(!varkey %in% tiers_ur_mono_filt$varkey) 
  # Count common pops
  pc_rare_filt <- tiered_variants_rare_filt %>%
    count_common_pops(AF = 0.001)
  print(pc_rare_filt)

  tiers_full_filt <- tiers_ur_mono_filt %>%
    filter(!varkey %in% tiered_variants_common$varkey)
  nvar_ur_rec <- length(unique(tiers_full_filt$varkey))

  cat(paste0(nvar_ur_mono - nvar_ur_rec, " variants removed after not passing recessive AF filter!\n"))
  # Which pops were common variants found in?
  
  pc_common_filt <- tiered_variants_common %>%
    filter(!varkey %in% tiers_full_filt$varkey) %>%
    count_common_pops()
  print(pc_common_filt)
  return(tiers_full_filt)
}

## Tiered variants sometimes occur multiple times (e.g. mulitple MOI, get only distinct cPAVs per ind, taking the highest tier where relevant)
tiers_replicates_class_filter <- function(tiers) {
  tiers %>%
    dplyr::select(participant_id, varkey, varid, tier, genotype) %>%
    distinct() %>%
    mutate(n = 1) %>%
    pivot_wider(names_from = tier, values_from = n) %>%
    mutate(tier = case_when(TIER1 == 1 ~ "TIER1",
                            TIER2 == 1 & is.na(TIER1) ~ "TIER2",
                            TIER3 == 1 & is.na(TIER1) & is.na(TIER2) ~ "TIER3")) %>%
    dplyr::select(participant_id, varkey, varid, tier, genotype)
}

retiers_replicates_class_filter <- function(tiers, filter_tier0 = FALSE) {
  tiers %<>%
    dplyr::select(-tier) %>%
    rename(tier = retier) %>%
    dplyr::select(participant_id, varkey, varid, tier, genotype) %>%
    distinct() %>%
    mutate(n = 1) %>%
    pivot_wider(names_from = tier, values_from = n) %>%
    mutate(tier = case_when(TIER2 == 1 ~ "TIER2",
                            TIER3 == 1 & is.na(TIER2) ~ "TIER3",
                            TIER0 == 1 & is.na(TIER2) & is.na(TIER3) ~ "TIER0")) %>%
    dplyr::select(participant_id, varkey, varid, tier, genotype)
  if (filter_tier0) {
    tiers %<>%
      filter(tier != "TIER0")
  }
  return(tiers)
}

retier_variants <- function(tiering_data, panels_applied, karyotypes, so_terms,
                            panel_genes, confidence_score = 3, keep_old_tier12 = FALSE) {

  ## Distinct tiers for variants
  tiering_data %>%
    tiers_replicates_class_filter() %>%
    mutate(tier = ifelse(tier == "TIER1", "TIER2", tier)) %>%
    group_by(tier) -> tiered_variants
  
  cat("Total count of variants (prior to retiering):\n")
  tiered_variants %>%
    tally() -> count_tiers
  print(count_tiers)
   
  ## Simplify tiering data to only include relevant SO terms + necessary annotations
  ## e.g. Consequence associated with specific HGNC feature
  tiering_data_simple <- tiering_data %>%
    filter(so_term %in% so_terms | is.na(so_term)) %>%
    distinct() %>%
    left_join(panels_applied, by = "participant_id",
              relationship = "many-to-many")

  ## Simplify PanelApp panels to only include those with specific confidence score (3 = "Green Genes")
  panel_genes_simple <- panel_genes %>%
    filter(confidence >= confidence_score) %>%
    dplyr::select(genomic_feature_hgnc, panel_name,
                  mode_of_inheritance) %>%
    rename(mode_of_inheritance_panel_app = mode_of_inheritance) %>%
    distinct() %>%
    mutate(mode_of_inheritance_panel_app = ifelse(is.na(mode_of_inheritance_panel_app),
                                                  "",
                                                  mode_of_inheritance_panel_app))

  ## Join data
  tiering_data_panels_simple <- tiering_data_simple %>%
    left_join(panel_genes_simple,
              by = c("genomic_feature_hgnc", "panel_name")) %>%
    left_join(karyotypes, by = "participant_id") %>%
    mutate(mode_of_inheritance = ifelse(segregation_pattern == "deNovo",
                                        "deNovo",
                                        mode_of_inheritance))

  ## Match MOI with PanelApp gene panel MOI (from Rare Disease Analysis Guide)
  ## TIER2 if MOI does match, TIER0 variants if MOI does not match, TIER3 if gene not in panel
  out_data <- tiering_data_panels_simple %>%
    mutate(retier = case_when(
      mode_of_inheritance == "biallelic" &
      mode_of_inheritance_panel_app %in% c(
        "BIALLELIC, autosomal or pseudoautosomal",
        "BOTH monoallelic and biallelic (but BIALLELIC mutations cause a more SEVERE disease form), autosomal or pseudoautosomal",
        "BOTH monoallelic and biallelic, autosomal or pseudoautosomal",
        "", "Other - please specifiy in evaluation comments"
      ) ~ "TIER2",
      mode_of_inheritance == "monoallelic" &
      mode_of_inheritance_panel_app %in% c(
        "MONOALLELIC, autosomal or pseudoautosomal, imprinted status unknown",
        "BOTH monoallelic and biallelic (but BIALLELIC mutations cause a more SEVERE disease form), autosomal or pseudoautosomal",
        "BOTH monoallelic and biallelic, autosomal or pseudoautosomal",
        "MONOALLELIC, autosomal or pseudoautosomal, maternally imprinted (paternal allele expressed)",
        "MONOALLELIC, autosomal or pseudoautosomal, NOT imprinted",
        "MONOALLELIC, autosomal or pseudoautosomal, paternally imprinted (maternal allele expressed",
        "X-LINKED: hemizygous mutation in males, biallelic mutations in females",
        "X-LINKED: hemizygous mutation in males, monoallelic mutations in females may cause disease (may be less severe, later onset than males)",
        "", "Other - please specifiy in evaluation comments"
      ) ~ "TIER2",
       mode_of_inheritance == "deNovo" &
       mode_of_inheritance_panel_app %in% c(
        "MONOALLELIC, autosomal or pseudoautosomal, imprinted status unknown",
        "BOTH monoallelic and biallelic (but BIALLELIC mutations cause a more SEVERE disease form), autosomal or pseudoautosomal",
        "BOTH monoallelic and biallelic, autosomal or pseudoautosomal",
        "MONOALLELIC, autosomal or pseudoautosomal, maternally imprinted (paternal allele expressed)",
        "MONOALLELIC, autosomal or pseudoautosomal, NOT imprinted",
        "MONOALLELIC, autosomal or pseudoautosomal, paternally imprinted (maternal allele expressed",
        "X-LINKED: hemizygous mutation in males, biallelic mutations in females",
        "X-LINKED: hemizygous mutation in males, monoallelic mutations in females may cause disease (may be less severe, later onset than males)",
        "", "Other - please specifiy in evaluation comments"
      ) ~ "TIER2",
      mode_of_inheritance == "monoallelic_not_imprinted" &
      mode_of_inheritance_panel_app %in% c(
        "MONOALLELIC, autosomal or pseudoautosomal, imprinted status unknown",
        "BOTH monoallelic and biallelic (but BIALLELIC mutations cause a more SEVERE disease form), autosomal or pseudoautosomal",
        "BOTH monoallelic and biallelic, autosomal or pseudoautosomal",
        "MONOALLELIC, autosomal or pseudoautosomal, maternally imprinted (paternal allele expressed)",
        "MONOALLELIC, autosomal or pseudoautosomal, NOT imprinted",
        "MONOALLELIC, autosomal or pseudoautosomal, paternally imprinted (maternal allele expressed)",
        "X-LINKED: hemizygous mutation in males, biallelic mutations in females",
        "X-LINKED: hemizygous mutation in males, monoallelic mutations in females may cause disease (may be less severe, later onset than males)",
        "", "Other - please specifiy in evaluation comments"
      ) ~ "TIER2",
      mode_of_inheritance == "monoallelic_maternally_imprinted" &
      mode_of_inheritance_panel_app %in% c(
        "MONOALLELIC, autosomal or pseudoautosomal, maternally imprinted (paternal allele expressed)",
        "", "Other - please specifiy in evaluation comments"
      ) ~ "TIER2",
      mode_of_inheritance == "monoallelic_paternally_imprinted" &
      mode_of_inheritance_panel_app %in% c(
        "MONOALLELIC, autosomal or pseudoautosomal, paternally imprinted (maternal allele expressed)",
        "", "Other - please specifiy in evaluation comments"
      ) ~ "TIER2",
      mode_of_inheritance == "mitochondrial" &
      mode_of_inheritance_panel_app %in% c(
        "MITOCHONDRIAL",
        "", "Other - please specifiy in evaluation comments"
      ) ~ "TIER2",
      mode_of_inheritance == "xlinked_monoallelic" &
      karyotype %in% c("XX", "XXX") &
      mode_of_inheritance_panel_app %in% c(
        "X-LINKED: hemizygous mutation in males, monoallelic mutations in females may cause disease (may be less severe, later onset than males)",
        "", "Other - please specifiy in evaluation comments"
      ) ~ "TIER2",
      mode_of_inheritance == "xlinked_biallelic" &
      karyotype %in% c("XX", "XXX") &
      mode_of_inheritance_panel_app %in% c(
        "X-LINKED: hemizygous mutation in males, biallelic mutations in females",
        "", "Other - please specifiy in evaluation comments"
      ) ~ "TIER2",
      is.na(mode_of_inheritance_panel_app) ~ "TIER3",
      TRUE ~ "TIER0")) %>%
  distinct()

  if (keep_old_tier12) {
    out_data %<>%
      mutate(retier = ifelse(tier %in% c("TIER1", "TIER2") & retier %in% c("TIER3", "TIER0"),
                             "TIER2",
                             retier))
  }

  ## Distinct tiers for variants (after retiering TIER3)
  out_data %>%
    retiers_replicates_class_filter() %>%
    mutate(tier = ifelse(tier == "TIER1", "TIER2", tier)) %>%
    group_by(tier) -> retiered_variants
  
  cat("Total count of variants (after retiering):\n")
  retiered_variants %>%
    tally() -> count_retiers
  print(count_retiers)

  return(out_data)
}

annot_filter <- function(annots, algo) {
  annot_score <- rlang::sym(paste0(algo, "_score"))
  annot_pred <- rlang::sym(paste0(algo, "_pred"))
 
  #### In cases where multiple variants map to canonical transcript...
  #### Take the most most pathogenic variant (max for AM/pp2, min for SIFT)
  annots %<>%
    mutate(varkey = paste0(CHROM, ":",
                           POS, "_", REF, "_", ALT),
           PrimateAI_score = 0) %>%
    dplyr::select(varkey, {{annot_score}}, {{annot_pred}}) %>%
    filter({{annot_score}} != "." & {{annot_pred}} != ".") %>%
    mutate({{annot_score}} := as.numeric({{annot_score}})) %>%
    distinct()
  if (algo %in% c("AM", "Polyphen2_HVAR")) {
    annots %<>%
      group_by(varkey) %>%
      filter({{annot_score}} == max({{annot_score}})) %>%
      ungroup()
  } else {
    annots %<>%
      group_by(varkey) %>%
      filter({{annot_score}} == min({{annot_score}})) %>%
      ungroup()
  }
  annots %>%
    dplyr::select(-{{annot_score}}) %>%
    mutate({{annot_pred}} := case_when(grepl("likely_benign|T|B", {{annot_pred}}) ~ "Tolerated",
                                       grepl("likely_pathogenic|D", {{annot_pred}}) ~ "Deleterious",
                                       grepl("ambiguous|P", {{annot_pred}}) ~ "Ambiguous"))
}

cadd_filter <- function(annots, score_cutoff = 30) {

  annots %>%
    mutate(varkey = paste0(CHROM, ":",
                           POS, "_", REF, "_", ALT)) %>%
    dplyr::select(varkey, CADD_PHRED) %>%
    filter(CADD_PHRED != ".") %>%
    mutate(CADD_PHRED = as.numeric(CADD_PHRED)) %>%
    distinct() %>%
    group_by(varkey) %>%
    filter(CADD_PHRED == max(CADD_PHRED)) %>%
    ungroup() %>%
    mutate(CADD_pred = ifelse(CADD_PHRED >= score_cutoff,
                              "Deleterious", "Tolerated")) %>%
    dplyr::select(varkey, CADD_pred)
}

pai_filter <- function(annots, score_cutoff = 0.821) {
  ## Pathogenicity
  annots %>%
    mutate(varkey = paste0(chromosome, ":",
                           position, "_", reference, "_", alternate)) %>%
    dplyr::select(varkey, PAI3D_score) %>%
    mutate(PAI3D_score = as.numeric(PAI3D_score)) %>%
    distinct() %>%
    group_by(varkey) %>%
    filter(PAI3D_score == max(PAI3D_score)) %>%
    ungroup() %>%
    mutate(PAI3D_pred = ifelse(PAI3D_score >= score_cutoff,
                               "Deleterious", "Tolerated")) %>%
    dplyr::select(varkey, PAI3D_pred)
}

configure_age_date <- function(meta) {
  elapsed_months <- function(end_date, start_date) {
    ed <- as.POSIXlt(end_date)
    sd <- as.POSIXlt(start_date)
    12 * (ed$year - sd$year) + (ed$mon - sd$mon)
  }
  meta %>%
    mutate(age = as.numeric(format(delivery_date, format = "%Y")) - year_of_birth,
           date_months = as.numeric(elapsed_months(delivery_date, min(meta$delivery_date))))
}

configure_diversity_statistics <- function(meta, roh, missing_gnomad_pavs) {
  meta %>%
    left_join(missing_gnomad_pavs %>%
                mutate(missing_pavs_count = HOM_ALT_CT + HET_CT,
                       missing_pavs_hets_count = HET_CT,
                       missing_pavs_homs_count = HOM_ALT_CT) %>%
                dplyr::select(platekey,
                              missing_pavs_count,
                              missing_pavs_hets_count,
                              missing_pavs_homs_count),
              by = "platekey") %>%
    left_join(roh, by = "participant_id") %>%
    mutate(all_variants_count = illumina_indels_all + illumina_snvs_all,
           #### Het/Hom Ratio is separated by INDELS/SNVs - find Het/Hom ratio for all variants
           all_variants_het_hom_ratio = ((illumina_indel_het_hom_ratio * illumina_indels_all) +
                                        (illumina_snv_het_hom_ratio * illumina_snvs_all)) /
                                        (illumina_indels_all + illumina_snvs_all),
           c_roh_mb = c_roh / 1000000,
           a_roh_mb = a_roh / 1000000,
           samtools_reads_mapped_percentage = samtools_reads_mapped /
                                              (samtools_reads_mapped + samtools_reads_unmapped)) %>%
    dplyr::select(-c(samtools_reads_mapped, samtools_reads_unmapped,
                     illumina_snvs_all, illumina_indels_all, illumina_indel_het_hom_ratio,
                     illumina_snv_het_hom_ratio, file_path, c_roh, a_roh))
}

configure_phenotype <- function(meta, phenotypes, phenotypes_nejm, panel_stats) {
  meta %>%
  left_join(phenotypes %>%
              group_by(participant_id) %>%
              mutate(phenotypes_count = n()) %>%
              ungroup() %>%
              mutate(phenotype = ifelse(phenotypes_count > 1,
                                        "Multiple phenotypes",
                                        normalised_specific_disease)) %>%
              dplyr::select(participant_id, phenotypes_count, phenotype) %>%
              distinct(),
            by = "participant_id") %>%
  distinct() %>%
  left_join(phenotypes_nejm, by = "phenotype") %>%
  left_join(panel_stats, by = "participant_id") %>%
  mutate(unique_genes_panels_applied_mb = unique_genes_panels_applied_bp / 1000000) %>%
  dplyr::select(-unique_genes_panels_applied_bp )
}

configure_penetrance <- function(meta, tiers) {
  meta %>%
    left_join(tiers %>%
                dplyr::select(participant_id, penetrance) %>%
                distinct() %>%
                group_by(participant_id, penetrance) %>%
                tally() %>%
                pivot_wider(names_from = penetrance, values_from = n) %>%
                mutate(penetrance = case_when(is.na(incomplete) ~ "complete",
                       TRUE ~ "incomplete")) %>%
                dplyr::select(participant_id, penetrance),
              by = "participant_id") %>%
    mutate(penetrance = ifelse(family_group_type == "Singleton", "complete", penetrance),
           penetrance = ifelse(is.na(penetrance), "complete", penetrance)) %>%
    dplyr::select(-biological_relationship_to_proband)
}

configure_affection <- function(meta) {
  meta %>%
    mutate(only_proband_affected = case_when(
            father_affected == "Yes" | mother_affected == "Yes" | full_brothers_affected > 0 | full_sisters_affected > 0 ~ "No",
            father_affected == "No" & mother_affected == "No" & full_brothers_affected == 0 & full_sisters_affected == 0 ~ "Yes",
            TRUE ~ "Unknown")) %>%
    dplyr::select(-c(father_affected, mother_affected, full_brothers_affected, full_sisters_affected))
}

meta_filter <- function(meta, assignments) {
  meta %>%
    filter(platekey %in% assignments$platekey,
           !is.na(all_variants_count))
}

relate_subsample_trios <- function(df, exclude_gia = "Remaining participants",
                                   nsamples = 10,
                                   kinship_coefficient = 0.0442, kinship_matrix = kinship_matrix) {
  df %>%
    filter(participant_type == "Proband" & family_group_type == "Trio with Mother and Father" & !participant_genetic_category %in% exclude_gia) %>%
    unrelated_filter(kinship_matrix = kinship_matrix,
                     kinship_coefficient = kinship_coefficient) %>% 
    group_by(participant_genetic_category) %>%
    slice_sample(n = nsamples) %>%
    ungroup() %>%
    rename(sample = platekey) %>%
    mutate(population = gsub(" ", "_", participant_genetic_category),
           group = population,
           sex = NA_character_) %>%
    dplyr::select(sample, population, group, sex)
}

summarise_plink_scount <- function(plink_scount) {
  plink_scount %>%
    rename(platekey = `#IID`) %>%
    group_by(platekey) %>%
    summarise(HOM_REF_CT = sum(HOM_REF_CT),
              HOM_ALT_CT = sum(HOM_ALT_CT),
              HET_CT = sum(HET_CT),
              NONSNP_CT = sum(NONSNP_CT),
              SINGLETON_GT = sum(SINGLETON_CT),
              MISSING_CT = sum(MISSING_CT))
}

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

unrelated_filter <- function(df, kinship_matrix, kinship_coefficient = 0.0442, show_rel = FALSE) {
  kin_df <- fread(kinship_matrix) %>%
    filter(
      IID1 %in% df$platekey &
      IID2 %in% df$platekey &
      KINSHIP > kinship_coefficient
    ) ## Heuristic
  if (nrow(kin_df) == 0)
    stop("No related pairs of individuals in the dataframe.\n")

  kin_df %>% group_by(IID1) %>% tally() %>% rename(IID = IID1) -> occ1
  kin_df %>% group_by(IID2) %>% tally() %>% rename(IID = IID2) -> occ2
  occ <- rbind(occ1, occ2) %>% group_by(IID) %>% summarise(n = sum(n))
  kin_df %<>% left_join(occ %>% rename(IID1 = IID, n1 = n), by = "IID1")
  kin_df %<>% left_join(occ %>% rename(IID2 = IID, n2 = n), by = "IID2") %>%
    mutate(n_mean = (n1 + n2) / 2) %>%
    arrange(desc(n_mean))

  cat("Attempting to remove minumum numbers of samples for unrelated cohort (heuristic)...\n")
  remove_vec <- c()
  for (row in 1:nrow(kin_df)) {
    ###
    print(row)
    id1 <- as.character(kin_df[row, 1])
    id2 <- as.character(kin_df[row, 2])
    if (id1 %in% remove_vec || id2 %in% remove_vec) {
    } else {
      id1_occ <- kin_df$n1[row]
      id2_occ <- kin_df$n2[row]
      if (id1_occ == id2_occ) {
      set.seed(123)
      remove_vec <- c(remove_vec, sample(c(id1, id2), 1))
      } else if (id1_occ > id2_occ) {
        remove_vec <- c(remove_vec, id1)
      } else if (id2_occ > id1_occ) {
        remove_vec <- c(remove_vec, id2)
      }
    }
  }
  cat(paste(
    "Number of relatives removed: ",
    length(remove_vec), "\n", sep = "")
  )
  df_norel <- df %>%
    filter(!platekey %in% remove_vec)

  if (show_rel == FALSE) {
    return(df_norel)
  } else if (show_rel == TRUE) {
    return(list(df_norel, remove_vec))
  }
}

configure_ethnicity_labels <- function(df) {
  df %>%
    mutate(participant_ethnic_category := case_when(
      is.na(participant_ethnic_category) ~ "Not Known or Not Stated",
      participant_ethnic_category == "Other Ethnic Groups: Chinese" |
      participant_ethnic_category == "Asian, Asian British or Asian Welsh: Chinese" ~
      "Asian or Asian British: Chinese",
      participant_ethnic_category == "Asian, Asian British or Asian Welsh: Bangladeshi" ~
      "Asian or Asian British: Bangladeshi",
      participant_ethnic_category == "Asian, Asian British or Asian Welsh: Indian" ~
      "Asian or Asian British: Indian",
      participant_ethnic_category == "Asian, Asian British or Asian Welsh: Pakistani" ~
      "Asian or Asian British: Pakistani",
      participant_ethnic_category == "Asian, Asian British or Asian Welsh: Other Asian" ~
      "Asian or Asian British: Any other Asian background",
      participant_ethnic_category == "Black, Black British, Black Welsh, Caribbean or African: African" ~
      "Black or Black British: African",
      participant_ethnic_category == "Black, Black British, Black Welsh, Caribbean or African: Caribbean" ~
      "Black or Black British: Caribbean",
      participant_ethnic_category == "Black, Black British, Black Welsh, Caribbean or African: Other Black" ~
      "Black or Black British: Any other Black background",
      participant_ethnic_category == "White: English, Welsh, Scottish, Northern Irish or British" ~
      "White: British",
      participant_ethnic_category == "White: Gypsy or Irish Traveller" ~ "White: Any other White background",
      participant_ethnic_category == "White: Roma" ~ "White: Any other White background",
      participant_ethnic_category == "Other ethnic group: Arab" ~ "Other Ethnic Groups: Any other ethnic group",
      participant_ethnic_category == "Not Stated" ~
      "Not Known or Not Stated",
      participant_ethnic_category == "White: Other White" ~
      "White: Any other White background",
      participant_ethnic_category == "Mixed or Multiple ethnic groups: White and Asian" ~
      "Mixed: White and Asian",
      participant_ethnic_category == "Mixed or Multiple ethnic groups: White and Black African" ~
      "Mixed: White and Black African",
      participant_ethnic_category == "Mixed or Multiple ethnic groups: White and Black Caribbean" ~
      "Mixed: White and Black Caribbean",
      participant_ethnic_category == "Mixed or Multiple ethnic groups: Other Mixed or Multiple ethnic groups" ~
      "Mixed: Any other mixed background",
      participant_ethnic_category == "Other ethnic group: Any other ethnic group" ~
      "Other Ethnic Groups: Any other ethnic group",
      TRUE ~ participant_ethnic_category)
    ) %>%
    separate(participant_ethnic_category,
      into = c("participant_ethnic_category_broad", "participant_ethnic_category_narrow"),
      sep = ": ", remove = FALSE,
      extra = "drop", fill = "right") %>%
    mutate(participant_ethnic_category_narrow = case_when(participant_ethnic_category_broad == "Not Known or Not Stated" ~ "Not Known or Not Stated",
                                                          TRUE ~ participant_ethnic_category_narrow),
           ethn_label = paste0(participant_ethnic_category_broad, "\n", participant_ethnic_category_narrow))
}

configure_genetic_labels <- function(df, assignments) {
  df %>%
    left_join(assignments %>%
                mutate(participant_genetic_category_broad = case_when(Assignment %in% c("Africa East", "Africa West", "Africa South") ~ "AFR",
                                                                      Assignment %in% c("Africa North", "Middle East") ~ "MID",
                                                                      Assignment %in% c("Philippines", "Asia East", "Japan") ~ "EAS",
                                                                      Assignment %in% c("Bangladesh", "Sri Lanka", "Pakistan") ~ "SAS",
                                                                      Assignment %in% c("Europe East", "Europe South West", "Europe North West", "Italy") ~ "EUR",
                                                                      Assignment %in% c("Ashkenazi") ~ "ASJ",
                                                                      TRUE ~ "Remaining participants"),
                       participant_genetic_category = ifelse(is.na(Assignment) | Assignment %in% c("Japan", "Finland", "South America"), "Remaining participants", Assignment),
                       participant_genetic_category = ifelse(participant_genetic_category == "Italy", "Mediterranean", participant_genetic_category)),
    by = "platekey")
}

#########################################################################
#########################################################################
###                                                                   ###
###                        PARTICIPANT DATA                           ###
###                                                                   ###
#########################################################################
#########################################################################

##-------------------------------------------------------------------------
## Meta Data
##-------------------------------------------------------------------------

##### Participant meta data (LabKey)
meta_data <- read_csv(paste0(out_dir, "labkey_data/meta_data.csv.gz"))

##-------------------------------------------------------------------------
## Missing PAVs from gnomAD data
##-------------------------------------------------------------------------

####### Prive UK Biobank genetically inferred ancestry (GIA) assignments
ukbb_assignments_only_data <- read_csv(paste0(out_dir, "ancestry_data/ukbb_assignments.csv")) %>%
  dplyr::select(platekey, Assignment)

#### Missing data
aggv2_miss <- fread(paste0(out_dir, "pavs_gnomad_miss_data/aggV2_all_gnomAD_missing_pavs.scount.gz"))
aggv2_miss_summary <- aggv2_miss %>%
  summarise_plink_scount() %>%
  mutate(dataset = "100kGP Probands")

#####################################################################
###### THESE FILES CANNOT BE MADE AVAILABLE TO ALL IN THE RE ########
#####################################################################

if (args[2] == "COVID_ACCESS") {
  cat("COVID-19 cohort access. Curating PAV missingness data from the aggV2 and COVID-19 cohort data files.\n")
  cat(paste0("Specifically, the files: ", out_dir, "pavs_gnomad_miss_data/aggCOVID_all_gnomAD_missing_pavs.scount.gz must be available!"))
  aggcovid_miss <- fread(paste0(out_dir, "pavs_gnomad_miss_data/aggCOVID_all_gnomAD_missing_pavs.scount.gz"))
  aggcovid_miss_summary <- aggcovid_miss %>%
    summarise_plink_scount() %>%
    mutate(dataset = "COVID-19 Cohort")

} else if (args[2] != "COVID_ACCESS") {
  cat(paste0("No COVID-19 cohort access\n", "Skipping the manipulation or combination of any COVID-19 cohort data."))
}

####################################################################
###### ---------------------------------------------------- ########
####################################################################

##-------------------------------------------------------------------------
## Combine stuff with meta data
##-------------------------------------------------------------------------

###### Tiering data (LabKey)
tiering_data <- read_tsv(paste0(out_dir, "labkey_data/tiering_data.csv.gz")) %>%
  add_varkey()

####### gnomADv3 GIA assignments
gnomad_assignments_only_data <- read_tsv(paste0(out_dir, "ancestry_data/gnomad_assignments.tsv"),
                                         col_names = FALSE) %>%
  dplyr::select(X1, X3) %>%
  rename(platekey = X1, participant_gnomad_category = X3) ## sort this out

##### Phenotypes
phenotypes_data <- read_csv(paste0(out_dir, "labkey_data/phenotypes_data.csv.gz"))
phenotypes_nejm_data <- read_csv("../public_data/phenotypes_nejm.csv")

##### Panel Meta Data
panels_meta_data <- read_csv(paste0(out_dir, "panel_app_data/panel_app_meta_data.csv.gz"))

#### ROH data
roh_data <- read_csv(paste0(out_dir, "roh_data/roh_data.csv.gz"))

### Combine all data...
meta_combined_data <- meta_data %>%
  configure_age_date() %>%
  configure_ethnicity_labels() %>%
  configure_genetic_labels(assignments = ukbb_assignments_only_data) %>%
  left_join(gnomad_assignments_only_data, by = "platekey") %>%
  configure_phenotype(phenotypes = phenotypes_data,
                      phenotypes_nejm = phenotypes_nejm_data,
                      panel_stats = panels_meta_data) %>%
  configure_penetrance(tiers = tiering_data) %>%
  configure_affection() %>%
  configure_diversity_statistics(missing_gnomad_pavs = aggv2_miss_summary,
                                 roh = roh_data) %>%
  meta_filter(assignments = ukbb_assignments_only_data)

probands_meta_combined_data <- meta_combined_data %>%
  filter(participant_type == "Proband")
nrow(probands_meta_combined_data)
### 29,425 families / probands

### Can only combine the COVID missing pavs data if they have COVID_ACCESS
if (args[2] == "COVID_ACCESS") { 
  ## Combine
  combined_missing_data <- rbind(aggv2_miss_summary %>% filter(platekey %in% probands_meta_combined_data$platekey),
                               aggcovid_miss_summary) %>%
  configure_genetic_labels(assignments = ukbb_assignments_only_data)

} else if (args[2] != "COVID_ACCESS") {
  ## Combine
  combined_missing_data <- rbind(aggv2_miss_summary %>% filter(platekey %in% probands_meta_combined_data$platekey)) %>%
    configure_genetic_labels(assignments = ukbb_assignments_only_data)
}

##-------------------------------------------------------------------------
## cPAVs data from tiering_data
##-------------------------------------------------------------------------

#### gnomADv2 AFs for tiered variants
tiering_frequency_gnomadv2_data <- fread(paste0(out_dir, "gnomad_afs_data/gnomad_v2_exomes_tiering.afs.gz"), fill = TRUE) %>%
  filter(filter == "PASS") %>%
  mutate(chromosome = gsub("chr", "", chromosome),
         AC_oth_gnomad_v2 = as.integer(AC_oth_gnomad_v2),
         AF_nfe_gnomad_v2 = AC_nfe_gnomad_v2 / AN_nfe_gnomad_v2,
         AF_fin_gnomad_v2 = AC_fin_gnomad_v2 / AN_fin_gnomad_v2,
         AF_asj_gnomad_v2 = AC_asj_gnomad_v2 / AN_asj_gnomad_v2,
         AF_afr_gnomad_v2 = AC_afr_gnomad_v2 / AN_afr_gnomad_v2,
         AF_sas_gnomad_v2 = AC_sas_gnomad_v2 / AN_sas_gnomad_v2,
         AF_eas_gnomad_v2 = AC_eas_gnomad_v2 / AN_eas_gnomad_v2,
         AF_oth_gnomad_v2 = AC_oth_gnomad_v2 / AN_oth_gnomad_v2,
         AF_amr_gnomad_v2 = AC_amr_gnomad_v2 / AN_amr_gnomad_v2
         ) %>%
  add_varkey(varid = FALSE) %>%
  dplyr::select(varkey, contains("AF_"))

#### gnomADv3 genomes for modernisation
tiering_frequency_gnomadv3_data <- fread(paste0(out_dir, "gnomad_afs_data/gnomad_v3_genomes_tiering.afs.gz"), fill = TRUE) %>%
  filter(filter == "PASS") %>%
  mutate(chromosome = gsub("chr", "", chromosome),
         AF_nfe_gnomad_v3 = AC_nfe_gnomad_v3 / AN_nfe_gnomad_v3,
         AF_fin_gnomad_v3 = AC_fin_gnomad_v3 / AN_fin_gnomad_v3,
         AF_asj_gnomad_v3 = AC_asj_gnomad_v3 / AN_asj_gnomad_v3,
         AF_afr_gnomad_v3 = AC_afr_gnomad_v3 / AN_afr_gnomad_v3,
         AF_sas_gnomad_v3 = AC_sas_gnomad_v3 / AN_sas_gnomad_v3,
         AF_eas_gnomad_v3 = AC_eas_gnomad_v3 / AN_eas_gnomad_v3,
         AF_oth_gnomad_v3 = AC_oth_gnomad_v3 / AN_oth_gnomad_v3,
         AF_amr_gnomad_v3 = AC_amr_gnomad_v3 / AN_amr_gnomad_v3
         ) %>%
  add_varkey(varid = FALSE) %>%
  dplyr::select(varkey, contains("AF_"))

#### Tiered variants frequency
tiering_frequency_data <- fread(paste0(out_dir, "labkey_data/tiering_frequency_data.csv.gz")) %>%
  add_varkey(varid = FALSE) %>%
  dplyr::select(varkey, uk10k_twinsuk, uk10k_alspac, gel_gl_6628_af) %>%
  mutate_at(vars(contains("_")), ~replace_na(as.numeric(.), 0)) %>%
  group_by(varkey) %>%
  summarise(uk10k_twinsuk = max(uk10k_twinsuk),
            uk10k_alspac = max(uk10k_alspac),
            gel_gl_6628_af = max(gel_gl_6628_af)) %>%
  tiers_new_freq_join(new_freq_data = tiering_frequency_gnomadv2_data) %>%
  tiers_new_freq_join(new_freq_data = tiering_frequency_gnomadv3_data)

#### Panel App statistics (new panels)
panel_app_stats_data <- fread("../public_data/panel_app_stats_data.csv.gz") %>%
  rename(genomic_feature_hgnc = gene_symbol) %>%
  ### Incorrect MOI for this gene
  mutate(mode_of_inheritance = ifelse(genomic_feature_hgnc == "RP1L1",
                                      "BOTH monoallelic and biallelic, autosomal or pseudoautosomal",
                                      mode_of_inheritance))

### PanelApp panels applied individual data
panels_applied_data <- fread(paste0(out_dir, "labkey_data/panels_applied_data.csv.gz"))

##### Karyotype data
karyotype_data <- meta_data %>%
  mutate(karyotype = case_when(is.na(karyotype) & participant_phenotypic_sex == "Female" ~ "XX",
                               is.na(karyotype) & participant_phenotypic_sex == "Male" ~ "XY",
                               TRUE ~ karyotype)) %>%
  dplyr::select(participant_id, karyotype)

### Specific SO terms from RD Pipeline
so_terms_vec <- c("SO:0001893", "SO:0001574", "SO:0001575", "SO:0001587", "SO:0001589", "SO:0001578",
                  "SO:0001582", "SO:0002012", "SO:0001889", "SO:0001821", "SO:0001822", "SO:0001583", 
                  "SO:0001630", "SO:0001626")

###### Only keep unique cPAVs identified per proband
###### Tier is the maximum (TIER1 > TIER2 > TIER3)
tiering_probands_retier_af_filter_data <- tiering_data %>%
  filter(participant_type == "Proband" &
         participant_id %in% meta_combined_data$participant_id) %>%
  tiers_harmonize_afs_filter(freq_data = tiering_frequency_data) %>%
  retier_variants(panels_applied = panels_applied_data,
                  panel_genes = panel_app_stats_data,
                  karyotypes = karyotype_data,
                  so_terms = so_terms_vec,
                  confidence = 3, keep_old_tier12 = TRUE) %>%
  retiers_replicates_class_filter()

# 36856 variants removed after not passing dominant AF filter!
# 4892 variants removed after not passing recessive AF filter!

##-------------------------------------------------------------------------
## Deleteriousness annotations (AlphaMissense, PrimateAI-3D, CADDv1.6)
##-------------------------------------------------------------------------

#### AlphaMissense predictions for tiered SNPs
#### 1,463,319 Tiered variants (AlphaMissense)
am_data <- fread(paste0(out_dir, "deleterious_predictions_data/alpha_missense_scores_tiered.tsv.gz")) %>%
  rename(AM_score = am_pathogenicity, AM_pred = am_class) %>%
  annot_filter(algo = "AM")

#### VEP annotations for tiered SNPs (CADD)
vep_data <- fread(paste0(out_dir, "deleterious_predictions_data/cadd_scores_tiered.tsv.gz"))

#### CADD score cutoff >= 30
cadd_data <- vep_data %>%
  cadd_filter()

#####################################################################
###### THESE FILES CANNOT BE MADE AVAILABLE TO ALL IN THE RE ########
#####################################################################

if (args[1] == "PAI3D_ACCESS") {
  cat("PrimateAI-3D access. Curating deleteriousness predictions from PrimateAI-3D scores.\n")
  cat(paste0("Specifically, the files: ", out_dir, "deleterious_predictions_data/primate_ai_3d_scores_tiered.tsv.gz must be available!"))
  #### PrimateAI-3D predictions for tiered SNPs
  pai_data <- fread(paste0(out_dir, "deleterious_predictions_data/primate_ai_3d_scores_tiered.tsv.gz")) %>%
    pai_filter()
} else if (args[1] != "PAI3D_ACCESS") {
  cat(paste0("No PrimateAI-3D access. Curating deleteriousness predictions from PrimateAI predictions (annotated using VEP).\n"),
      "NOTE: All plots and files will still say PrimateAI-3D and results will differ from the publication!!\n")
  #### PrimateAI scores
  pai_data <- vep_data %>%
    annot_filter(algo = "PrimateAI") %>%
    rename(PAI3D_pred = PrimateAI_pred)
}

####################################################################
###### ---------------------------------------------------- ########
####################################################################

### Combine data for export...
missense_data <- full_join(am_data, pai_data, by = "varkey")
predictions_data <- full_join(cadd_data, missense_data, by = "varkey")

##-------------------------------------------------------------------------
## Relate trees data
##-------------------------------------------------------------------------

relate_ss_inds <- relate_subsample_trios(meta_combined_data)

#########################################################################
#########################################################################
###                                                                   ###
###                          WRITE DATA                               ###
###                                                                   ###
#########################################################################
#########################################################################

data_dir <- paste0(out_dir, "combined_data/")
dir.create(data_dir)
cat(paste0("Writing data to: ", data_dir, "\n"))

#####################################################################
###### THESE FILES CANNOT BE MADE AVAILABLE TO ALL IN THE RE ########
#####################################################################

if (args[2] == "COVID_ACCESS") {
  
  cat(paste0("COVID-19 cohort access. Writing unrelated individuals to FAM file: ", data_dir, "aggCOVID_unrelated_pops.fam"))
  cat(paste0("If so, please ensure that the file: ", kinship_matrix, " also contains individuals from the aggCOVID dataset."))
  ##-------------------------------------------------------------------------
  ## COVIDv5 samples
  ##-------------------------------------------------------------------------

  ### <3rd degree relatives removed
  unrelated_individuals_covid_fam <- combined_missing_data %>%
    dplyr::select(platekey, participant_genetic_category, dataset) %>%
    unrelated_filter(kinship_matrix = kinship_matrix,
                     kinship_coefficient = 0.0884) %>%
    filter(dataset == "COVID-19 Cohort") %>%
    group_by(participant_genetic_category) %>%
    mutate(n_in_pop = n()) %>%
    ungroup() %>%
    filter(n_in_pop >= 100 & participant_genetic_category != "Remaining participants") %>%
    dplyr::select(-n_in_pop, dataset) %>%
    mutate(father_id = 0, mother_id = 0, sex_code = 0, pheno_val = 1,
           population = gsub(" ", "_", participant_genetic_category),
           platekey2 = platekey) %>%
    dplyr::select(platekey, platekey2, father_id, mother_id, sex_code, pheno_val, population)

    ### aggCOVIDv5 unrelated individuals poplabels (FAM)
    write.table(unrelated_individuals_covid_fam, file.path(paste0(data_dir, "aggCOVID_unrelated_pops.fam")),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

### Number of relatives removed: 1019
} 

####################################################################
###### ---------------------------------------------------- ########
####################################################################

### Combined tiering frequency data with gnomADv3, gnomADv2 etc.
fwrite(tiering_frequency_data,
       file.path(paste0(data_dir, "combined_tiering_frequency_data.csv.gz")))

### Retiered data filtered by updated AFs and using updated PanelApp panels
fwrite(tiering_probands_retier_af_filter_data ,
       file.path(paste0(data_dir, "filtered_retiering_probands_data.csv.gz")))

### Combined dataframe of deleteriousness predictions (CADD, ALPHAMISSENSE, PRIMATEAI-3D)
fwrite(predictions_data,
       file.path(paste0(data_dir, "combined_predictions_data.csv.gz")))

### Meta data with all statistics used for GLMs
fwrite(meta_combined_data,
       file.path(paste0(data_dir, "combined_meta_data.csv.gz")))

### Combined missingness data for aggCOVID + aggV2
fwrite(combined_missing_data,
       file.path(paste0(data_dir, "combined_missing_pavs_gnomad_nocovid.csv.gz")))

### Pop labels for running relate.
write.table(relate_ss_inds, file.path(paste0(data_dir, "relate_ss_poplabels.tsv")),
            quote = FALSE, sep = " ", row.names = FALSE)


###########################################################################
###########################################################################
###                                                                     ###
###                         END OF SCRIPT                               ###
###                                                                     ###
###########################################################################
###########################################################################

