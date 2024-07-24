#!/usr/bin/env Rscript

##-------------------------------------------------------------------------
## Configs + packages

source("../config/config.R")

library(tidyverse)
library(data.table)
library(parallel)
library(Rlabkey)

# LabKey settings
labkey.setDefaults(baseUrl = baseUrl)
labkey.setDebugMode(FALSE)
labkey.setWafEncoding(FALSE)

###########################################################################
###########################################################################
###                                                                     ###
###                           FUNCTIONS                                 ###
###                                                                     ###
###########################################################################
###########################################################################

## General SQL query function (required for parallelisation)
sqlquery <- function(sql) {
  labkey.executeSql(
    sql = sql,
    maxRows = 1000000000,
    folderPath = labkey,
    schemaName = "lists",
    colNameOpt = "rname")
}

pull_tiering <- function(sql = NULL, labkey) {
  #Tiering data table is problematically large.
  #Loop through the table to do this
  #Go chrom by chrom
  chroms <- c(as.character(seq(1:22)), 'X', 'Y', 'MT')
  nCores <- detectCores() - 1
  outList <- mclapply(chroms, function(x) {
    sql <-  paste0("SELECT
      t.participant_id,
      t.rare_diseases_family_id,
      t.participant_type,
      t.tier,
      t.penetrance,
      t.segregation_pattern,
      t.mode_of_inheritance,
      t.phenotype,
      t.chromosome,
      t.position,
      t.reference,
      t.alternate,
      t.genotype,
      t.genomic_feature_hgnc,
      t.ensembl_id,
      t.so_term,
      t.consequence_type,
      t.father_affected,
      t.mother_affected,
      t.full_brothers_affected,
      t.full_sisters_affected,
      FROM tiering_data AS t
      WHERE t.chromosome = ('", x, "') AND t.assembly = 'GRCh38'")
    sqlquery(sql)
  }, mc.cores = nCores) ## add consequence type as optional parameter
  names(outList) <- chroms
  out <- Reduce("rbind", outList) %>%
    as_tibble() %>%
    mutate(participant_id = as.numeric(participant_id))
  return(out)
}

pull_tiering_frequency <- function(sql = NULL, labkey) {
  if (is.null(sql)) {
    sql <- paste0("SELECT
    tf.chromosome,
    tf.position,
    tf.reference,
    tf.alternate,
    tf.uk10k_twinsuk,
    tf.uk10k_alspac,
    tf.gnomad_exomes_nfe,
    tf.gnomad_exomes_afr,
    tf.gnomad_exomes_eas,
    tf.gnomad_exomes_amr,
    tf.gnomad_exomes_asj,
    tf.gnomad_exomes_oth,
    tf.gnomad_exomes_fin,
    tf.gnomad_genomes_afr,
    tf.gnomad_genomes_nfe,
    tf.gnomad_genomes_fin,
    tf.gnomad_genomes_eas,
    tf.gnomad_genomes_oth,
    tf.gnomad_genomes_amr,
    tf.gel_gl_6628_af,
    tf.gel_gl_6628_mt_af,
    FROM tiered_variants_frequency as tf
    WHERE tf.assembly = 'GRCh38'")
  }
  out <- sqlquery(sql) %>%
    as_tibble()
  return(out)
}

mixedrank <- function(x) order(gtools::mixedorder(x))

get_sorted_positions <- function(tiers) {
  variants_tab <- tiers %>%
    dplyr::select(chromosome, position,
                  reference, alternate) %>%
    distinct() %>%
    mutate(chromosome = paste0("chr", chromosome))
  variants_tab %>%
    filter(!chromosome %in% c("chrM", "chrX")) %>%
    dplyr::arrange(position) %>%
    dplyr::arrange(mixedrank(chromosome)) -> autosomes_tab
  variants_tab %>%
    filter(chromosome == "chrX") %>%
    dplyr::arrange(position) -> chrX_tab
  variants_tab %>%
    filter(chromosome == "chrM") %>%
    dplyr::arrange(position) -> chrM_tab
  rbind(autosomes_tab, chrX_tab, chrM_tab) %>%
    unite("allelle", reference:alternate, sep = ",", remove = TRUE)
}

##-------------------------------------------------------------------------
## Participant Metadata Pull -
##-------------------------------------------------------------------------

pull_participant <- function(sql = NULL, labkey) {
  if (is.null(sql)) {
    sql <- paste0("SELECT
    p.participant_id,
    p.participant_type,
    p.rare_diseases_family_id,
    p.year_of_birth,
    p.mother_affected,
    p.father_affected,
    p.full_brothers_affected,
    p.full_sisters_affected,
    p.fathers_ethnic_category_other,
    p.fathers_other_relevant_ancestry,
    p.mothers_ethnic_category_other,
    p.mothers_other_relevant_ancestry,
    p.participant_phenotypic_sex,
    p.handling_gmc_trust,
    p.participant_ethnic_category,
    p.biological_relationship_to_proband,
    FROM participant as p
    WHERE p.programme = 'Rare Diseases'")
  }
  out <- sqlquery(sql) %>%
    as_tibble() %>%
    mutate(participant_id = as.numeric(participant_id),
           rare_diseases_family_id = as.character(rare_diseases_family_id))
  return(out)
}

pull_bioinfo <- function(sql = NULL, labkey) {
  if (is.null(sql)) {
    sql <- paste0("SELECT
    agg.participant_id,
    agg.karyotype,
    agg.illumina_autosome_mean_coverage,
    agg.samtools_error_rate,
    agg.samtools_reads_mapped,
    agg.samtools_reads_unmapped,
    agg.illumina_snvs_all,
    agg.illumina_indels_all,
    agg.illumina_snv_het_hom_ratio,
    agg.illumina_indel_het_hom_ratio,
    FROM aggregate_gvcf_sample_stats as agg")
  }
  out <- sqlquery(sql) %>%
    as_tibble() %>%
    mutate(participant_id = as.numeric(participant_id))
  return(out)
}

pull_fp_platekey <- function(sql = NULL, labkey) {
  if (is.null(sql)) {
    sql <- paste0("SELECT
      fp.participant_id,
      fp.platekey,
      fp.file_path,
      fp.delivery_date,
      FROM genome_file_paths_and_types as fp
      WHERE fp.file_sub_type = 'Standard VCF'
      AND fp.type IN ('rare disease germline', 'experimental germline', 'rare disease unknown', 'cancer germline')")
  }
  out <- sqlquery(sql) %>%
    as_tibble() %>%
    mutate(participant_id = as.numeric(participant_id)) %>%
    group_by(participant_id) %>%
    filter(delivery_date == max(delivery_date)) %>%
    ungroup()
  return(out)
}

pull_rd_family <- function(sql = NULL, labkey) {
  if (is.null(sql)) {
    sql <- paste0("SELECT
      fam.rare_diseases_family_id,
      fam.family_group_type,
      fam.multiple_monogenic_likely_code,
      FROM rare_diseases_family as fam")
  }
  out <- sqlquery(sql) %>%
    as_tibble() %>%
    mutate(rare_diseases_family_id = as.character(rare_diseases_family_id))
  return(out)
}

#### Combine SQL datasets
pull_meta <- function(labkey) {
  ### Pull meta data
  cat("Pulling participant table data...\n")
  participant <- pull_participant(labkey = labkey)
  cat("Pulling bioinformatics data...\n")
  bioinfo <- pull_bioinfo(labkey = labkey)
  fp_platekey <- pull_fp_platekey(labkey = labkey)
  cat("Pulling rare diseases family meta data...\n")
  rd_family <- pull_rd_family(labkey = labkey)

  cat("Joining meta data...\n")
  meta <- participant %>%
    ### Pull datasets
    left_join(bioinfo, by = "participant_id") %>%
    left_join(fp_platekey, by = "participant_id") %>%
    left_join(rd_family, by = "rare_diseases_family_id") %>%
    drop_na(file_path)

  return(meta)
}

##-------------------------------------------------------------------------
## Phenotype data pull
##-------------------------------------------------------------------------

pull_phenotypes <- function(sql = NULL, labkey) {
  if (is.null(sql)) {
    sql <- paste0("SELECT
    d.participant_id,
    d.normalised_specific_disease,
    d.normalised_disease_sub_group,
    d.normalised_disease_group
    FROM rare_diseases_participant_disease as d")
  }
  out <- sqlquery(sql) %>%
    as_tibble() %>%
    mutate(participant_id = as.numeric(participant_id))
  return(out)
}

##-------------------------------------------------------------------------
## Panels data pull
##-------------------------------------------------------------------------

pull_panels <- function(sql = NULL, labkey) {
  if (is.null(sql)) {
    sql <- paste0("SELECT
    pan.participant_id,
    pan.panel_name,
    FROM panels_applied as pan")
  }
  out <- sqlquery(sql) %>%
    as_tibble() %>%
    mutate(participant_id = as.numeric(participant_id)) %>%
    ## A few names with non-standard characters require adjustment
    mutate(panel_name = case_when(panel_name == "Anophthalmia/microphthamia" ~ "Anophthalmia",
                                  panel_name == "Significant early-onset obesity +/- other endocrine features and short stature" ~ "Severe early-onset obesity",
                                  panel_name == "Renal tract calcification (or Nephrolithiasis/nephrocalcinosis)" ~ "Nephrolithiasis",
                                  panel_name == "Amyotrophic lateral sclerosis/motor neuron disease" ~ "Amyotrophic lateral sclerosis",
                                  panel_name == "Congenital hearing impairment (profound/severe)" ~ "Congenital hearing impairment",
                                  TRUE ~ panel_name),
           panel_sub = gsub("/", "", gsub(" ", "_", panel_name)))
  return(out)
}

##-------------------------------------------------------------------------
## Exit Questionnaire Pull -
##-------------------------------------------------------------------------

pull_exit <- function(sql = NULL, labkey) {
  if (is.null(sql)) {
    sql <- paste0("SELECT
    ex.participant_id,
    ex.case_solved_family,
    ex.confirmation_decision,
    ex.event_date,
    ex.assembly,
    ex.chromosome,
    ex.acmg_classification,
    ex.position,
    ex.reference,
    ex.alternate,
    FROM gmc_exit_questionnaire as ex")
  }
  out <- sqlquery(sql) %>%
    as_tibble() %>%
    mutate(participant_id = as.numeric(participant_id),
      reference = na_if(reference, "NA"),
      alternate = na_if(alternate, "NA"),
      chromosome = gsub("chr", "", chromosome)) %>%
    distinct()
  return(out)
}

###########################################################################
###########################################################################
###                                                                     ###
###                        RUN DATA PULL                                ###
###                                                                     ###
###########################################################################
###########################################################################

## Project started with v17 data release..
labkey <- labkey_v17
cat("Pulling tiering data...\n")
cat("May take a while (10 minutes).\n")
tiering_data <- pull_tiering(labkey = labkey)

cat("Pulling tiering frequency...\n")
tiering_frequency_data <- pull_tiering_frequency(labkey = labkey)

cat("Making sorted positions file from tiering data...\n")
tiering_variants_data <- get_sorted_positions(tiers = tiering_frequency_data)

meta_data <- pull_meta(labkey = labkey)

cat("Pulling disease data...\n")
phenotypes_data <- pull_phenotypes(labkey = labkey)

cat("Pulling panels data...\n")
panels_applied_data <- pull_panels(labkey = labkey)

## Use most up to date exit Q.
labkey <- labkey_v18
cat("Pulling exit data...\n")
exit_data <- pull_exit(labkey = labkey)

#########################################################################
#########################################################################
###                                                                   ###
###                          WRITE DATA                               ###
###                                                                   ###
#########################################################################
#########################################################################

data_dir <- paste0(out_dir, "labkey_data/")

## Create directory.
dir.create(file.path(data_dir))
cat(paste0("Writing data to: ", data_dir, "\n"))

### Original tiering data (unfiltered)
fwrite(tiering_data,
       file.path(paste0(data_dir, "tiering_data.csv.gz")))

### Some tiering frequency data information from the original pipeline
fwrite(tiering_frequency_data,
       file.path(paste0(data_dir, "tiering_frequency_data.csv.gz")))

### Tiering sorted variants (for VEP pull + extracting tiering variants from AF data)
write.table(tiering_variants_data,
            file.path(paste0(data_dir, "tiering_sorted_variants.tsv")),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

### Exit data for ACMG criteria + diagnoses
fwrite(exit_data,
       file.path(paste0(data_dir, "exit_data.csv.gz")))

### General participant meta data
fwrite(meta_data,
       file.path(paste0(data_dir, "meta_data.csv.gz")))

### Data on participant rare disease phenotypes
fwrite(phenotypes_data,
       file.path(paste0(data_dir, "phenotypes_data.csv.gz")))

### Data on the original panels applied to probands.
fwrite(panels_applied_data,
       file.path(paste0(data_dir, "panels_applied_data.csv.gz")))

###########################################################################
###########################################################################
###                                                                     ###
###                         END OF SCRIPT                               ###
###                                                                     ###
###########################################################################
###########################################################################