#!/usr/bin/env Rscript

##-------------------------------------------------------------------------
## Configs + packages

source("../config/config.R")

library(tidyverse)
library(data.table)
library(bignspr)
library(tools)

###########################################################################
###########################################################################
###                                                                     ###
###                           FUNCTIONS                                 ###
###                                                                     ###
###########################################################################
###########################################################################

### Projection function
project_genomes <- function(bedfile, all_freq, loadings, correction, tmpDir) {
  ## Match alleles
  bim <- sub_bed(bedfile, ".bim") %>%
    bigreadr::fread2(select = c(1, 4:6),
                   col.names = c("chr", "pos", "a1", "a0")) %>%
    mutate(beta = 1, chr = as.character(chr))

  matched <- bim %>%
    snp_match(all_freq[1:4])

  ## Get genotypes from bedfile (WARNING: CREATES TEMPORARY FILE)
  random_number <- paste0(as.character(sample.int(100, 10)), collapse = "")
  backingfile <- paste0(tmpDir, "temp", random_number)
  rds <- snp_readBed2(
    bedfile,
    backingfile = backingfile,
    ind.col = matched$`_NUM_ID_.ss`
  )
  obj.bigsnp <- snp_attach(rds)

  ## Rapidly impute missing genotypes
  ## exclude those with lots of missing genotypes
  G <- obj.bigsnp$genotypes
  five_perc_miss <- 0.05 * nrow(obj.bigsnp$fam)

  nb_na <- big_counts(G)[4, ]
  ind <- which(nb_na < five_perc_miss)
  cat(paste0(length(ind), " SNPs remaining after filtering on >5% missingness\n"))

  G2 <- snp_fastImputeSimple(G)

  ### Project samples onto PC space
  PROJ <- as.matrix(loadings[matched$`_NUM_ID_`[ind], -(1:4)])
  all_proj <- big_prodMat(G2, sweep(PROJ, 2, correction / 2, '*'),
    ind.col = ind,
    # scaling to get G if beta = 1 and (2 - G) if beta = -1
    center = 1 - matched$beta[ind],
    scale = matched$beta[ind]
  ) %>%
  data.frame() %>%
  drop_na() %>%
  as.matrix()
  rownames(all_proj) <- obj.bigsnp$fam[, c(2)]

  ### Create matrix cross product
  X <- crossprod(PROJ,
    as.matrix(all_freq[matched$`_NUM_ID_`[ind], -(1:4)])
  )
  system(glue::glue("rm {backingfile}.bk"))
  system(glue::glue("rm {backingfile}.rds"))
  ### Return labelled matrix without missing samples
  return(list(X, all_proj))
}

### Assign participants and get approximate FST function
assign_euclidean <- function(X, all_proj, group, threshold) {
  # Get pop centres from matrix crossprod
  all_centers <- t(X)

  ### squared distance to centres to assign
  max_sq_dist <- max(dist(all_centers)^2) / 0.16

  ## Find square distance to centres from projection
  all_sq_dist <- apply(all_centers, 1, function(one_center) {
    rowSums(sweep(all_proj, 2, one_center, '-')^2)
  })
  all_sq_dist <- all_sq_dist / max_sq_dist

  if (is.null(dim(all_sq_dist))) {
    ind <- which.min(all_sq_dist)
    cluster <- if (isTRUE(all_sq_dist[ind] < threshold)) names(ind) else NA
  } else {
    cluster <- apply(all_sq_dist, 1, function(x) {
        ind <- which.min(x)
        if (isTRUE(x[ind] < threshold)) ind else NA
      })
  }
  if (!all(is.na(cluster))) {
    cluster <- group[cluster]
  }
  col <- colnames(all_sq_dist)
  return(cbind(data.frame(
    platekey = rownames(all_proj),
    Assignment = cluster
    ), data.frame(all_sq_dist) %>%
      setNames(paste("Approximate FST to", col, sep = " "))
  ))
}

#### Proportional genetic similarity function
assign_admixture <- function(X, all_proj) {
  # Generate inputs for QP
  cp_X_pd <- Matrix::nearPD(crossprod(X), base.matrix = TRUE)
  Amat <- cbind(1, diag(ncol(X)))
  bvec <- c(1, rep(0, ncol(X)))
  grp_fct <- factor(group, unique(group))

  # solve a QP for each projected individual
  all_res <- apply(all_proj, 1, function(y) {
    quadprog::solve.QP(
      Dmat = cp_X_pd$mat,
      dvec = crossprod(y, X),
      Amat = Amat,
      bvec = bvec,
      meq  = 1
    )$sol %>%
    tapply(grp_fct, sum) %>%
    round(7)
  })
  all_res <- t(all_res)

  all_prop <- all_res %>%
    as.data.frame() %>%
    setNames(paste("Proportional similarity to", colnames(.), sep = " ")) %>%
    mutate(platekey = rownames(.), .before = colnames(.)[1])

  return(all_prop)
}

prop_adjust_assignment <- function(df, prop, assign_col, pops_vec) {
  assign_col_sym <- rlang::sym(assign_col)
  columns_props <- grepl("Proportional", colnames(df))
  df %>% mutate(
    max = apply(df[, ..columns_props] %>% dplyr::select(contains(pops_vec)), 1, max),
    max_pop = pops_vec[apply(df[, ..columns_props] %>% dplyr::select(contains(pops_vec)), 1, which.max)],
    {{assign_col_sym}} := case_when(
        max >= prop & max_pop == {{assign_col_sym}} ~ {{assign_col_sym}},
        max >= prop & is.na({{assign_col_sym}}) ~ max_pop,
        TRUE ~ {{assign_col_sym}}
    )) %>%
    dplyr::select(-c(max, max_pop))
}

###########################################################################
###########################################################################
###                                                                     ###
###                           DATA INGEST                               ###
###                                                                     ###
###########################################################################
###########################################################################

## PLINK Bed file (HQ snps)
bedfile <- paste0(hq_snps_plink, ".bed")
cat(paste0("If you have access to COVID-19 cohort, ensure that the file: ", hq_snps_plink, " also contains these individuals!"))

### Projection
### Read in PC loadings and reference alleles frequencies for GRCh38
all_freq_ukbb <- bigreadr::fread2(prive_ukbb_build38_afs) %>%
  dplyr::select(-rsid) %>%
  mutate(chr = as.character(chr))

loadings_ukbb <- bigreadr::fread2(prive_ukbb_build38_loadings) %>%
  dplyr::select(-rsid) %>%
  mutate(chr = as.character(chr))

## We will apply the correction for PC shrinkage calculated by Prive et al. 2022
n_pcs <- 16
correction_full <- c(1, 1, 1, 1.008, 1.021, 1.034, 1.052, 1.074, 1.099,
                     1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443)
correction <- correction_full[1:n_pcs]

## As in Prive et al. 2021
group <- colnames(all_freq_ukbb)[-(1:4)]
group[group %in% c("Scandinavia", "United Kingdom", "Ireland")]   <- "Europe North West"
group[group %in% c("Europe South East", "Europe North East")] <- "Europe East"

###########################################################################
###########################################################################
###                                                                     ###
###                           ANALYSIS RUN                              ###
###                                                                     ###
###########################################################################
###########################################################################

### Run projection
system.time(
  projection_files <- project_genomes(
    bedfile = bedfile,
    all_freq = all_freq_ukbb,
    loadings = loadings_ukbb,
    correction = correction,
    tmpDir = tmpDir
))

### Gather reference (X) projection and aggV2 (all_proj)
X <- projection_files[[1]]
all_proj <- projection_files[[2]]

### Get FST to references
system.time(
  fst_aggV2 <- assign_euclidean(
    X = X,
    all_proj = all_proj,
    group = group,
    threshold = 0.002
  )
)

### Run to get proportional genetic similarity
system.time(
  mixtures_aggV2 <- assign_admixture(
    X = X,
    all_proj = all_proj
  )
)

### Merge FST and proportional genetic similarity info
final_table <- fst_aggV2 %>%
  left_join(mixtures_aggV2, by = "platekey")

#### Write the adjusted assignments + projected PCA
prive_assignments <- final_table %>%
  prop_adjust_assignment(prop = 0.8, assign_col = "Assignment")

prive_assignments_proj <- cbind(prive_assignments, all_proj)

#########################################################################
#########################################################################
###                                                                   ###
###                          WRITE DATA                               ###
###                                                                   ###
#########################################################################
#########################################################################

data_dir <- paste0(out_dir, "ancestry_data/")

## Create directory.
dir.create(file.path(data_dir))
cat(paste0("Writing data to: ", data_dir, "\n"))

##### Write full data
prive_assignments_proj %>%
  write.csv(
    file.path(paste0(data_dir, "ukbb_assignments.csv")),
    row.names = FALSE
  )

###########################################################################
###########################################################################
###                                                                     ###
###                         END OF SCRIPT                               ###
###                                                                     ###
###########################################################################
###########################################################################