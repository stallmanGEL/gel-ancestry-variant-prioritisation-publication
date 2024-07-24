#!/usr/bin/env Rscript

##-------------------------------------------------------------------------
## Configs + packages

source("../config/config.R")

library(tidyverse)
library(bigsnpr)
library(tools)
library(data.table)
library(magrittr)
library(umap)

###########################################################################
###########################################################################
###                                                                     ###
###                           FUNCTIONS                                 ###
###                                                                     ###
###########################################################################
###########################################################################


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

###########################################################################
###########################################################################
###                                                                     ###
###                           DATA INGEST                               ###
###                                                                     ###
###########################################################################
###########################################################################

### Read in bed file and get only meta RD programme inds
(obj.bed <- bed(paste0(hq_snps_plink, ".bed")))

meta_combined_data <- read_csv(paste0(out_dir, "combined_data/meta_combined_data.csv.gz"))

### Get min list of unrelated inds (heuristic)
rel_list <- unrelated_filter(meta_combined_data,
    kinship_matrix = kinship_matrix,
    show_rel = TRUE)
ind.rel <- match(rel_list[[2]], obj.bed$fam$sample.ID)
ind.norel <- match(rel_list[[1]]$platekey, obj.bed$fam$sample.ID)

###########################################################################
###########################################################################
###                                                                     ###
###                           ANALYSIS RUN                              ###
###                                                                     ###
###########################################################################
###########################################################################

## Run PCA main on unrelated inds
obj.svd <- bed_autoSVD(obj.bed, ind.row = ind.norel, k = 20)

pcs_matrix <- matrix(NA, nrow(obj.bed), ncol(obj.svd$u))
pcs_matrix[ind.norel, ] <- predict(obj.svd)

### Project outliers (none here) and unrelated individuals
proj <- bed_projectSelfPCA(obj.svd, obj.bed,
                           ind.row = ind.rel,
                           ncores = 1)
pcs_matrix[ind.rel, ] <- proj$OADP_proj

## Get only RD programme inds
pcs_df <- as.data.frame(pcs)
pcs_df$platekey <- obj.bed$fam$sample.ID
pcs_df_rd <- na.omit(pcs_df)

## Run UMAP on top 16 PCs
umap_rd <- umap(pcs_df_rd[, 1:16] %>% as.matrix(),
                n_neighbours = 15,
                min_dist = 0.4)

pca_umap_df_rd <- cbind(pcs_df_rd[, c(1:16, 21)], umap_rd[[1]])
colnames(pca_umap_df_rd) <- c(paste0("PC", 1:16),
                       "platekey",
                       paste0("UMAP", 1:2))

#########################################################################
#########################################################################
###                                                                   ###
###                          WRITE DATA                               ###
###                                                                   ###
#########################################################################
#########################################################################

## Analysis data
data_dir <- paste0(analysis_dir, "pca_analysis/")

## Create directory.
dir.create(file.path(data_dir))
cat(paste0("Writing analyses to: ", data_dir, "\n"))

## Write PC/UMAP and PC loadings
write.csv(pca_umap_df_rd,
          paste0(data_dir, "pca_umap.csv"),
          row.names = FALSE)

###########################################################################
###########################################################################
###                                                                     ###
###                         END OF SCRIPT                               ###
###                                                                     ###
###########################################################################
###########################################################################