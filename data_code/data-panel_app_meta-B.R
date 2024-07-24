#!/usr/bin/env Rscript

##-------------------------------------------------------------------------
## Configs + packages


source("../config/config.R")

library(tidyverse)
library(data.table)

###########################################################################
###########################################################################
###                                                                     ###
###                           FUNCTIONS                                 ###
###                                                                     ###
###########################################################################
###########################################################################

## Statistics for panel data...
get_panel_meta_data <- function(panels_applied, panel_app_stats,
                                confidence_score = 3) {
  ### Participant IDs
  participant_ids <- unique(panels_applied_data$participant_id)
  inds_panels_list <- vector(mode = "list",
                             length = length(participant_ids))
  panel_app_stats_filt <- panel_app_stats %>%
    filter(confidence >= confidence_score)

  counter <- 1
  for (participant_index in 1:length(participant_ids)) {
    participant <- participant_ids[participant_index]
    cat(paste0(counter, ": ", participant, "\n"))

    panels <- panels_applied %>%
      filter(participant_id == participant) %>%
      pull(panel_name) %>%
      unique()

    panels_applied_genes_bp <- panel_app_stats_filt %>%
      filter(panel_name %in% panels) %>%
      distinct()

    inds_panels_list[[counter]] <- data.frame(
      participant_id = participant,
      panels_applied_count = length(panels),
      unique_genes_panels_applied_count = nrow(panels_applied_genes_bp),
      unique_genes_panels_applied_bp = sum(panels_applied_genes_bp$bp)
    )
    counter <- counter + 1
  }
  inds_panel_data <- do.call("rbind", inds_panels_list)
  return(inds_panel_data)
}

###########################################################################
###########################################################################
###                                                                     ###
###                          DATA INGEST                                ###
###                                                                     ###
###########################################################################
###########################################################################

### PanelApp panels applied individual data
panels_applied_data <- fread(paste0(out_dir, "labkey_data/panels_applied_data.csv.gz"))

### PanelApp stats data for new panels
panel_app_stats_data <- fread("../public_data/panel_app_stats_data.csv.gz")

###########################################################################
###########################################################################
###                                                                     ###
###                          DATA CREATE                                ###
###                                                                     ###
###########################################################################
###########################################################################

## Panel Meta Data
cat("Calculating individual statistics based on PanelApp panels (may take a while)...\n")
panel_app_meta_data <- panels_applied_data %>%
  get_panel_meta_data(panel_app_stats = panel_app_stats_data,
                      confidence_score = 3)

###########################################################################
###########################################################################
###                                                                     ###
###                          WRITE DATA                                 ###
###                                                                     ###
###########################################################################
###########################################################################

data_dir <- paste0(out_dir, "panel_app_data/")

## Create directory.
dir.create(file.path(data_dir))
cat(paste0("Writing data to: ", data_dir, "\n"))

fwrite(panel_app_meta_data,
       file.path(paste0(data_dir, "panel_app_meta_data.csv.gz")))

###########################################################################
###########################################################################
###                                                                     ###
###                         END OF SCRIPT                               ###
###                                                                     ###
###########################################################################
###########################################################################