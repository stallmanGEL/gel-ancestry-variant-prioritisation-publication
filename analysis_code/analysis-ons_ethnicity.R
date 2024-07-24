#!/usr/bin/env Rscript

##-------------------------------------------------------------------------
## Configs + packages

source("../config/config.R")

library(tidyverse)
library(data.table)
library(magrittr)
#library(rstatix)

###########################################################################
###########################################################################
###                                                                     ###
###                           FUNCTIONS                                 ###
###                                                                     ###
###########################################################################
###########################################################################

run_rd_vs_ons_prop_test <- function(ons_rd_age_sex_ethn_matrix) {
  prop_test_list <- list()
  COUNTER <- 1
  age_bin_vec <- ons_rd_age_sex_ethn_matrix  %>%
    pull(age_bin) %>%
    unique()
  ethnicity_vec <- ons_rd_age_sex_ethn_matrix %>%
    pull(participant_ethnic_category) %>%
    unique()

  ons_rd_age_ethn_refmt <- ons_rd_age_sex_ethn_matrix %>%
      group_by(participant_ethnic_category, age_bin, group) %>%
      summarise(num_inds = sum(num_inds))

  for (bin in age_bin_vec) {
    ons_rd_agebin_ethn_refmt <- ons_rd_age_ethn_refmt %>%
      filter(age_bin == bin) %>%
      group_by(group) %>%
      mutate(num_inds_remainder = sum(num_inds),
             num_inds_remainder = num_inds_remainder - num_inds)

    for (ethnicity in ethnicity_vec) {
      print(COUNTER)
      xtab <- ons_rd_agebin_ethn_refmt %>%
        filter(participant_ethnic_category == ethnicity) %>%
        dplyr::select(-c(participant_ethnic_category, age_bin)) %>%
        column_to_rownames(var = "group")

      prop_test_list[[COUNTER]] <- prop_test(xtab, detailed = TRUE) %>%
        mutate(age_bin = bin,
               participant_ethnic_category = ethnicity) %>%
        rename(`100KGP Probands` = estimate1,
               `ONS 2021` = estimate2)
      COUNTER <- COUNTER + 1
    }
  }
  return(do.call("rbind", prop_test_list) %>%
           rename(p.value = p) %>%
           mutate(p.adjust = p.value * COUNTER,
                  estimate = (`100KGP Probands` - `ONS 2021`) / `ONS 2021`,
                  lower = conf.low / `ONS 2021`,
                  upper = conf.high / `ONS 2021`))
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

###########################################################################
###########################################################################
###                                                                     ###
###                       DATA INGES + WRANGLING                        ###
###                                                                     ###
###########################################################################
###########################################################################

meta_combined_data <- read_csv(paste0(out_dir, "combined_data/combined_meta_data.csv.gz"))

## Number of 100k Rare Disease probands collected in England per sex, ethnicity, and age bin
## age at date of ONS census (21/03/2021) given participant's age at delivery date
rd_age_sex_ethn_matrix <- meta_combined_data %>%
  filter(participant_type == "Proband" &
         !handling_gmc_trust %in% c("Wales", "Scotland", "Northern Ireland") &
         participant_ethnic_category != "Not Known or Not Stated") %>%
  mutate(age = age + round(as.numeric(as.Date("2021-03-21") -
                                      as.Date(delivery_date)) / 365.25)) %>%
  group_by(age, participant_ethnic_category, participant_phenotypic_sex) %>%
  summarise(num_inds = n()) %>%
  mutate(group = "100KGP Probands") %>%
  ungroup() %>%
  configure_ethnicity_labels()

## Number of people in England (ONS 2021) per sex, ethnicity, and age bin
ons_age_sex_ethn_matrix <- fread("../public_data/england_census_sex_age_ethn_2021.tsv") %>%
  mutate(age = ifelse(age == "100 or over", 100, age),
         age = as.numeric(age)) %>%
  replace(is.na(.), 0) %>%
  pivot_longer(!age,
               names_to = "ethn_sex",
               values_to = "num_inds") %>%
  mutate(ethn_sex = case_when(grepl("Male", ethn_sex, fixed = TRUE) ~ gsub(" Male", "_Male", ethn_sex),
                              grepl("Female", ethn_sex, fixed = TRUE) ~ gsub(" Female", "_Female", ethn_sex)),
         group = "ONS 2021") %>%
  separate(ethn_sex,
           into = c("participant_ethnic_category",
                    "participant_phenotypic_sex"),
           sep = "_", remove = TRUE,
           extra = "drop", fill = "right") %>%
  configure_ethnicity_labels()

#### Combined dataset
ons_rd_age_sex_ethn_matrix <- rbind(ons_age_sex_ethn_matrix,
                                    rd_age_sex_ethn_matrix) %>%
  mutate(age_bin = case_when(age = between(age, 0, 9) ~ "Age 0 to 9",
                             age = between(age, 10, 19) ~ "Age 10 to 19",
                             age = between(age, 20, 29) ~ "Age 20 to 29",
                             age = between(age, 30, 39) ~ "Age 30 to 39",
                             age = between(age, 40, 49) ~ "Age 40 to 49",
                             age = between(age, 50, 59) ~ "Age 50 to 59",
                             age = between(age, 60, 69) ~ "Age 60 to 69",
                             age >= 70 ~ "Age 70 or over")) %>%
  dplyr::select(-age) %>%
  group_by(age_bin, participant_ethnic_category,
           participant_phenotypic_sex, group) %>%
  summarise(num_inds = sum(num_inds)) %>%
  pivot_wider(names_from = "participant_ethnic_category", values_from = num_inds) %>%
  replace(is.na(.), 0) %>%
  pivot_longer(!c(age_bin, participant_phenotypic_sex, group),
               names_to = "participant_ethnic_category", values_to = "num_inds") %>%
  group_by(participant_ethnic_category, group) %>%
  mutate(proportion_sex = ifelse(participant_phenotypic_sex == "Female",
                                 (num_inds / sum(num_inds)) * 100,
                                 -1 * ((num_inds / sum(num_inds)) * 100)))

###########################################################################
###########################################################################
###                                                                     ###
###                           ANALYSIS RUN                              ###
###                                                                     ###
###########################################################################
###########################################################################

#### Run test and create table
# prop_test_df <- ons_rd_age_sex_ethn_matrix %>%
#  run_rd_vs_ons_prop_test()

#########################################################################
#########################################################################
###                                                                   ###
###                          WRITE DATA                               ###
###                                                                   ###
#########################################################################
#########################################################################

## Analysis data
data_dir <- paste0(analysis_dir, "ons_ethnicity_analysis/")

## Create directory.
dir.create(file.path(data_dir))
cat(paste0("Writing analyses to: ", data_dir, "\n"))

## ONS age sex ethnicity matrix
fwrite(ons_rd_age_sex_ethn_matrix,
       paste0(data_dir, "ons_100kGP_probands_age_sex_ethnicity.csv.gz"))

###########################################################################
###########################################################################
###                                                                     ###
###                         END OF SCRIPT                               ###
###                                                                     ###
###########################################################################
###########################################################################
