# `analysis_code` directory
Code required to perform analyses in the paper **Missing genetic diversity impacts variant prioritisation for rare disorders** as input for the `plotting_code/plots.R` script.

Code can be run  in any order, and assumes all files created and found in the **Genomics England Research Environment** as specified in `data_code/` exist (specifically those mentioned below).

## Output files
The output files of all scripts in this directory (as specified in _Code Description_) will appear in the `analysis_dir` directory specified in `config/config.R`. 

These files can also be found in the `TBC` directory in the **Genomics England Research Environment**.

---

## Code Description
```
Rscript analysis-ons_ethnicity.R
```
Will look for the following files created as specified in `data_code/`:  
- `combined_data/combined_meta_data.csv.gz`

As well as the file:
- `../public_data/england_census_sex_age_ethn_2021.tsv`

Will then compare ethnicity proportions between these two files.  
**Output:**

- `ons_ethnicity_analysis/ons_100kGP_probands_age_sex_ethnicity.csv.gz`

```
Rscript analysis-pca.R
```
Will look for the following files created as specified in `data_code/`:  
- `combined_data/combined_meta_data.csv.gz`

As well as the `kinship_matrix` and `hq_snps_plink` files as specified in `config/config.R`. Will perform PCA on unrelated individuals and then project relatives onto this PC space. Will then perform UMAP with top 16 PCs.   
**Output:**

- `pca_analysis/pca_umap.csv`

```
Rscript analysis-glm.R
```
Will look for the following files created as specified in `data_code/`:  
- `combined_data/combined_meta_data.csv.gz`
- `combined_data/combined_tiering_frequency_data.csv.gz`
- `combined_data/combined_predictions_data.csv.gz`
- `combined_data/filtered_retiering_probands_data.csv.gz`
- `labkey_data/exit_data.csv.gz`

Will then run a variety of Generalised Linear Models (GLMs) as shown in the paper and shown in Supplementary Tables. Statistics for these models (and other files for plotting) will then be found in:

- `glm_analysis/`

```
Rscript analysis-covid_filter.R
```
Will look for the following files created as specified in `data_code/`
- `covid_afs_data/aggCOVID_tiering.acount.gz`
- `combined_data/combined_tiering_frequency_data.csv.gz`
- `combined_data/filtered_retiering_probands_data.csv.gz`
- `labkey_data/exit_data.csv.gz`
- `combined_data/combined_meta_data.csv.gz`

Will then generate statistics based on cPAVs/GPcPAVs that are common in different GIA groups (but e.g. ultra-rare elsewhere).   
**Output:**
- `covid_analysis/aggCOVID_faf95_filter_tiering.csv`