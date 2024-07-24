# `data_code` directory
Code required to generate data used in the paper **Missing genetic diversity impacts variant prioritisation for rare disorders** as input for the scripts in the `analysis_code` directory and/or the `plotting_code` directory.

## Code structure
Scripts are written to be run sequentially.   
**A.** `data-labkey-A.R`, `data-ukbb_assign.R`, `data-aggV2_gnomad_sites-A.sh`, and `data-missing_gnomad_pavs-A.sh` must be run first.  
**B.** Followed by `data-delet_preds-B.sh`,  `data-gnomad_afs-B.sh`, `data-gnomad_assign-B.py`, `data-panel_app_meta-B.R`, `data-roh-B.R`  
**C.** Followed by `data-combine_filter-C.R`.  
**D.** And finally `data-relate-D.sh`.

___

## Output files
The output files of all scripts in this directory (as specified in _Code Description_) will appear in the output directory specified in the relevant configuration files in the `config` directory.

- For `.R` scripts: `out_dir` as specifified in `config/config.R`
- For `.py` scripts: `out_dir` as specified in `config/config.py`
- For `.sh` scripts: `OUT_DIR` as specified in `config/config.sh`

These files can also be found in the `TBC` directory in the **Genomics England Research Environment**.

---

## Code Description

### **Group A**

```
Rscript data-labkey-A.R
```
Will pull relevant meta data, tiering info and diagnostic info for 100,000 Genomes Project (100kGP) participants from `LabKey` data tables.   
**Output files:** 
- `labkey_data/meta_data.csv.gz`
- `labkey_data/tiering_data.csv.gz`
- `labkey_data/phenotypes_data.csv.gz`
- `labkey_data/exit_data.csv.gz`
- `labkey_data/tiering_frequency_data.csv.gz`
- `labkey_data/panels_applied_data.csv.gz`
- `labkey_data/tiering_sorted_variants.tsv`

```
Rscript data-ukbb_assign-A.R
```  

Will use the the `hq_snps_plink` (.bed, .bim, .fam files for all aggV2  / aggV2 and aggCOVIDv5 participants) alongside `prive_ukbb_build38_afs` and `prive_ukbb_build38_loadings` (reference data from the UK Biobank generated from [Prive et al. 2021]()) specified in `config/config.R` to assign individuals to genetically-inferred ancestry (GIA) groups.    
**Output files:** 
- `ancestry_data/ukbb_assignments.csv`

```
bash data-aggV2_gnomad_sites-A.sh
```  
Will loop through aggV2 `plink2` .pgen files specified in `config/aggV2.file_list` and subset to `GNOMAD_ANC_ALLELES` (TSV file of sites provided by gnomADv3 used for reference population assignment) as specified in `config/config.sh` then merge into a single .bed, .bim, .fam fileset.   
**Output files:** 
- `aggV2_gnomad_sites_data/aggV2_all_chr_gnomad_anc.bed`
- `aggV2_gnomad_sites_data/aggV2_all_chr_gnomad_anc.bim`
- `aggV2_gnomad_sites_data/aggV2_all_chr_gnomad_anc.fam`

```
bash data-missing_gnomad_pavs-A.sh aggV2
```  
Will loop through aggV2 files specified in `config/aggV2.file_list`, extract PAVs based on `VEP 105` annotations, and then compare that with variants in `GNOMAD_V3_VCF` and `GNOMAD_V2_VCF` (gnomAD VCF files) as specified in `config/config.sh`. All variants that do not appear in either are then counted per individual.   
**Output files:** 
- `pavs_gnomad_miss_data/aggV2_all_gnomAD_missing_pavs.scount.gz`

### **Group B**

```
bash data-delet_preds-B.sh
```  

Will use the positions specified the file `labkey_data/tiering_sorted_variants.tsv` and loop through the `VEP 105` annotated VCFs and pull `CADD` scores at those specified sites. It will also extract `AlphaMissense` scores from the `AM_SCORES` file specified in `config/config.sh`.    
**Output files:** 
- `deleterious_predictions_data/cadd_scores_tiered.tsv.gz`
- `deleterious_predictions_data/alpha_missense_scores_tiered.tsv.gz`

```
bash data-gnomad_afs-B.sh
```  
Will use the positions specified the file `labkey_data/tiering_sorted_variants.tsv` and extract per population AN/ACs at these sites from `GNOMAD_V3_VCF` and `GNOMAD_V2_VCF`.    
**Output files:** 
- `gnomad_afs_data/gnomad_v2_exomes_tiering.afs.gz`
- `gnomad_afs_data/gnomad_v3_genomes_tiering.afs.gz`

```
python data-gnomad_assign-B.py
```  
Will use `aggV2_gnomad_sites_data/aggV2_all_chr_gnomad_anc.bed`, `aggV2_gnomad_sites_data/aggV2_all_chr_gnomad_anc.bim`, and `aggV2_gnomad_sites_data/aggV2_all_chr_gnomad_anc.fam` files alongside `gnomad_v3_loadings` and `gnomad_v3_onnx_rf` (PC loadings and ONNX random forest ancestry assignment models from gnomADv3) as specified in `config/config.py` to assign all 100kGP individuals to ancestry groups.    
**Output files:** 
- `ancestry_data/gnomad_assignments.tsv`

```
Rscript data-panel_app_meta-B.R
```  

Will use the per participant `PanelApp` panels applied as specified in `labkey_data/panels_applied_data.csv.gz` alongside the new gene panels found in `../public_data/panel_app_stats_data.csv.gz` (shared alongside this repo) and calculate some meta-data for each participant.    
**Output files:** 
- `panels_app_data/panel_app_meta_data.csv.gz`

```
Rscript data-roh-B.R
```  

Will use file paths specified in `labkey_data/meta_data.csv.gz` and then generate some Runs of Homozygosity (ROH) statistics based on bed files available for the majority of 100kGP participants.     
**Output files:** 
- `roh_data/roh_data.csv.gz`

### **Group C**

```
Rscript data-combine_filter-C.R
```  

Will take many of the above files and combine and/or filter them in various ways to be used for the files in `analysis_code` (as well as scripts in **Group D**).  
**Output files:**
- `combined_data/combined_meta_data.csv.gz`
- `combined_data/combined_tiering_frequency_data.csv.gz`
- `combined_data/combined_predictions_data.csv.gz`
- `combined_data/combined_missing_pavs_gnomad.csv.gz`
- `combined_data/filtered_retiering_probands_data.csv.gz`
- `combined_data/relate_ss_poplabels.tsv`

### **Group D**

```
bash data-relate-D.sh
```  

Will take a variety of files specified in `config/config.sh` (e.g. phased VCFs, accessibility masks, genetic maps) alongside the individuals specified in `combined_data/relate_ss_poplabels.tsv` and then estimate within and between group coalescence rates from ARGs generated from these files.  
**Output files:**
- `relate_data/phased_panel_relate_subset_anc.pairwise.coal`

Note: This will take some time. Parallelising across chromosomes is advised if running (although the script shown loops through chromosomes).

____

## Code requiring elevated permissions
Certain code to generate the files used in the paper can be run using the above scripts with elevated access priveleges (see `../README.md`). If you have relevant access, these files can be generated in the following way.

**Group A**

```
bash data-missing_gnomad_pavs-A.sh aggCOVID
```  
Assumes you have access to the files as specified in `config/aggCOVID_v5.file_list`. If so...  
**Output files:**
- `pavs_gnomad_miss_data/aggCOVID_all_gnomAD_missing_pavs.scount.gz`

**Group B**

```
bash data-delet_preds-B.sh PAI3D_ACCESS
``` 

Assumes you have access to the `PAI3D_SCORES` file as specified in `config/config.sh`. If so, the following files will also be generated:
- `deleterious_predictions_data/primate_ai_3d_scores_tiered.tsv.gz`

**Group C**

```
Rscript data-combine_filter-C.R PAI3D_ACCESS COVID_ACCESS
``` 

Assumes that you have generated both of the above files (`pavs_gnomad_miss_data/aggCOVID_all_gnomAD_missing_pavs.scount.gz` and `deleterious_predictions_data/primate_ai_3d_scores_tiered.tsv.gz`), and hence will incorporate them into the output files (`combined_data/combined_missing_pavs_gnomad.csv.gz`, `combined_data/combined_predictions_data.csv.gz`) respectively. It will also generate the file:
- `combined_data/aggCOVID_unrelated_pops.fam`


**Group D**

```
bash covid_afs-D.sh
```

This file takes the `combined_data/aggCOVID_unrelated_pops.fam` and loops through the .pgen files in `config/aggCOVID_v5.file_list` to generating per GIA group Allele Counts and Allele Numbers, then combines the output. It also takes the sites in `labkey_data/tiering_sorted_variants.tsv` and subsets this to only the sites in this file.  
**Output:**
- `covid_afs_data/aggCOVID_full.acount.gz`
- `covid_afs_data/aggCOVID_tiering.acount.gz`