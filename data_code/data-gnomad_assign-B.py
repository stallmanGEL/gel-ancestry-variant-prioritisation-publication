#!/usr/bin env python3

##-------------------------------------------------------------------------
## Configs + packages

### SCRIPT ADAPTED FROM GNOMAD (https://github.com/broadinstitute/gnomad_qc/blob/main/gnomad_qc/example_notebooks/ancestry_classification_using_gnomad_rf.ipynb)

import onnx
import hail as hl
from gnomad.sample_qc.ancestry import apply_onnx_classification_model, assign_population_pcs

read_if_exists = True

## Config
exec(open('../config/config.py').read())

#########################################################################
#########################################################################
###                                                                   ###
###                       DATA & VARIABLES                            ###
###                                                                   ###
#########################################################################
#########################################################################

# Check that this is what was used in v3
v3_num_pcs = 20

# Check that this is what was used in v3
v3_min_prob = 0.8

# Plink files (from data-aggV2_gnomad_sites-A.sh)
plink_to_project = out_dir + "aggV2_gnomad_sites_data/aggV2_all_chr_gnomad_anc"

# Load ONNX model
with hl.hadoop_open(gnomad_v3_onnx_rf, "rb") as f:
    v3_onx_fit = onnx.load(f)

# Load loadings Hail table
v3_loading_ht = hl.read_table(gnomad_v3_loadings)

# Read MT and loadings
mt = hl.import_plink(bed=plink_to_project + ".bed",
                     bim=plink_to_project + ".bim",
                     fam=plink_to_project + ".fam",
                     reference_genome="GRCh38")

#########################################################################
#########################################################################
###                                                                   ###
###                          RUN ANALYSIS                             ###
###                                                                   ###
#########################################################################
#########################################################################

mt = mt.checkpoint(
    plink_to_project + ".mt", 
    overwrite=not read_if_exists, 
    _read_if_exists=read_if_exists
)

# Project new genotypes onto loadings.
v3_pcs_ht = hl.experimental.pc_project(
    mt.GT, v3_loading_ht.loadings, v3_loading_ht.pca_af,
)

# Checkpoint PC projection results.
v3_pcs_ht = v3_pcs_ht.checkpoint(
    plink_to_project + "_pca.ht", 
    overwrite=not read_if_exists, 
    _read_if_exists=read_if_exists
)

ht, model = assign_population_pcs(
    v3_pcs_ht,
    pc_cols=v3_pcs_ht.scores[:v3_num_pcs],
    fit=v3_onx_fit,
    min_prob=v3_min_prob,
    apply_model_func = apply_onnx_classification_model,
)

#########################################################################
#########################################################################
###                                                                   ###
###                          DATA EXPORT                              ###
###                                                                   ###
#########################################################################
#########################################################################

ht.export(out_dir + "ancestry_data/gnomad_assignments.tsv")

###########################################################################
###########################################################################
###                                                                     ###
###                         END OF SCRIPT                               ###
###                                                                     ###
###########################################################################
###########################################################################
