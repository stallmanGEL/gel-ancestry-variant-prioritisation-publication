### Static wildcards (THESE CAN BE CHANGED)
OUT_DIR=/pgen_int_work/BRS/dd/sam/data/gel-ancestry-variant-prioritisation-publication/data
ANALYSIS_DIR=/pgen_int_work/BRS/dd/sam/data/gel-ancestry-variant-prioritisation-publication/analysis
TMP_DIR=/re_scratch/dd/sam/temp

### GEL data resources
PHASED_AGGV2_VCF=/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/phased_data/genotypes/phased_panel_chr${CHROM}.vcf.gz

### Public data resources
AM_SCORES=/public_data_resources/AlphaMissense/AlphaMissense_hg38.tsv.gz
GNOMAD_V2_VCF=/public_data_resources/gnomad/2.1.1/grch38_lift_over/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
GNOMAD_V3_VCF=/public_data_resources/gnomad/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.${CHROM}.vcf.bgz
GENETIC_MAP=/public_data_resources/recomb_hg38/genetic_map_GRCh38_chr${CHROM}.txt
ANCESTRAL_GENOME=/public_data_resources/homo_sapiens_ancestor_hg38/homo_sapiens_ancestor_${CHROM}.fa

#### These scores are not publicly available in the Genomics England Research Environment
PAI3D_SCORES=/pgen_int_work/BRS/dd/sam/data/PrimateAI-3D_scores.csv.gz

### Research Registry data folder.
PATH_TO_RELATE=/pgen_int_work/BRS/dd/sam/software/relate/
TGP_MASK=/pgen_int_work/BRS/dd/sam/software/accessibility_mask/tgp_pilot_mask_GRCh38_chr${CHROM}.fasta.gz
GNOMAD_ANC_ALLELES=/pgen_int_work/BRS/dd/data/gnomad_ancestry_random_forest/v3/gnomad_v3_ancestry_alleles.tsv