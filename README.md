# gel-ancestry-variant-prioritisation-publication
Code required to generate data, replicate analysis and generate plots from the paper **Missing genetic diversity impacts variant prioritisation for rare disorders**.

LINK TBC

## Repository structure
The repository is structured into 5 folders. Further information can be found in the relevant `README` files in each directory.
1. `public_data` contains publicly available data files relevant for certain analyses.
2. `config` contains configuration files and paths to flat files in the **Genomics England Research Environment**.
3. `data_code` contains code for generating data used for running analysis. **NOTE: Certain files required to run these scripts exactly as in the paper may be unavailable to external researchers (see below).**
4. `analysis_code` contains code for running analyses. Uses files output from `data_code`.
5. `plotting_code` contains code for creating plots found in the paper and Supplement. Uses files output from `analysis_code` and `data_code`.

___

## Genomics England Research Environment
To run code in this repository, and to gain access to the relevant datasets and files as specified in `config` you must first gain access to the **Genomics England Research Environment**. Information on how to gain access can be found in the [Genomics England Research Environment User Guide](https://re-docs.genomicsengland.co.uk/welcome/).  

Once logged into the **Genomics England Research Environment**, this repo can be found in the directory: `TBC` along with all files generated from the scripts found in `data_code` and `analysis_code`.

**NOTE:** To run these analyses, please ssh into the Genomics England [Double Helix HPC](https://re-docs.genomicsengland.co.uk/hpc/). Please see the relevant documentation for guidance on job submission on Double Helix.

---

## Requirements

### R
All `.R` code was run using `R 4.3.3`. This can be loaded on the [Double Helix HPC](https://re-docs.genomicsengland.co.uk/hpc/) using `module load R/4.3.3`.   

Packages required for running all `.R` code in this repo can be found in `r_requirements.txt`. Please see the relevant documentation on how to [Work with R packages on Double Helix](https://re-docs.genomicsengland.co.uk/r_packages/).

### Python
All `.py` code was run using `Python 3.9.0` via `conda`. 

Packages required for running `data_code/data-gnomad_assign-B.py` can be found in `python_requirements.txt`. Please see the relevant documentation on how to [Work with Python packages and personal conda environments on Double Helix](https://re-docs.genomicsengland.co.uk/hpc_conda/).

---

## Data availability notes for external researchers

If you are an external researcher with access to the **Genomics England Research Environment** and wish to run certain scripts in the `data_code` folder in this respository, please note that certain datasets and files are currently **_unavailable_** to external researchers.   

Note, however, files generated from these datasets are available in the relevant subfolders of the `TBC` directory in the **Genomics England Research Environment**, enabling full replication of the results. Furthermore, permission arguments can be set so that all code in the `data_code` can be run without access to these original datasets (however, outputs will differ from that in the paper.)

#### The COVID-19 Cohort Dataset

Data from the [COVID-19 Genomics Cohort](https://re-docs.genomicsengland.co.uk/covid5/) was used as per the _Identifying cPAVs classified as ultra-rare that are common in a diverse reference database from the UK population_ secton in the paper.  

This dataset is currently not accessible to researchers on the **Genomics England Research Environment** and are instead hosted for research use on [CloudOS](https://re-docs.genomicsengland.co.uk/cloudos/). 

Any code that requires access to this dataset (as specified in the `data_code` directory `README`) will run without them unless a `COVID_ACCESS` command line flag is added after the relevant scripts, which will then require that you have permissions to access to the relevant files.

#### PrimateAI-3D scores

`PrimateAI-3D` scores require a [Signed academic license agreement](https://primateai3d.basespace.illumina.com/download#:~:text=Please%20click%20here%20to%20accept%20the%20academic%20license%20agreement) to access, and are therefore not universally available in the **Genomics England Research Environment**.  

Any code that requires access to these scores (as specified in the `data_code` directory `README`) will run without them unless a `PAI3D_ACCESS` command line flag is added after the relevant scripts, which will then require that you have permissions to access to the relevant files.

---

## Contact

Please contact samuel.tallman@genomicsengland.co.uk for any further information regarding code usage.
