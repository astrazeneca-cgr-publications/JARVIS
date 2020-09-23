# JARVIS

## Overview
- **JARVIS** ("Junk" Annotation genome-wide Residual Variation Intolerance Score): a comprehensive deep learning framework to prioritise non-coding variants in whole genomes, using human-lineage purifying selection features and primary sequence context.

- gwRVIS (genome-wide Residual Variation Intolerance Score): genome-wide intolerance to variation score


<br><br>



## Installation instruction for JARVIS/gwRVIS modules

- Python dependencies
```
conda create -n jarvis python=3.7 r r-devtools r-tidyverse 
conda config --add channels bioconda
conda config --add channels conda-forge
conda activate jarvis  


conda install --file requirements.txt  
```

- R dependencies 
> install.packages(c("glm2", "glmnet", "lmridge", "plotmo", "pRoc", "ggplot2", "ggridges", "RColorBrewer")) 
> devtools::install_github("thomasp85/patchwork")


<br><br>

## Run
- Instructions to generate the JARVIS and gwRVIS scores are available in the [README.md](modules/README.md) file within `modules`.
- Subsequent sub-folders may also contain their own README files with instructions to run them independently or for ad-hoc analyses.
- Other folders and theri sub-folders (such as `ensembl/`, `gnomad/` and `other_datasets/`) are accompanied with README files and scripts to download and pre-process any other required datafiles that are not available in the JARVIS GitHub repositoy.


<br><br>

## Data availability
genome-wide **JARVIS** and **gwRVIS** scores will be publicly available upon publication in a peer-reviewed journal.


<br><br>
