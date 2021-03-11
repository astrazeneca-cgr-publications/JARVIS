# JARVIS

## Overview
- **JARVIS** ("Junk" Annotation genome-wide Residual Variation Intolerance Score): a comprehensive deep learning framework to prioritise non-coding variants in whole genomes, using human-lineage purifying selection features and primary sequence context.

- gwRVIS (genome-wide Residual Variation Intolerance Score): genome-wide intolerance to variation score

|Publication: |
| :---- |
|[Prioritizing non-coding regions based on human genomic constraint and sequence context with deep learning](https://www.nature.com/articles/s41467-021-21790-4). <br/>
Vitsios et al., __Nature Communications__, March 8, 2021 https://doi.org/10.1038/s41467-021-21790-4  |

<br>

![](misc/JARVIS-DNN-network.jpg?raw=true)

<br>



## Installation instructions for JARVIS/gwRVIS modules

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
- Other folders and their sub-folders (such as `ensembl/`, `gnomad/` and `other_datasets/`) are accompanied with README files and scripts to download and pre-process any other required datafiles that are not available in the JARVIS GitHub repositoy.


<br><br>

## Data availability
**JARVIS** and **gwRVIS** scores, across the whole genome, are publicly available at the following location:
[http://jarvis.public.cgr.astrazeneca.com](http://jarvis.public.cgr.astrazeneca.com)

All scores have been generated based on the **hg19** human assembly version.

<br><br>
