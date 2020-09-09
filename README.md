# JARVIS

- **JARVIS** ("Junk" Annotation genome-wide Residual Variation Intolerance Score): a comprehensive deep learning framework to prioritise non-coding variants in whole genomes, using human-lineage purifying selection features and primary sequence context.

- gwRVIS (genome-wide Residual Variation Intolerance Score): genome-wide intolerance to variation score


**Data availability**: genome-wide JARVIS and gwRVIS scores will be publicly available upon publication in a peer-reviewed journal.





## Installation instruction for JARVIS/gwRVIS modules
```
conda create -n jarvis python=3.7 r r-devtools r-tidyverse 
conda config --add channels bioconda
conda config --add channels conda-forge
conda activate jarvis  


conda install --file requirements.txt  
```

# R dependencies 
> install.packages(c("glm2", "glmnet", "lmridge", "plotmo", "pRoc", "ggplot2", "ggridges", "RColorBrewer")) 
> devtools::install_github("thomasp85/patchwork")
