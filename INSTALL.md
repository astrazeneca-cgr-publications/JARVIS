Installation instruction for gwRVIS/JARVIS
------------------------------------------

conda create -n jarvis python=3.7 r r-devtools r-tidyverse 
conda config --add channels bioconda
conda config --add channels conda-forge
conda activate jarvis  


conda install --file requirements.txt  

R: > install.packages(c("glm2", "glmnet", "lmridge", "plotmo", "pRoc", "ggplot2", "ggridges", "RColorBrewer")) 
> devtools::install_github("thomasp85/patchwork")
