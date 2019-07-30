# TODO LIST:
> Run logistic regression between intergenic and ucnes (or other classes) to see if gwRVIS can distinguish them efficiently.

-- Add CV ROC curves to "scores_benchmarking/run_clinvar_benchmarking.py":
Follow code from 'https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html'





> Try just a few more versions of the logistic regression for gwRVIS calculations, using bins and/or other features (how about negative binomial regression?)

> Add annotation for CpG islands using pre-built Ensembl tracks for hg19.

> Add methylation data for each window. Think about other data types too.

> Add CNNs and RNNs with Keras functional API


-----------------
# DONE:
- Make sure mut_rate, cpg etc. are calculated for the very last window too.
