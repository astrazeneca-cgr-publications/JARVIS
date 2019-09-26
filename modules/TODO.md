# TODO:
# 'Major'
> Add CNNs and RNNs with Keras functional API.
    See from 'deep_learn_module/':
        - prepare_data.py
	- functional_nn_models.py
	- train_nn_modely.py


> Train on UTRs -> Predict on Intergenic: to justify if I can use all UTR variants for the JARVIS training to then present it as a non-coding variant prioritisation score for intergenic regions too.
> Combine HGMD and ClinVar into a single dataset



# 'Minor'
> Try just a few more versions of the logistic regression for gwRVIS calculations, using bins and/or other features (how about negative binomial regression?)

> Add annotation for CpG islands using pre-built Ensembl tracks for hg19.

> Add methylation data for each window. Think about other data types too.


# 'NOTES'
- Currently using ClinVar variants, allowing multiple alleles at the same loci and counting these as distinct mutation events ----

-----------------
# DONE:
- Make sure mut_rate, cpg etc. are calculated for the very last window too. [DONE]
- Run logistic regression between intergenic and ucnes (or other classes) to see if gwRVIS can distinguish them efficiently. [DONE]
- Add CV ROC curves to "scores_benchmarking/run_clinvar_benchmarking.py":
Follow code from 'https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html' [DONE]
- Jarvis-variant_classification cannot read hgmd as input at the moment (seems to be fixed for clinvar pathogenic input) [DONE] - Fixed
- Compare pathogenic variants sets between ClinVar and HGMD (see overlap, etc.) [DONE]: They are highly overlapping
