# CURRENT:

• Add aggregate (avg.) metrics for sensitivity, precision, etc. in "jarvis/deep_learn_raw_seq/train_nn_model.py"

• Validate JARVIS on non-coding benign vs pathogenic variants from denovodb (when trained with all of ClinVar and not using denovodb as a control set)
	- Look at: 'benchmark_against_denovo_db.py'
	- Consider aggregating denovo-db and Clinvar pathogenic for the training with CV to see results!
 

• Use the same CV batches for benchmarking all other scores (if possible)

• Get metrics of Sensitivity (low number of FN) and Specifity (low number of FP) to see if I can use it safely on the whole genome without risking to predict too many False Positives.

• Re-run JARVIS training with Random Forest too, to extract feature importance and for benchmarking comparison.


• Submit jobs for jarvis training with clinvar/hgmd pathogenic sets and clinvar/denovodb(non_SSC)/topmed_uniq benign variants

• Train JARVIS on UTR and predict on intergenic


- Build JARVIS model across all non-coding regions - and compare against other scores trained on all same regions as well (ideally with same batches in CV!)


- Create sets of 'pathogenic'/'benign' based on their gERP++ scores or based on the CADD datasets



[!] Attention:
- cnn1_brnn1 on intergenic gave avg. AUC=0.5


Minor:
- Check https://github.com/slundberg/shap for Feature Importance from CNNs.
- Add a conservation score into Jarvis (when not trained based on gERP-labelling)
- Compare AF of variants that are both in ClinVar and HGMD and of those that are only in HGMD.
- Fix tidyverse installation: make_ggridges_plot.py uses gather() that requires this dependency - see if I can use some other function to discard this dependency
- (Maybe) add n-repeated k-fold CV for Jarvis (and clinvar-variant classifier) - https://stats.stackexchange.com/questions/82546/how-many-times-should-we-repeat-a-k-fold-cv
- Run grid-search for optimal nn (dnn / cnn) architecture selection -- need to create static splits of the batches for this task to be comparable.
• Default  pos_neg_ratio = 1/1; Add this to config file as a parameter
• Create clean instance of model for each fold in 'jarvis/variant_classification/classifiers.py' (--redundant as currently only Logistic Regression and Random Forest are to be used in this module)


DONE:
- Get random sample of control/benign variants that is up to x10 (?) larger than the pathogenic set, for the genomic class or set of genomic classes under consideration! [DONE]
- Combine structured and sequence data [DONE]
- Download Ryan s new control datasets and build new sets of control variants per score [DONE]
- test_and_evaluate_model -- need Debugging... - See what happens when using RandomForest for ClinVar classification, if I manage to predict successfully (check Sensitivity, Precision, etc.) and see what I can do for the DNN to not overfit the data. [DONE] -- FIXED! (with utr case): with pos/neg ratio = 1/1. 
- For performance comparison between structured, sequences and both get a good run from sequences or structured and use the exact same CV batches to assess performance! [DONE]

- Optimse cnn+dnn network for predictions with 'both' types of features ('structured' and 'sequences') [DONE]
- Make sure I fix imbalance everywhere  [DONE]
- Set RF as default classifier for JARVIS in 'jarvis/variant_classification/run_variant_classification.py' [DONE] -- the DNN (FeedF and CNNs are in the 'jarvis/deep_learn_raw_seq/' module
- Remove any overlapping benign variants from pathogenic file before classification [DONE]
