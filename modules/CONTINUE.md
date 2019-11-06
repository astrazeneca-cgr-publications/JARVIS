*** FINAL TO_DO ***

- Conservation (phastcons46way primates -- extract most/least conserved nt) 
	 tabix /projects/cgr/users/kclc950/gwRVIS/other_datasets/conservation/phastCons46way_primates/bed/chr16.phastCons46way.primates.high_conf_regions.bed.gz -B intergenic.full_table.bed  
	 Use bed files in gwrvis_distribution (e.g. gwrvis.ucne.mutually_excl.bed) to subset most/least conserved regions by genomic class
	- Label positive/negative data points as those that are most-conserved_&_most-intolerant vs least-conserved_&_most-tolerant


- Convert score annotations to hg38 .......


- Run original Topmed (without liftover) with most recent annotation (hg38) and see if I get similar clinvar classifications with gnomad-3


- New gnomad version 3 (use updated annotations - liftunder probably not an option)

- topmed (? Try gnomad 3.0 first and see if I ll stick with that eventually): try with hgmd/denovo/topmed_uniq; Then run jarvis for each of these too


- Sliding window (during Review)


# "CURRENT":
- Check status of: "./submit_all_jarvis_jobs.sh conf/config.yaml 1" and re-plot performance metrics with "jarvis/performance_estimation/process_performance_metrics.py"  - Introns failed (now re-running for structured and then sequences and both)


- In ML benchmarking, remove NAs instead of imputing them with 0.5 (and count how many were discarded)

- Create sets of 'most/least conserved' (as a proxy for 'pathogenic'/'benign') based on their gERP++ (or other) scores or based on the CADD datasets
	Current --> Sort bed files by conservation (phastcons46way in pimates)
	- Intersect with "full_gwrvis_and_regulatory_features.All_genomic_classes.tsv"


- Play with MAF when using the new gnomAD release (r3.0)



• Validate JARVIS on non-coding benign vs pathogenic variants from denovodb (when trained with all of ClinVar and not using denovodb as a control set)
	- Look at: 'benchmark_against_denovo_db.py' [DONE - but sample size is small and maybe not the best case for benchmarking]
	- Consider aggregating denovo-db and Clinvar pathogenic for the training with CV to see results!
 

• Train JARVIS-sequences using 3kb around the variant instead of the gwRVIS 3kb window



• Submit jobs for jarvis training with clinvar/hgmd pathogenic sets and clinvar/denovodb(non_SSC)/topmed_uniq benign variants

• Train JARVIS on UTR and predict on intergenic





- Run everything again with TOPMED this time instead of gnomAD


[!] Attention:
- cnn1_brnn1 on intergenic gave avg. AUC=0.5
	- Play a bit more with some LSTM architectures (and probably better with a much smaller window size)



Minor:
- Check https://github.com/slundberg/shap for Feature Importance from CNNs.
- Add a conservation score into Jarvis (when not trained based on gERP-labelling)
- Compare AF of variants that are both in ClinVar and HGMD and of those that are only in HGMD.
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
• Re-run JARVIS training with Random Forest too, to extract feature importance and for benchmarking comparison. [DONE]
• Get metrics of Sensitivity (low number of FN) and Specifity (low number of FP) to see if I can use it safely on the whole genome without risking to predict too many False Positives. [DONE]
• Use the same CV batches for benchmarking all other scores (if possible) [DONE] -- Not possible, as not all scores will be defined for the same types of variants. So, I will just train them with independent CV batches, that will represent their generalised predictive power.

• Complete: Add metrics for sensitivity, precision, etc. in "jarvis/deep_learn_raw_seq/train_nn_model.py" [DONE]
	- Also test n-repeated in CV train_nn_model.py [DONE]
• Add aggr. (avg.) metrics for sensitivity, precision, etc. in "variant_classification/run_variant_classification.py" [DONE]
- Combine metrics from the two modules and plot them. [DONE]
- Fix tidyverse installation: make_ggridges_plot.py uses gather() that requires this dependency - see if I can use some other function to discard this dependency [DONE]
- (Maybe) add n-repeated k-fold CV for Jarvis (and clinvar-variant classifier) - https://stats.stackexchange.com/questions/82546/how-many-times-should-we-repeat-a-k-fold-cv [DONE]
- Add DANN and DANQ in final benchmarkings [DONE: DANQ does not provide their scores...] 
- Need to fix submit_all_jarvis_jobs.sh: I need to run for structured first only (get new fixed CV batch) and then in another run for sequences and both. [DONE]
- Build JARVIS model across all non-coding regions - and compare against other scores trained on all same regions as well (ideally with same batches in CV!) [DONE]
- Sensitivity analysis with all classes [DONE - optimal win: 3kb] 
- Create hg38/ annotation based on hg19 [DONE]
- Change all hg19 hardcoded references... [DONE]
- Update Ensembl annotations for chromatin structure, methylation, etc. (Look for 'Monocytes_CD14plus' in scripts) [DONE]
- Convert variant annotations to hg38 ....... [DONE]
