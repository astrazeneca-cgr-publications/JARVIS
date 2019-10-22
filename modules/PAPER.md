Add:

- Results from score bencmarking, using either clinvar or hgmd as pathogenic, and clinvar, denovodb or topmed for benign! Already gives some very good resutls, either for gwRVIS alone or in combination with CADD, justifying the existence of a lot more information.

	-- Check also results with Random Forest (instead of DNN) [DONE: Jarvis with Neural Nets performs better both for UTRs and intergenic, so can also show it and show the feature importance; gwRVIS is on top]

- Results jarvis/variant_classification/run_variant_classification.py: OK (for intergeneic mainly, and a bit for UTRs)

- CNN2_FC2 performs the best (CNN3_FC2 is slightly worse)

- JARVIS for UTR with structured already ranks pretty highly. Hopefully, by combining both structured and sequences it will exceed CADD performace too. -- It reaches almost the same level; can try a bit deeper CNN

- Check results at 'denovodb_benchmarking-ssc/ .. /intergenic.all_scores.pdf'
