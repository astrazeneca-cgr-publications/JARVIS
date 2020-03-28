# Get phastCons scores per chr
sbatch ./intersect_phastcons_per_chr.sh

# Integrate JARVIS, gwRVIS and phastCons
sbatch ./submit_win_grouping.sh

# Analyse most-intolerant/least-conserved or least-intolerant/most-conserved (GO, Pathway enrichment, etc.)
analyse_wins_by_metrics.py
