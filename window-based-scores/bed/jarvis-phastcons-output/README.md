zcat < mean_scores_by_window.bed.gz | tail -n+2 | cut -d',' -f2-7 | sed 's/,/\t/g' > mean_scores_by_window.bed
