library(pROC)
options(warn=-1)

args = commandArgs(trailingOnly=TRUE)
input_dir = '/projects/cgr/users/kclc950/JARVIS/other_datasets/structural-variants_single-nt/delong_test'




calc_roc_for_score <- function(score, pathogenic_class){

	# y_labels
	y_lab = vector()
	res = readLines(paste(input_dir, '/SV-', score, '.', pathogenic_class, '.y_label_lists.txt', sep=''))

	cnt = 1
	for(line in res){
		y_lab = c(y_lab, as.numeric(strsplit(res[cnt], ", ")[[1]]))
		cnt = cnt + 1
	}

	# y_probas
	y_prob = vector()
	res = readLines(paste(input_dir, '/SV-', score, '.', pathogenic_class, '.y_proba_lists.txt', sep=''))

	
	cnt = 1
	for(line in res){
		y_prob = c(y_prob, as.numeric(strsplit(res[cnt], ", ")[[1]]))
		cnt = cnt + 1
	}

	cur_roc = roc(y_lab, y_prob, quiet=T)

	return(cur_roc)
}




# ================   MAIN   ================
base_scores = c('jarvis')
scores = c('ncER_10bp', 'linsight', 'gwrvis', 'jarvis', 'cadd', 'eigen', 'orion')
all_scores = unique(c(scores, base_scores))

pathogenic_classes = c('promoter', 'copy_gain', 'lof', 'dup_lof', 'inv_span', 'dup_partial', 'utr', 'intronic')


# Iterate through each pathogenic class
for(pathogenic_class in pathogenic_classes){

	output_file = paste(input_dir, '/', pathogenic_class, '.DL-delong-results.txt', sep='')
	sink(output_file)

	roc_data = vector(mode="list", length=length(all_scores))
	names(roc_data) = all_scores


	idx = 1
	for(score in all_scores){

		#print(score)

		cur_roc = calc_roc_for_score(score, pathogenic_class)

		roc_data[[idx]] = cur_roc
		#print(roc_data[[idx]])
		idx = idx + 1
	}



	for(score1 in base_scores){
		cat(paste('\n\n\n>>> [Base score]:', score1))
		for(score2 in scores){
			if(score2 == score1){ next }

			cur_pval = ''

			cat(paste('\n', score2, '\n'))	
			
			# check if score2's auc is less than score1's auc
			res = roc.test(roc_data[[score2]], roc_data[[score1]],
			      quiet=T, method='delong', alternative='less')
			print(paste(score2, ' < ', score1, ' : ', res$p.value, sep=''))

			# check if score1's auc is less than score2's auc
			res = roc.test(roc_data[[score1]], roc_data[[score2]],
			      quiet=T, method='delong', alternative="less")
			print(paste(score1, ' < ', score2, ' : ', res$p.value, sep=''))

			# two-sided
			#res = roc.test(roc_data[[score1]], roc_data[[score2]], 
			#      quiet=T, method='delong')
			

		}
	}

}
