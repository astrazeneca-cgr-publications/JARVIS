library(pROC)
options(warn=-1)

args = commandArgs(trailingOnly=TRUE)
input_dir = args[1]

output_file = paste(input_dir, 'DL-delong-results.txt', sep='/')
sink(output_file)


calc_roc_for_score <- function(score){

	# y_labels
	y_lab = vector()
	res = readLines(paste(input_dir, '/', score, '.y_label_lists.txt', sep=''))

	cnt = 1
	for(line in res){
		y_lab = c(y_lab, as.numeric(strsplit(res[cnt], ", ")[[1]]))
		cnt = cnt + 1
	}

	# y_probas
	y_prob = vector()
	res = readLines(paste(input_dir, '/', score, '.y_proba_lists.txt', sep=''))

	
	cnt = 1
	for(line in res){
		y_prob = c(y_prob, as.numeric(strsplit(res[cnt], ", ")[[1]]))
		cnt = cnt + 1
	}

	cur_roc = roc(y_lab, y_prob, quiet=T)

	return(cur_roc)
}



base_scores = c('JARVIS-structured.intergenic_utr_lincrna_ucne_vista', 'JARVIS-sequences.intergenic_utr_lincrna_ucne_vista', 'JARVIS-both.intergenic_utr_lincrna_ucne_vista')
scores = c('ncER_10bp', 'cdts', 'linsight', 'gwrvis', 'jarvis', 'cadd', 'dann', 'phyloP46way', 'phastCons46way', 'orion')

#base_scores = c('gwrvis')
#scores = c('gwrvis', 'cadd', 'phyloP46way', 'phastCons46way', 'orion', 'gwrvis+cadd')


all_scores = unique(c(scores, base_scores))
roc_data = vector(mode="list", length=length(all_scores))
names(roc_data) = all_scores


idx = 1
for(score in all_scores){
	cur_roc = calc_roc_for_score(score)

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

#roc.test(roc_data[['ncER_10bp']], roc_data[['linsight']], method='delong')
#roc.test(roc_data[['ncER_10bp']], roc_data[['cadd']], method='delong')
#roc.test(roc_data[['jarvis']], roc_data[['cadd']], method='delong')
#roc.test(roc_data[['linsight']], roc_data[['jarvis']], method='delong')
