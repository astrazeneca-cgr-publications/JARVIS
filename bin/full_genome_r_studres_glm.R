library(glm2)
library(MASS)
library(glmnet)
library(lmridge)
library(plotmo)
library(pROC)
library(e1071)


args = commandArgs(trailingOnly=T)

out_dir = args[1]
filter_outliers_before_regression = toupper(args[2])
all_variants_upper_thres = as.numeric(args[3])
win_len = as.numeric(args[4])
kmer = as.numeric(args[5])
chr_type = args[6]


gwrvis_dir = paste(out_dir, 'gwrvis_scores', sep='/')
input_path = paste(out_dir, "/tmp/", sep='')
files = NULL
if( chr_type == 'autosomal'){
	files = list.files(path = input_path, pattern = "^Xy.")
} else if( chr_type == 'sex'){
	files = list.files(path = input_path, pattern = "^Xy.chrX")
}
print(files)

total_entries_per_chr = list()
total_df = data.frame()
all_chrs = vector()

for(file in files){

	chr = file
	chr = gsub("Xy.chr", '', chr)
	chr = gsub(".txt", '', chr)

	# process X chromosome separately if input chr_type == autosomal
	if(chr == 'X' & chr_type == 'autosomal'){ next }

	all_chrs = c(all_chrs, chr)
	
	
	df = read.table(paste(input_path, file, sep='/'), sep=',', header=T)
	rownames(df) = paste('chr', chr, '_', df$idx, sep='')
	df$idx = NULL
	print(head(df, 20))
	total_entries_per_chr[chr] = nrow(df)

	if(nrow(total_df) > 0){
		total_df = rbind(total_df, df)
	} else{
		total_df = df
	}
	
	print(nrow(total_df))
}


# Remove windows that had no variant data in the original VCF file
placeholder_val = -1
total_df = total_df[ total_df$mut_rate != placeholder_val, ]
print(nrow(total_df))


# ========== Auxiliary functions ==========
is_outlier <- function(vec, z_thres=3.5){          
	#vec = vec[ !is.nan(vec) ]         
	median = median(na.omit(vec))        
	#print(paste('median:', median))          

	diff = sqrt((vec - median)**2)         
	med_abs_deviation = median(na.omit(diff))         
	modified_z_score = 0.6745 * diff / med_abs_deviation         

	#print(head(modified_z_score))          
	return(modified_z_score > z_thres) 
} 

convert_pval_to_phred <- function(pval){
	phred = -10 * log10(as.numeric(pval))
	return(phred)
}

scale_and_center <- function(x){
	x = as.numeric(x)
	return( (x - mean(x)) / sd(x) ) 
}
# =========================================



# [Currently flag is False: probably not use it in production too] exclude outliers before fitting regression model
if(as.logical(filter_outliers_before_regression)){
	print('filtering before regression...')
	total_df = total_df[ total_df$all_variants <= all_variants_upper_thres, ]
}

## ---------------------------- Original ----------------------------
regr_model = glm(formula = y ~ all_variants, data=total_df) 

# Get RMSE (Root Mean Square Error)
error = residuals(regr_model)
RMSE = sqrt(mean(error^2))
print(paste('RMSE:', RMSE))

# Get studentised residuals
stud_res = studres(regr_model)

total_df = cbind(total_df, 'stud_res' = as.numeric(stud_res))
print(head(total_df))



# ===== BETA ====== 
# (Calulate Watterson estimator: Ne is the effective population size and μ is the per-generation mutation rate of the population of interest)
# Ne = 10000
# total_df$theta = 4 * Ne * total_df$mut_rate
# =================

## ===== BETA ===== (scale total number of variants by mutation rate)
# final_df$all_variants = final_df$all_variants / final_df$mut_rate
# =============================================================

# ===== BETA ===== 
# Ridge regression (L2 regularization)
# y = total_df$y
# X = total_df[ , colnames(total_df) %in% c('all_variants', 'mut_rate', 'gc_content')]
# X = as.matrix(X)
# ================

# ===== BETA =====
# (Scale and center features)
# total_df$mut_rate = scale_and_center(total_df$mut_rate)
# total_df$all_variants = scale_and_center(total_df$all_variants)
# total_df$gc_content = scale_and_center(total_df$gc_content)
# total_df$y = scale_and_center(total_df$y)
# ================





# ========== BETA Regression Models ==========
#regr_model = glm(formula = y ~ all_variants + mut_rate + gc_content + cpg, data=total_df) 

#total_df$weighted_bins_sum = total_df[, 'bin_1'] + (1/5) * total_df[, 'bin_2'] + (1/10) * total_df[, 'bin_3'] + 
#					(1/50) * total_df[, 'bin_4'] + (1/200) * total_df[, 'bin_5'] + 
#					(1/1000) * total_df[, 'bin_6']
#regr_model = glm(formula = y ~ weighted_bins_sum + mut_rate, data=total_df) 
#regr_model = lmridge(y ~ weighted_bins_sum + mut_rate, data=total_df, K=-0.2, scaling=c("sc")) # ridge regression 
#regr_model = lmridge(y ~ all_variants + mut_rate + all_ac + all_af + gc_content, data=total_df, K=-0.2, scaling=c("sc")) # ridge regression 


# > Multivariate Regressin
#regr_model = lm(cbind(bin_1, bin_2, bin_3, bin_4, bin_5, bin_6) ~ all_variants + mut_rate + common_div_by_all_var_ratio + gc_content, data=total_df)
#regr_model = lm(cbind(bin_1, bin_2, bin_3, bin_4, bin_5, bin_6) ~ all_variants , data=total_df)


# > With cross  validation 
#regr_model = cv.glmnet(X, Y, alpha=0, standardize=TRUE, type.measure='auc') # ridge regression (alpha=0) - can't get residuals
#regr_model = cv.glmnet(X, Y, alpha=1, standardize=TRUE, keep=T, type.measure = "mse") # lasso (alpha=1) - can't get residuals

# ======== BEST ======== :
# regr_model = lmridge(y ~ all_variants + mut_rate + cpg + gc_content, data=total_df, K=0.5, scaling=c("sc")) # ridge regression 
# regr_model = lmridge(y ~ all_variants + cpg, data=total_df, K=-0.2, scaling=c("sc")) # ridge regression 


# --- SUPPORT VECTOR REGRESSION ---
# - Tuning:
#tuneResult = tune(svm, y ~ all_variants + mut_rate + gc_content, data=total_df, ranges=list(epsilon = seq(0,1,0.1), cost = 2^(2:9)))
#print(tuneResult)
#pdf(paste(out_dir, '/svm_tuning_regr_model.pdf', sep=''))
#plot(tuneResult)
#dev.off()

#regr_model = svm(y ~ all_variants + mut_rate + gc_content, data=total_df)
#predictedY = predict(regr_model, total_df)
#error = total_df$y - predictedY
# ----------------------------------


# =====  BETA: Getting residuals =====
# > Get (studentised?) residuals from lmridge
# stud_res = residuals(regr_model) # <<- with lmridge


# > Multivariate:
#stud_res_sum = stud_res[, 'bin_1'] + (1/5) * stud_res[, 'bin_2'] + (1/10) * stud_res[, 'bin_3'] + 
#					(1/50) * stud_res[, 'bin_4'] + (1/200) * stud_res[, 'bin_5'] + 
#					(1/1000) * stud_res[, 'bin_6']
#total_df = cbind(total_df, 'stud_res' = as.numeric(stud_res_sum), row.names=rownames(stud_res))




annotate_windows_by_tolerance <- function(final_df){

	final_df['annot'] = 'other'
	print(head(final_df))
	inferred_chr = rownames(final_df)[1]
	inferred_chr = gsub('chr', '', inferred_chr)
	inferred_chr = gsub('_.*', '', inferred_chr)
	print(inferred_chr)

	bottom_prob = .01
	up_prob = 1- bottom_prob

	extremes = quantile(final_df$stud_res, c(bottom_prob, up_prob), na.rm=T)
	print(extremes)

	#print(nrow(final_df[ final_df$stud_res <= extremes[1], ]))
	#print(nrow(final_df[ final_df$stud_res > extremes[1] & final_df$stud_res < extremes[2], ]))
	#print(nrow(final_df[ final_df$stud_res >= extremes[2], ]))

	final_df[ final_df$stud_res <= extremes[1], "annot" ] = 'intolerant'
	final_df[ final_df$stud_res >= extremes[2], "annot" ] = 'tolerant'

	final_df = final_df[ order(final_df$stud_res), ]
	#print(head(final_df, 10))
	#print(tail(final_df, 10))

	toler_colors = c('#de2d26', '#3182bd', '#bdbdbd')
	names(toler_colors) = c('intolerant', 'tolerant', 'other')
	final_df$annot[ final_df$annot == 'intolerant' ] = '#de2d26'
	final_df$annot[ final_df$annot == 'tolerant' ] = '#3182bd'
	final_df$annot[ final_df$annot == 'other' ] = '#bdbdbd'
	#print(head(final_df))

	if( chr_type == 'autosomal'){
		png(paste(out_dir, '/scatter/autosomal-regression_scatterplot.png', sep=''), height = 10, width = 10, units = 'in', res = 600)
		#plot(final_df$X, final_df$y, xlab='All Variants', ylab='Common variants', pch=16, cex=0.8, col=final_df$annot)
		plot(final_df$mut_rate, final_df$y, xlab='All Variants', ylab='Common variants', pch=16, cex=0.8, col=final_df$annot)
		abline(regr_model, col='#31a354', lwd=2)
		dev.off()
	}
}
#annotate_windows_by_tolerance(total_df)


# Unfold studentised residuals for each chromosome
unfold_studres_from_each_chr <- function(final_df){

	for(chr in all_chrs){
		print(paste('chr', chr, sep=':'))

		total_entries = as.numeric(total_entries_per_chr[chr])
		print(total_entries)

		chr_idx = paste('chr', chr, '_', sep='')

	
		print('---------')
		df = final_df[ grepl(chr_idx, rownames(final_df)), ]
		print(head(df))

		df$row_names = rownames(df)
		df$row_names = gsub(chr_idx, '', df$row_names)
		print(head(df))	

		df = df[ order(as.numeric(df$row_names)), ]
		print(head(df, 10))
		print(tail(df, 10))

		stud_res = as.numeric(df$stud_res)
		names(stud_res) = df$row_names
		print(head(stud_res))



		# ========== Fill-in windows that have no associated values ==========
		merged_scores = stud_res
		last_idx = as.numeric(df[ nrow(df), 'row_names']) 

		if( last_idx != nrow(df)-1 ){ 
			print("Dealing with chromosomes that have windows with no values...")
			print(paste('[last_idx != nrow(df)] - now filling missing values at chr:', chr)) 
			print(paste('last_idx:', last_idx, 'nrow(df)-1', nrow(df)-1))
			print(paste('total_entries:', total_entries))
			print(paste('typeof - total_entries:', str(total_entries)))
			
			zero_win_indexes = setdiff(seq(0, total_entries - 1), df$row_names)
			#print(zero_win_indexes)

			tmp_zeros_vec = rep('NaN', length(zero_win_indexes))
			names(tmp_zeros_vec) = zero_win_indexes
			#print(tmp_zeros_vec)

			merged_scores = c(merged_scores, tmp_zeros_vec)
			merged_scores = merged_scores[ order(as.numeric(names(merged_scores))) ]
			#print( tail( head(merged_scores, 17176), 5 ))
			print(head(merged_scores, 10))
			print(tail(merged_scores, 10))
		}
		# ====================================================================

		#cur_regr_model = glm(formula = y ~ all_variants + mut_rate + gc_content, data=df) # regression for current chromosome only

		dir.create(file.path(out_dir, '/scatter'), showWarnings = FALSE)
		pdf(paste(out_dir, '/scatter/regression_scatterplot_chr', chr, '.pdf', sep=''), height = 10, width = 10)

		plot(df$all_variants, df$y, xlab='All Variants', ylab='Common variants', pch=16, cex=1)
		#abline(cur_regr_model, col='#31a354', lwd=2)
		dev.off()

		cat(merged_scores, file=paste(gwrvis_dir, '/studres.chr', chr, '.txt', sep=''), fill=F)
	}
}
unfold_studres_from_each_chr(total_df)
