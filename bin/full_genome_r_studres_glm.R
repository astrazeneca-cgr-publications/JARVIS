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
chr_type = args[5]


rvis_dir = paste(out_dir, 'rvis_scores', sep='/')
input_path = paste(out_dir, "/tmp/", sep='')
files = NULL
if( chr_type == 'autosomal'){
	files = list.files(path = input_path, pattern = "^Xy.")
} else if( chr_type == 'sex'){
	files = list.files(path = input_path, pattern = "^Xy.chrX")
}
print(files)

total_entries_list = list()
total_df = data.frame()
all_chrs = vector()

for(file in files){

	chr = file
	chr = gsub("Xy.chr", '', chr)
	chr = gsub(".txt", '', chr)
	print(file)
	# Xy.chr10.txt
	print(chr)

	# process X chromosome separately if input chr_type == autosomal
	if(chr == 'X' & chr_type == 'autosomal'){
		next
	}
	all_chrs = c(all_chrs, chr)
	
	
	df = read.table(paste(input_path, file, sep='/'), sep=',', header=T)
	rownames(df) = paste('chr', chr, '_', df$idx, sep='')
	df$idx = NULL
	print(head(df, 20))
	total_entries_list[chr] = nrow(df)

	if(nrow(total_df) > 0){
		total_df = rbind(total_df, df)
	} else{
		total_df = df
	}
	
	print(nrow(total_df))
}


# BETA - [OPTION]: exclude regions with '0' variants from the regression fitting analysis (so that it doesn't skew the fit)
#print(nrow(total_df))
#total_df = total_df[ total_df$X != 0, ]

print(head(total_df))
print(dim(total_df))




if( all_variants_upper_thres == -1){
	all_variants_upper_thres = win_len 
}
print(paste('all_variants_upper_thres:', all_variants_upper_thres))


print(nrow(total_df))
## Remove windows that had no variant data in the original VCF file
placeholder_val = -1
#total_df = total_df[ total_df$X != placeholder_val, ]
total_df = total_df[ total_df$mut_rate != placeholder_val, ]
print(nrow(total_df))






is_outlier <- function(vec, z_thres=3.5){          
	#vec = vec[ !is.nan(vec) ]         
	median = median(na.omit(vec))        
	#print(paste('median:', median))          
	#print(head(vec))         

	diff = sqrt((vec - median)**2)         
	#print(head(diff))          
	med_abs_deviation = median(na.omit(diff))         
	modified_z_score = 0.6745 * diff / med_abs_deviation         

	#print(head(modified_z_score))          
	return(modified_z_score > z_thres) 
} 




# [Conditionally] exclude outliers before fitting regression model
if(as.logical(filter_outliers_before_regression)){
	print('filtering before regression...')
	total_df = total_df[ total_df$all_variants <= all_variants_upper_thres, ]
}

# ===== BETA ====== (normalise mut_rate and gc_content)
total_df$mut_rate = total_df$mut_rate / ( 7 * 3000) # ( kmer * win_len )
#total_df$gc_content = total_df$gc_content / 3000 # ( win_len )

Ne = 10000
total_df$theta = 4 * Ne * total_df$mut_rate
# =============================================================



## ===== BETA ===== (scale total number of variants by mutation rate)
# final_df$all_variants = final_df$all_variants / final_df$mut_rate
# =============================================================




# == BETA == Ridge regression (L2 regularization)
Y = total_df$y
print(str(Y))
X = total_df[ , colnames(total_df) %in% c('all_variants', 'mut_rate', 'gc_content')]
print(str(X))
X = as.matrix(X)
print(str(X))
# ==========


convert_pval_to_phred <- function(pval){
	phred = -10 * log10(as.numeric(pval))
	
	return(phred)
}


## ===== NEW =====
# Convert p-values to phred scores
#total_df$mut_rate = sapply(total_df$mut_rate, function(x) convert_pval_to_phred(x))


# Scale and center features
scale_and_center <- function(x){
	x = as.numeric(x)
	return( (x - mean(x)) / sd(x) ) 
}
#total_df$mut_rate = scale_and_center(total_df$mut_rate)
#total_df$all_variants = scale_and_center(total_df$all_variants)
#total_df$gc_content = scale_and_center(total_df$gc_content)
#total_df$y = scale_and_center(total_df$y)

## {{ ... or range standardise (0 to 1) - NOT giving good results
range01 <- function(x){
	return( (x-min(x)) / (max(x)-min(x)) )
}
##total_df$mut_rate = range01(total_df$mut_rate)
##total_df$all_variants = range01(total_df$all_variants)
##total_df$gc_content = range01(total_df$gc_content)
##total_df$y = range01(total_df$y)
## ... }}



#total_df$weighted_bins_sum = total_df[, 'bin_1'] + (1/5) * total_df[, 'bin_2'] + (1/10) * total_df[, 'bin_3'] + 
#					(1/50) * total_df[, 'bin_4'] + (1/200) * total_df[, 'bin_5'] + 
#					(1/1000) * total_df[, 'bin_6']
total_df$weighted_bins_sum = total_df[, 'bin_1'] + (1/5) * total_df[, 'bin_2'] + (1/10) * total_df[, 'bin_3'] + 
					(1/50) * total_df[, 'bin_4'] + (1/200) * total_df[, 'bin_5'] + 
					(1/1000) * total_df[, 'bin_6']

# beta:
#bins_distance = total_df$theta - total_df$weighted_bins_sum



## ===================>>>>>>>> Original <<<<<<<<<=====================
#regr_model = glm(formula = y ~ all_variants, data=total_df) 
#regr_model = glm(formula = y ~ all_variants + mut_rate + gc_content + cpg, data=total_df) 

# GOOD:
#regr_model = glm(formula = y ~ weighted_bins_sum + mut_rate, data=total_df) 
#regr_model = lmridge(y ~ weighted_bins_sum + mut_rate, data=total_df, K=-0.2, scaling=c("sc")) # ridge regression 



print(head(total_df))

# cols: all_ac       all_af all_variants  y    mut_rate gc_content    theta
# BETA
#regr_model = glm(formula = y ~ all_variants + theta, data=total_df) 
#regr_model = glm(formula = all_variants ~ theta + y, data=total_df) 

#regr_model = lmridge(y ~ all_variants + mut_rate + all_ac + all_af + gc_content, data=total_df, K=-0.2, scaling=c("sc")) # ridge regression 


# BETA: Multivariate Regressin
#regr_model = lm(cbind(bin_1, bin_2, bin_3, bin_4, bin_5, bin_6) ~ all_variants + mut_rate + common_div_by_all_var_ratio + gc_content, data=total_df)
#regr_model = lm(cbind(bin_1, bin_2, bin_3, bin_4, bin_5, bin_6) ~ all_variants , data=total_df)



# LINEAR REGRESSION
#regr_model = cv.glmnet(X, Y, alpha=0, standardize=TRUE, type.measure='auc') # ridge regression (alpha=0) - can't get residuals
#regr_model = cv.glmnet(X, Y, alpha=1, standardize=TRUE, keep=T, type.measure = "mse") # lasso (alpha=1) - can't get residuals

#regr_model = lmridge(y ~ all_variants + mut_rate + gc_content, data=total_df, K=-0.2, scaling=c("sc")) # ridge regression 
# ======== BEST ======== :
#regr_model = lmridge(y ~ all_variants + mut_rate + cpg + gc_content, data=total_df, K=0.5, scaling=c("sc")) # ridge regression 
regr_model = lmridge(y ~ all_variants + cpg, data=total_df, K=-0.2, scaling=c("sc")) # ridge regression 


#regr_model = lm(y ~ all_variants + mut_rate + gc_content, data=total_df)


## --- SUPPORT VECTOR REGRESSION ---
# tuning:
#tuneResult = tune(svm, y ~ all_variants + mut_rate + gc_content, data=total_df, ranges=list(epsilon = seq(0,1,0.1), cost = 2^(2:9)))
#print(tuneResult)
#pdf(paste(out_dir, '/svm_tuning_regr_model.pdf', sep=''))
#plot(tuneResult)
#dev.off()
#stop()
# ----------------------------------

#regr_model = svm(y ~ all_variants + mut_rate + gc_content, data=total_df)
#predictedY = predict(regr_model, total_df)
#error = total_df$y - predictedY



# Get RMSE (Root Mean Square Error)
error = residuals(regr_model)
RMSE = sqrt(mean(error^2))
print("==============================================")
print(paste('RMSE:', RMSE))
print("==============================================")
#stop()


# >>>>>>> ROC - Only for classification <<<<<<<
#prob <- predict(regr_model, type = c("response"))
#plot(roc(total_df$y, prob), print.auc = TRUE)


#pdf(paste(out_dir, '/Summary_regr_model.pdf', sep=''))
#plot(regr_model)
#dev.off()

## >>>>>>> DBG:
#resid = resid(regr_model) # lm - residuals
#print('--------')
#print(resid)
#print('--------')
#stop()
### <<<<<<<

# Lasso - Results 
if(0){
	cv.lasso = regr_model
	pdf(paste(out_dir, '/lasso_regr_model.pdf', sep=''))
	plot(cv.lasso)
	plot(cv.lasso$glmnet.fit, xvar="lambda", label=TRUE)
	plotres(cv.lasso)
	dev.off()
	print(cv.lasso$lambda.min) 
	print(cv.lasso$lambda.1se)
	print(coef(cv.lasso, s=cv.lasso$lambda.min))
}

#mod = regr_model
## Scaled Coefficients 
#print(mod$coef)


## Re-Scaled Coefficients 
#print(coef(mod)) 

## ridge predicted values 
#print(head(predict(mod)))

## ridge residuals 
#print(head(residuals(mod)))

##ridge and VIF trace -- too SLOW
#pdf(paste(out_dir, '/lmridge_regr_model.pdf', sep=''))
#plot(mod)
#dev.off()

## ridge VIF values
#print(head(vif(mod)))

## ridge Var-Cov matrix 
#print(vcov(mod))

## ridge biasing parameter by researchers 
#print(kest(mod)) 

## ridge fitted values 
#print(head(fitted(mod)))

## ridge statistics 1 
#print(rstats1(mod)) 
## ridge statistics 2 
#print(rstats2(mod))

#print(summary(regr_model))

## ******* Original: *******
#stud_res = studres(regr_model)

# ====================  BETA ================================
## ** Alternatives instead of stud-res - Beta **
stud_res = residuals(regr_model) # <<- with lmridge


# BETA - multivariate:
#stud_res_sum = stud_res[, 'bin_1'] + (1/5) * stud_res[, 'bin_2'] + (1/10) * stud_res[, 'bin_3'] + 
#					(1/50) * stud_res[, 'bin_4'] + (1/200) * stud_res[, 'bin_5'] + 
#					(1/1000) * stud_res[, 'bin_6']
#total_df = cbind(total_df, 'stud_res' = as.numeric(stud_res_sum), row.names=rownames(stud_res))


total_df = cbind(total_df, 'stud_res' = as.numeric(stud_res))
print(head(total_df))


# BETA - multivariate:
#stud_res = stud_res_sum



## -> TO-DO: check that multiplying factor is not zero before multiplying
## GOOD:
total_neg_df = subset(total_df, stud_res < 0)

## State-of-the-art:
#total_neg_df$stud_res = ( (total_neg_df$stud_res) * (total_neg_df$bin_1 + 1) * (total_neg_df$bin_2 + 1) / (total_neg_df$bin_6 + 1) ) / ((total_neg_df$y + 1) / (total_neg_df$all_variants + 1)) 

## == BETA:
#total_neg_df$stud_res = -log(abs(total_neg_df$stud_res))

# -> total_neg_df$stud_res = (total_neg_df$stud_res) * total_neg_df$bin_1 * total_neg_df$bin_2 / total_neg_df$bin_6 # * total_neg_df$bin_4 * total_neg_df$bin_5 * total_neg_df$bin_6

#total_neg_df$stud_res = total_neg_df$stud_res * (total_neg_df$bin_1 + 1) * (total_neg_df$bin_2 + 1) / (total_neg_df$bin_6 + 1)
#total_neg_df$stud_res = total_neg_df$stud_res * total_neg_df$bin_1 * total_neg_df$bin_2 / (total_neg_df$bin_6 + 1) * total_neg_df$theta <-- maybe?
print(head(total_neg_df))

print('========================')

total_pos_df = subset(total_df, stud_res >= 0)
## State-of-the-art:
#total_pos_df$stud_res = ( (total_pos_df$stud_res) / (total_pos_df$bin_1 + 1) / (total_pos_df$bin_2 + 1) * (total_pos_df$bin_6) ) / ((total_pos_df$y + 1) / (total_pos_df$all_variants + 1))

## == BETA:
#total_pos_df$stud_res = log(total_pos_df$stud_res + 1)

#total_pos_df$stud_res = 1 / total_pos_df$stud_res

#total_pos_df$stud_res = ( (total_pos_df$stud_res) / (total_pos_df$bin_1+1) / (total_pos_df$bin_2 + 1) * (total_pos_df$bin_6) ) * ((total_concat_df$y + 1) / (total_concat_df$all_variants + 1))

# -> total_pos_df$stud_res = total_pos_df$stud_res * (total_pos_df$bin_6 + 1) # total_pos_df$bin_1 * total_pos_df$bin_2 / total_pos_df$bin_5 # * total_pos_df$bin_4 * total_pos_df$bin_5 * total_pos_df$bin_6

#total_pos_df$stud_res = total_pos_df$stud_res / (total_pos_df$bin_1 + 1) / (total_pos_df$bin_2 + 1) * (total_pos_df$bin_6 + 1)
#total_pos_df$stud_res = total_pos_df$stud_res * total_pos_df$bin_1 * total_pos_df$bin_2 / (total_pos_df$bin_6 + 1) * total_pos_df$theta <-- maybe?
print(head(total_pos_df))



total_concat_df = rbind(total_pos_df, total_neg_df)
total_concat_df = total_concat_df[match(rownames(total_df), rownames(total_concat_df)),]



# common_div_by_all_var_ratio	y	ac	af	bin_1	bin_2	bin_3	bin_4	bin_5	bin_6	all_variants	mut_rate	gc_content	theta weighted_bins_sum	stud_res
# DBG:
###total_concat_df$stud_res = total_concat_df$stud_res + (total_concat_df$bin_6 - total_concat_df$bin_1) / total_concat_df$stud_res
#total_concat_df$stud_res = total_concat_df$stud_res * total_concat_df$

# BEST:
total_concat_df$stud_res = ( total_concat_df$stud_res * (total_concat_df$bin_1 + 1) * (total_concat_df$bin_2 + 1) / (total_concat_df$bin_6 + 1) ) / ((total_concat_df$y + 1) / (total_concat_df$all_variants + 1)) 


stud_res = total_concat_df$stud_res
print(max(stud_res))
print(min(stud_res))


# ==============================
#stud_res = stud_res * total_df$bin_1 * total_df$bin_2 / (total_df$bin_6 + 1) 
# or (not so good):
#stud_res = stud_res / total_df$y
# ==============================

# === BETA ===
#stud_res = stud_res[ , colnames(stud_res) == 'bin_1'] #+ 
	           #stud_res[ , colnames(stud_res) == 'bin_2'] + 
		   #stud_res[ , colnames(stud_res) == 'bin_3'] +
		   #stud_res[ , colnames(stud_res) == 'bin_4'] + 
		   #stud_res[ , colnames(stud_res) == 'bin_5'] + 
		   #stud_res[ , colnames(stud_res) == 'bin_6'] 
	           #(1/5) * stud_res[ , colnames(stud_res) == 'bin_2'] + 
		   #(1/10) * stud_res[ , colnames(stud_res) == 'bin_3'] +
		   #(1/50) * stud_res[ , colnames(stud_res) == 'bin_4'] + 
		   #(1/200) * stud_res[ , colnames(stud_res) == 'bin_5'] + 
		   #(1/1000) * stud_res[ , colnames(stud_res) == 'bin_6'] 

## = BETA
#stud_res = total_df$theta - total_df$weighted_bins_sum			
#print(head(stud_res))

# === GOOD BETAs ===
#stud_res = total_df$y - total_df$all_variants
# ----> ! 
#stud_res = total_df$y - total_df$weighted_bins_sum
#print(head(stud_res))



#N=15496
#stud_res = N / (total_df$all_ac) 
#print(head(total_df))
#stud_res = 2 * N / (total_df$all_af) 

#stud_res = predict(mod) # <<- with lmridge

#stud_res = fitted(mod)

#stud_res = resid(regr_model) # lm - residuals

#stud_res = error # svm
## ***************************************************************



final_df = total_concat_df
final_df$row_names = rownames(final_df)
#final_df = merge(total_df, stud_res_df, by=0, all=TRUE)
#colnames(final_df)[1] = 'row_names'

## =========================== BETA ===========================
#final_df$stud_res = final_df$stud_res * total_df$all_ac

#final_df$stud_res = sign(final_df$stud_res) * (final_df$stud_res^2) -- May be working well
print(head(final_df))

# =============================================================


if(0){
	final_df['annot'] = 'other'
	print(head(final_df))

	bottom_prob = .01
	up_prob = 1- bottom_prob

	extremes = quantile(final_df$stud_res, c(bottom_prob, up_prob), na.rm=T)
	print(extremes)
	print(nrow(final_df[ final_df$stud_res <= extremes[1], ]))
	print(nrow(final_df[ final_df$stud_res > extremes[1] & final_df$stud_res < extremes[2], ]))
	print(nrow(final_df[ final_df$stud_res >= extremes[2], ]))



	final_df[ final_df$stud_res <= extremes[1], "annot" ] = 'intolerant'
	final_df[ final_df$stud_res >= extremes[2], "annot" ] = 'tolerant'

	final_df = final_df[ order(final_df$stud_res), ]
	print(head(final_df, 10))
	print(tail(final_df, 10))


	toler_colors = c('#de2d26', '#3182bd', '#bdbdbd')
	names(toler_colors) = c('intolerant', 'tolerant', 'other')
	final_df$annot[ final_df$annot == 'intolerant' ] = '#de2d26'
	final_df$annot[ final_df$annot == 'tolerant' ] = '#3182bd'
	final_df$annot[ final_df$annot == 'other' ] = '#bdbdbd'
	print(head(final_df))


	if( chr_type == 'autosomal'){
		png(paste(out_dir, '/scatter/whole_genome-regression_scatterplot.png', sep=''), height = 10, width = 10, units = 'in', res = 600)
		#plot(final_df$X, final_df$y, xlab='All Variants', ylab='Common variants', pch=16, cex=0.8, col=final_df$annot)
		plot(final_df$mut_rate, final_df$y, xlab='All Variants', ylab='Common variants', pch=16, cex=0.8, col=final_df$annot)
		abline(regr_model, col='#31a354', lwd=2)
		dev.off()
	}
}



# Unfold studentised residuals for each chromosome
for(chr in all_chrs){
	print(paste('chr', chr, sep=':'))

	total_entries = as.numeric(total_entries_list[chr])
	print(total_entries)

	chr_idx = paste('chr', chr, '_', sep='')

	print(nrow(final_df))
	df = final_df[ grep(chr_idx, final_df$row_names), ]
	print(nrow(df))
	print(head(df))
	print(tail(df))

	df$row_names = gsub(chr_idx, '', df$row_names)
	#print(head(df))	

	df = df[ order(as.numeric(df$row_names)), ]
	# colnames(df)[ncol(df)] = 'stud_res'
	print(head(df, 20))
	print(tail(df, 10))


	stud_res = as.numeric(df$stud_res)
	names(stud_res) = df$row_names
	#print(head(stud_res))
	#print( tail( head(stud_res, 17176), 5 ))

	merged_scores = stud_res

	# <----- Potentially Redundant ------>
	last_idx = as.numeric(df[ nrow(df), 'row_names']) 
	if( last_idx != nrow(df)-1 ){ # deal with chromosomes that have windows with no values
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
	# <---------------------->
	#print(head(df))

	#cur_regr_model = glm(formula = y ~ all_variants + mut_rate + gc_content, data=df) # regression for current chromosome only

	dir.create(file.path(out_dir, '/scatter'), showWarnings = FALSE)
	#png(paste(out_dir, '/scatter/regression_scatterplot_chr', chr, '.png', sep=''), height = 10, width = 10, units = 'in', res = 600)
	pdf(paste(out_dir, '/scatter/regression_scatterplot_chr', chr, '.pdf', sep=''), height = 10, width = 10)

	plot(df$all_variants, df$y, xlab='All Variants', ylab='Common variants', pch=16, cex=1)
	#abline(cur_regr_model, col='#31a354', lwd=2)

	
	plot(df$gc_content, df$y, xlab='GC content', ylab='Common variants', pch=16, cex=1)
	#abline(cur_regr_model, col='#31a354', lwd=2)

	plot(df$mut_rate, df$y, xlab='Mutability rate', ylab='Common variants', pch=16, cex=1)
	#abline(cur_regr_model, col='#31a354', lwd=2)
	dev.off()

	cat(merged_scores, file=paste(rvis_dir, '/studres.chr', chr, '.txt', sep=''), fill=F)


	
	df_to_save = data.frame('row_names' = names(merged_scores), 'stud_res' = merged_scores)
	#df_to_save[ df_to_save$stud_res == 'NaN', 'stud_res'] = ''
	#print(head(df_to_save))
	write.table(df_to_save, paste(rvis_dir, '/rvis_scores_chr', chr, '.csv', sep=''), sep=',', row.names=F, col.names=F, quote=F)
}
