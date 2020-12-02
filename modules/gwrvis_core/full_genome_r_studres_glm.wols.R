library(glm2)
library(MASS)
library(glmnet)
library(lmridge)
library(plotmo)
library(pROC)




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



#args = commandArgs(trailingOnly=T)
args <- "out/topmed-heteroskedasticity-winlen_3000.MAF_0.001/"
out_dir = args[1]



# --- Static parameters ----
all_variants_upper_thres = -1	# as.numeric(args[2])
win_len = 3000		# as.numeric(args[3])
chr_type = 'autosomal'  # args[4]





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
all_chrs = seq(1,22)




total_df_file = paste(input_path, 'total_df.Xy.tsv', sep='/')

# Read files with common/all variants data per window / per chromosome
if( file.exists(total_df_file) ){

	print(">> Reading pre-compiled total_df.Xy.tsv file ...")
	total_df = read.table(total_df_file, header=T)
	#print(head(total_df))

	for(chr in all_chrs){
		tmp_df = total_df[ grepl(paste('chr', chr, '_', sep=''), rownames(total_df)), ]
		total_entries_per_chr[chr] = nrow(tmp_df)
	}

} else {
	all_chrs = vector()

	for(file in files){

		chr = file
		chr = gsub("Xy.chr", '', chr)
		chr = gsub(".txt", '', chr)

		# process X chromosome separately if input chr_type == autosomal
		if(chr == 'X' & chr_type == 'autosomal'){ next }

		print(paste('>chr', chr))
		all_chrs = c(all_chrs, chr)
		
		
		df = read.table(paste(input_path, file, sep='/'), sep=',', header=T)
		rownames(df) = paste('chr', chr, '_', df$idx, sep='')
		df$idx = NULL
		#print(head(df, 20))
		#print(tail(df, 20))
		total_entries_per_chr[chr] = nrow(df)

		if(nrow(total_df) > 0){
			total_df = rbind(total_df, df)
		} else{
			total_df = df
		}
		
		print(nrow(total_df))
	}

	write.table(total_df, file=total_df_file, quote=F)
	print(paste("Saved total_df to file:", total_df_file))

}



# Remove windows that had no variant data in the original VCF file
placeholder_val = -1
print(nrow(total_df))
total_df = total_df[ total_df$mut_rate != placeholder_val, ]
print(nrow(total_df))





# [Not used in production]: excludes windows with number of variants over a threshold
if(all_variants_upper_thres != -1){
	print('filtering before regression...')
	total_df = total_df[ total_df$all_variants <= all_variants_upper_thres, ]
}





## ============================= Main Analysis =============================
# y: common variants
cat("\n Running Linear Regression...\n")
# regr_model = glm(formula = y ~ all_variants, data=total_df) 
# regr_model = lm(formula = y ~ all_variants, data=total_df)
regr_model = lm(formula = y ~ all_variants, data=total_df, weights=1/all_variants) 



# Get RMSE (Root Mean Square Error)
error = residuals(regr_model)
RMSE = sqrt(mean(error^2))
print(paste('RMSE:', RMSE))

# Get studentised residuals
stud_res = studres(regr_model)




total_df = cbind(total_df, 'stud_res' = as.numeric(stud_res))

pdf('stud_res_weighted.pdf')
plot(total_df$all_variants, total_df$stud_res,pch='.')
dev.off()


# 
# varfunc.ols <- lm(log(stud_res^2) ~ log(all_variants), data=total_df)
# total_df$varfunc <- exp(varfunc.ols$fitted.values)
# 
# gls <-  lm(y ~ all_variants, weights = 1/sqrt(varfunc), data = total_df)
# 
# 
# total_df$stud_res <- studres(gls)
# 
# 




# ==========================================================================






# ============================ Post-processing and Plotting ============================
annotate_windows_by_tolerance <- function(final_df){

	final_df['annot'] = 'other'
	#print(head(final_df))

	bottom_prob = .01
	up_prob = 1- bottom_prob

	extremes = quantile(final_df$stud_res, c(bottom_prob, up_prob), na.rm=T)
	#print(extremes)

	final_df[ final_df$stud_res <= extremes[1], "annot" ] = 'intolerant'
	final_df[ final_df$stud_res >= extremes[2], "annot" ] = 'tolerant'
	final_df = final_df[ order(final_df$stud_res), ]


	toler_colors = c('#de2d26', '#3182bd', '#bdbdbd')
	names(toler_colors) = c('intolerant', 'tolerant', 'other')
	final_df$annot[ final_df$annot == 'intolerant' ] = '#de2d26'
	final_df$annot[ final_df$annot == 'tolerant' ] = '#3182bd'
	final_df$annot[ final_df$annot == 'other' ] = '#bdbdbd'
	#print(head(final_df))

	if( chr_type == 'autosomal'){
		png(paste(out_dir, '/scatter/autosomal-regression_scatterplot.png', sep=''), height = 10, width = 10, units = 'in', res = 600)
		#plot(final_df$X, final_df$y, xlab='All Variants', ylab='Common variants', pch=16, cex=0.8, col=final_df$annot)
		plot(final_df$all_variants, final_df$y, xlab='All Variants', ylab='Common variants', pch=16, cex=0.8, col=final_df$annot)
		abline(regr_model, col='#31a354', lwd=2)
		dev.off()
	}
}
cat("\n Plotting scatterplot with genome-wide regression fit ...\n")
annotate_windows_by_tolerance(total_df)




# Unfold studentised residuals for each chromosome
unfold_studres_from_each_chr <- function(final_df){

	for(chr in all_chrs){
		print(paste('chr:', chr))

		total_entries = as.numeric(total_entries_per_chr[chr])
		#print(total_entries)

		chr_idx = paste('chr', chr, '_', sep='')

	
		df = final_df[ grepl(chr_idx, rownames(final_df)), ]
		#print(head(df))

		df$row_names = rownames(df)
		df$row_names = gsub(chr_idx, '', df$row_names)
		#print(head(df))	

		df = df[ order(as.numeric(df$row_names)), ]
		#print(head(df, 10))
		#print(tail(df, 10))

		stud_res = as.numeric(df$stud_res)
		names(stud_res) = df$row_names
		#print(head(stud_res))



		# ========== Fill-in windows that have no associated values ==========
		merged_scores = stud_res
		last_idx = as.numeric(df[ nrow(df), 'row_names']) 

		if( last_idx != nrow(df)-1 ){ 
			#print("Dealing with chromosomes that have windows with no values...")
			#print(paste('[last_idx != nrow(df)] - now filling missing values at chr:', chr)) 
			#print(paste('last_idx:', last_idx, 'nrow(df)-1', nrow(df)-1))
			#print(paste('total_entries:', total_entries))
			
			zero_win_indexes = setdiff(seq(0, total_entries - 1), df$row_names)
			#print(zero_win_indexes)

			tmp_zeros_vec = rep('NaN', length(zero_win_indexes))
			names(tmp_zeros_vec) = zero_win_indexes
			#print(tmp_zeros_vec)

			merged_scores = c(merged_scores, tmp_zeros_vec)
			merged_scores = merged_scores[ order(as.numeric(names(merged_scores))) ]
			#print( tail( head(merged_scores, 17176), 5 ))
			#print(head(merged_scores, 10))
			#print(tail(merged_scores, 10))
		}
		# ====================================================================


		# > Plot regression fit per chromosome
		#dir.create(file.path(out_dir, '/scatter'), showWarnings = FALSE)
		#pdf(paste(out_dir, '/scatter/regression_scatterplot_chr', chr, '.pdf', sep=''), height = 10, width = 10)
		#plot(df$all_variants, df$y, xlab='All Variants', ylab='Common variants', pch=16, cex=1)
		#dev.off()

		cat(merged_scores, file=paste(gwrvis_dir, '/studres.chr', chr, '.txt', sep=''), fill=F)
	}
}
cat("\n Storing stud. residuals per chromosome ...\n")
unfold_studres_from_each_chr(total_df)

cat("\nDONE\n")
