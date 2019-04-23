import matplotlib 
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from datetime import datetime
import sys
import os
from subprocess import call
from custom_utils import create_out_dir, get_run_parameters
from random import shuffle
import subprocess
from pathlib import Path


startTime = datetime.now()

args = sys.argv
if len(args) != 3:
	print("[Error]: insufficient input arguments. Expected command call:\n> python full_genome_collapsed_chr_parser.py [chr] [config_file]")
	sys.exit()
 
chr = args[1]
config_file = args[2] #'config.log'

# read run parameters from config file and store into a dictionary
run_params = get_run_parameters(config_file)
print(run_params)


dataset = run_params['dataset']		# e.g. 'gnomad'
win_len = run_params['win_len']		# e.g. 250 (window length in nt)
all_variants_upper_thres = run_params['all_variants_upper_thres']	# e.g. 200 (filter out windows with more than 200 variants before fitting regression)
MAF_thres = run_params['MAF_thres']	# e.g. 0.0001 (Minor Allele Frequency)
filter_outliers_before_regression = run_params['filter_outliers_before_regression'] 	# e.g. True
generate_intermediate_plots = run_params['generate_intermediate_plots']		# e.g. False
variants_table_dir = run_params['variants_table_dir']


rscript_path = 'Rscript'
cluster_rscript_path = '/users/kclc950/packages/R-devel/bin/Rscript'
if os.path.exists(cluster_rscript_path):
	rscript_path = cluster_rscript_path

human_ref_genome_2bit = '../hg19/homo_sapiens_GRCh37_FASTA/hsa37.2bit'



# =========================================

data_dir = '../' + dataset + '/' + variants_table_dir
print('data_dir: *' + data_dir + '*')


simplified_vcf_file = data_dir + '/chr' + chr + '_' + dataset + '_table.txt.collapsed'


# create out_MAF{threshold} dir to store output results (plots and gwRVIS csv files)
out_dir = create_out_dir(config_file)
print(out_dir)


var_ratios_dir = out_dir + '/var_ratios'
if not os.path.exists(var_ratios_dir):     
	os.makedirs(var_ratios_dir, exist_ok=True)

rvis_dir = out_dir + '/rvis_scores'
if not os.path.exists(rvis_dir):     
	os.makedirs(rvis_dir, exist_ok=True)

plots_dir = out_dir + '/plots_per_chrom'
if not os.path.exists(plots_dir):     
	os.makedirs(plots_dir, exist_ok=True)

# create tmp/ dir to store intermediate results
tmp_dir = out_dir + '/tmp'
if not os.path.exists(tmp_dir):     
	os.makedirs(tmp_dir, exist_ok=True)

scatter_dir = out_dir + "/scatter" 
if not os.path.exists(scatter_dir):         
	os.makedirs(scatter_dir, exist_ok=True)
# ===============================================================


# >>>>>>>>>>>> Read mutability rates by 7- or 3-nucleotide into a data frame <<<<<<<<<<<<
kmer = 7 # or 3

# k-mer = 3
mut_3mer_matrix_file = '../mutability_matrices/mutation_rate_by_trinucleotide_matrix.txt'
mut_3mer_matrix = pd.read_csv(mut_3mer_matrix_file, sep='\t', header=0, index_col=False)

mut_3mer_matrix['sum'] = mut_3mer_matrix.loc[:, 'A'] + mut_3mer_matrix.loc[:, 'T'] + mut_3mer_matrix.loc[:, 'C'] + mut_3mer_matrix.loc[:, 'G']
#print(mut_3mer_matrix.head())
mut_3mer_matrix_dict = dict(zip(mut_3mer_matrix['trint'], mut_3mer_matrix['sum']))
#print(mut_3mer_matrix_dict)


# k-mer = 7
mut_7mer_matrix_file = '../mutability_matrices/heptamer_mutability_rates.processed_sums.txt'
mut_7mer_matrix = pd.read_csv(mut_7mer_matrix_file, sep=',', header=0, index_col=False)
#print(mut_7mer_matrix.head())
mut_7mer_weighted_sum_dict = dict(zip( mut_7mer_matrix['ref_seq'], mut_7mer_matrix['weighted_sum']))
mut_7mer_sum_dict = dict(zip( mut_7mer_matrix['ref_seq'], mut_7mer_matrix['sum']))
#print(mut_7mer_sum_dict)


mut_rate_dict = dict()
if kmer == 7:
	mut_rate_dict = mut_7mer_weighted_sum_dict
elif kmer == 3:
	mut_rate_dict = mut_3mer_matrix_dict



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



df = pd.read_table(simplified_vcf_file, low_memory=False)
df.columns = ['POS', 'REF', 'ALT', 'QUAL', 'AC', 'AF', 'AN', 'GQ']

# => NEW: filter variants with median Genotype Quality < 20
print('Filtering variants with median GQ < 20 ...')
print(df.shape)
df = df.loc[ df['GQ'] >= 20, :]
print(df.shape)
print(df.head())



# >> VCF Table Pre-Processing
df['LEN_DIFF'] = df['REF'].str.len() - df['ALT'].str.len()

# > Keep only SNVs - Filter out indels
print(df.shape)
df = df.loc[ df['LEN_DIFF'] == 0, :]
print(df.shape)


df = df[['POS', 'ALT', 'AC', 'AF', 'AN']]
print(df.info())


# cleanup rows with missing AF values (if any)
print(df.shape)
df.loc[ df.AF.astype(str) == '.', 'AF'] = -1
print(df.shape)

df['AF'] = df['AF'].astype(float)

# => NEW: remove variants with AF=0
df = df.loc[ df.AF != 0]
print(df.shape)



start_idx = df['POS'].iloc[0]
end_idx = df['POS'].iloc[-1]

chr_range = end_idx - start_idx + 1

# >> New code to allow for comparison of respective windows
chr_first_window_idx = int(start_idx / win_len)
first_win_offset = start_idx - (chr_first_window_idx * win_len)


total_num_windows = int(chr_range / win_len) + 1 # DEBUG: added '+ 1' --> to be tested: seems to be correct

print('Start index:', start_idx, ' | End index:', end_idx)
print('Chromosome range:', chr_range)
print('Num. of rows:', str(df.shape[0]))
print('First window index:', str(chr_first_window_idx))
print('First window offset:', str(first_win_offset) + ' nt')
print('Total genomic windows to scan:', total_num_windows)



# record start coordinate of first variant (exlcuding indels) at each chromosome
chr_start_coords_file = out_dir +'/chr_start_coords.txt'

## TO-DO: need to update other scripts that read from chr_start_coords.txt [DONE]
tmp_fh = open(chr_start_coords_file, 'a')
tmp_fh.write(chr + '\t' + str(start_idx) + '\t' + str(chr_first_window_idx) + '\n')
tmp_fh.close()



## ** Essential df pre-processing to speed-up parsing **
print('--------------------')
#df['POS'] = df['POS'] - start_idx # [DEPRECATED] in new version
#df['POS'] = df['POS'] - first_win_offset # [BETA] in new version
df['WIN'] = (df['POS'] / win_len).astype(int)


# Expected to be lower than 'total_num_windows' due to lack of variants or coverage inadequacies in some genomic regions - CHECKED OK
print('len(df[WIN].unique): ' + str(len(df['WIN'].unique())))
print(df.head(20))



placeholder_val = -1

def get_collapsed_counts(df, placeholder_val):

	print(df.head())
	print(df.tail())


	## fill in windows with no variants in the original VCF file with a placeholder value
	ac_collapsed_df = df.groupby(['WIN'])['AC'].agg('mean')
	# all_ac_df = pd.DataFrame(pd.Series(ac_collapsed_df))

	af_collapsed_df = df.groupby(['WIN'])['AF'].agg('mean')
	# all_af_df = pd.DataFrame(pd.Series(af_collapsed_df))
	# print(all_af_df.head())


	all_var_collapsed_df = df.groupby(['WIN']).agg(['count'])
	#all_variants_df = pd.DataFrame(pd.Series(all_var_collapsed_df.iloc[:, 0]), pd.Series(ac_collapsed_df), pd.Series(af_collapsed_df))
	all_variants_df = pd.concat([pd.Series(all_var_collapsed_df.iloc[:, 0]), pd.Series(ac_collapsed_df), pd.Series(af_collapsed_df)], axis=1)
	all_variants_df.columns = ['count', 'ac', 'af']
	print(all_variants_df.shape)
	print(all_variants_df.head())


	bins = [[1,1], [2,5], [6,10], [11,50], [51,200], [201,10000000]]
	#bins = [[1,5], [6,10], [11,50], [51,200], [201, 500], [501,10000000]]

	bin_cnt = 1
	for b in bins:
		low = b[0]
		high = b[1]

		print('===================================')
		print('bin: [', low, '-', high, ']')
		subset = df.loc[(df['AC'] >= low) & (df['AC'] <= high)]
		print(subset.head())

		grouped_subset = subset.groupby(['WIN'])['AC'].agg('count')
		print(type(grouped_subset))
		print(grouped_subset.head())
		print(grouped_subset.tail())

		all_variants_df = pd.concat([all_variants_df, grouped_subset.rename('bin_' + str(bin_cnt))], axis=1, )
		all_variants_df.fillna(0, inplace=True)
		print(all_variants_df.head())

		bin_cnt += 1


	all_variants_df.index.name = None

	## ===== Scale Allele Counts (either AC or plain numbers of variants) =====
	# > scale by AF
	# all_variants_df['count'] = all_variants_df['count'] * pd.Series(af_collapsed_df)

	# > scale by AC
	# all_variants_df['count'] = all_variants_df['count'] * pd.Series(ac_collapsed_df)

	print(all_variants_df.head())



	#first_index = all_variants_df.index[0]
	#print('first index: ' + str(first_index))
	all_variants_df.index -= chr_first_window_idx
	#print(all_variants_df.head())
	#print(all_variants_df.tail())

	print(all_variants_df.head())


	all_variant_indexes = all_variants_df.index
	print('all_variant_indexes:',all_variant_indexes[:20])
	print('all_variant_indexes:', len(all_variant_indexes))
	print('all_variant_indexes:',all_variant_indexes[-10:])


	all_window_indexes = np.arange(total_num_windows + 1)
	print('all_window_indexes:', all_window_indexes)
	zero_window_indexes = np.setdiff1d(all_window_indexes, all_variant_indexes)
	#print('zero_window_indexes:',zero_window_indexes)
	print('zero_window_indexes:', len(zero_window_indexes))


	tmp_zero_variants_df = pd.DataFrame(placeholder_val, index=zero_window_indexes, columns=['tmp_count'])
	zero_variants_df = pd.concat([tmp_zero_variants_df, tmp_zero_variants_df, tmp_zero_variants_df], axis=1)
	zero_variants_df.columns = ['count', 'ac', 'af']
	print(zero_variants_df.head())
	print(zero_variants_df.tail())


	concat_df = pd.concat([all_variants_df, zero_variants_df])
	concat_df = concat_df.sort_index()
	print(concat_df.head(20))
	print(concat_df.tail(10))
	print(concat_df.shape)

	return concat_df



## >>>> Calculate mutability rate, GC conente and CpG dinucleotides for each window <<<<
# Calculate aggreagate or average mutability rate, GC-dontent and CpG dinucleotides for each window
# and return a series with all window indexes and the calculated mut. rate, GC and GpG counts.
def get_expected_mutability_by_trimer_per_window(df, placeholder_val):
	
	agg_mut_rates_per_window = dict()
	gc_content_per_window = dict()
	cpg_per_window = dict()

	print('win_len:', win_len)

	valid_window_indexes = df.loc[:, 'WIN'].unique()
	print('Num of valid windows:', len(valid_window_indexes))	


	print(df.head())
	print(df.tail())



	cc = 0
	for win in valid_window_indexes:
		seq_start = win * win_len
		seq_end = seq_start + win_len - 1

		#print('Seq start:', seq_start, ' Seq end:', seq_end)

		cmd = './util/twoBitToFa ' + human_ref_genome_2bit + ':' + chr + ':' + str(seq_start) + '-' + str(seq_end) + " /dev/stdout | grep -v '>' | tr '\n' ' ' | sed 's/ //g'"
		#print(cmd)
		p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

		stdout, stderr = p.communicate()
		seq = str(stdout, "utf-8").rstrip()
		#print('Seq:', seq)


		# mutability rate
		aggregate_rate = 0
		for k in range(0, kmer):
		    sub_seqs = split_seq_into_same_size_susbseq(seq[k:], kmer)
		    aggregate_rate += sum([mut_rate_dict[s] for s in sub_seqs if 'N' not in s])
		#print(aggregate_rate)


		# CpG dinucleotides
		cpg = seq.count('CG')
		#print('CpG:', cpg)

		
		# gc content
		gc_content = seq.count('G') + seq.count('C')


		
		win_idx = win - chr_first_window_idx
		agg_mut_rates_per_window[ win_idx ] = aggregate_rate
		cpg_per_window[ win_idx ] = cpg
		gc_content_per_window[ win_idx ] = gc_content
		
		cc += 1
		if cc % 100 == 0:
			print('Current window:', cc , 'out of', len(valid_window_indexes))
		#if(cc > 30):
		#	break
		#	sys.exit()
		
	#print(agg_mut_rates_per_window)


	# mutability-rate
	agg_mut_rates_per_window_df = pd.DataFrame.from_dict(agg_mut_rates_per_window, orient='index')
	agg_mut_rates_per_window_df.columns = ['mut_rate']
	print(agg_mut_rates_per_window_df.head(10))
	
	# GC content
	gc_content_per_window_df = pd.DataFrame.from_dict(gc_content_per_window, orient='index')
	gc_content_per_window_df.columns = ['gc_content']
	print(gc_content_per_window_df.head(10))

	# CpG dinucleotides
	cpg_per_window_df = pd.DataFrame.from_dict(cpg_per_window, orient='index')
	cpg_per_window_df.columns = ['cpg']
	print(cpg_per_window_df.head(10))


	if len(list(set(agg_mut_rates_per_window_df.index) - set(gc_content_per_window_df.index))) != 0: 
		raise ValueError('agg_mut_rates_per_window_df.index != gc_content_per_window_df.index')
		sys.exit()

	overall_features_df = pd.concat([agg_mut_rates_per_window_df, gc_content_per_window_df, cpg_per_window_df], axis=1)
	print(overall_features_df.head())


	# make start-window-index: 0
	offset_valid_window_indexes = np.array(valid_window_indexes) - chr_first_window_idx


	all_window_indexes = np.arange(total_num_windows + 1)
	print('all_window_indexes:', len(all_window_indexes))
	zero_window_indexes = np.setdiff1d(all_window_indexes, offset_valid_window_indexes)
	print('zero_window_indexes:', len(zero_window_indexes))
	print('zero_window_indexes:',zero_window_indexes[:20], '...')
	print('valid_window_indexes:', offset_valid_window_indexes[:20], '...')	



	zero_variants_df = pd.DataFrame(placeholder_val, index=zero_window_indexes, columns=['mut_rate', 'gc_content', 'cpg'])
	print(zero_variants_df.head())
	print(zero_variants_df.tail())


	concat_df = pd.concat([overall_features_df, zero_variants_df])
	concat_df = concat_df.sort_index()
	print(concat_df.head(30))
	print(concat_df.tail(30))
	print(concat_df.shape)
		
	return concat_df



def split_seq_into_same_size_susbseq(seq, size):
	return [''.join(x) for x in zip(*[list(seq[z::size]) for z in range(size)])]





# Skip mutability rate calculations per window if they have already been calculated
overall_features_df_file = tmp_dir + '/chr' + chr + '.overall_features_df.csv'
overall_features_df = None



if Path(overall_features_df_file).exists():
	print("\n>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<\n")
	print("overall_features_df_file exists!")
	print("\n>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<\n")

	overall_features_df = pd.read_csv(overall_features_df_file, index_col=0)
else:
	overall_features_df = get_expected_mutability_by_trimer_per_window(df, placeholder_val)
	overall_features_df.to_csv(overall_features_df_file)


sys.exit()



# > Get all variants across all windows
#print('__all_variants_df__')
all_variants_df = get_collapsed_counts(df, placeholder_val)


#print(all_variants_df.shape)
#print(df.head(3))
#print(df.tail(3))


# ________ Get common variants across all windows _________
print('__common_variants_df__')
tmp_com_var_df = df.loc[ df['AF'] >= MAF_thres, ]


### >>>>>> BETA: numerical-cutoff for common variants calculation (instead of MAF) <<<<<<<
#ALLELE_NUM_CUTOFF = 1
#tmp_com_var_df = df.loc[ df['AC'] > ALLELE_NUM_CUTOFF, ]

### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

print(tmp_com_var_df.head(20))
print(tmp_com_var_df.tail(10))

print("\n\n---------------------------------------\n\n")
common_variants_df = get_collapsed_counts(tmp_com_var_df, 0)   # <<<<-------- BETA: try that with -1 as placeholder [DONE: signal vanishes]


#print(common_variants_df.shape)

print(common_variants_df.head())
print(common_variants_df.tail())




rvis_score_quotients = common_variants_df['count'] / all_variants_df['count']
rvis_score_quotients = rvis_score_quotients.fillna(-1)
rvis_score_quotients.columns = ['ratio']
print(rvis_score_quotients.head())

rvis_score_quotients.to_csv(var_ratios_dir + '/common_all_var_ratios_chr' + chr + '.csv')


# Create linear regression object 
#X = list(all_variants_df['count'])
all_variants = list(all_variants_df['count'])
all_ac = list(all_variants_df['ac'])
all_af = list(all_variants_df['af'])


### >>>>>> BETA: shuffle Xy <<<<<<<<
#shuffle(X)
### >>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<


common_variants = list(common_variants_df['count'])  # <-- y in regression

### >>>>>> BETA: shuffle Xy <<<<<<<<
#shuffle(y)
### >>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<

idx = list(all_variants_df.index)

#tmp_df = pd.DataFrame({'idx':idx, 'X': X, 'y':y})
#tmp_df = pd.DataFrame({'idx': idx, 'y': common_variants, 'all_variants': all_variants, 'all_ac': all_ac, 'all_af': all_af})
tmp_df = pd.DataFrame({'idx': idx, 'y': common_variants, 'common_div_by_all_var_ratio': rvis_score_quotients})
tmp_df = pd.concat([tmp_df, all_variants_df], axis=1)
tmp_df = tmp_df.rename(columns = {'count': 'all_variants'})

tmp_df = pd.concat([tmp_df, overall_features_df], axis=1)
print(tmp_df.head())
print(tmp_df.tail())

#tmp_df = tmp_df[ ['idx', 'X', 'y'] ]
#print(len(tmp_df))
#print(tmp_df.head(20))
#print(tmp_df.tail(20))



xy_file = tmp_dir + '/Xy.chr' + chr + '.txt'
tmp_df.to_csv(xy_file, index=False, line_terminator='\r\n')




## --deprecated
if generate_intermediate_plots:
	rvis_scores_arr = np.array(rvis_scores)
	print(min(rvis_scores_arr))
	print(max(rvis_scores_arr))

	# Plot RVIS scores for current chromosome
	f1 = plt.figure()
	plt.plot(rvis_scores_arr, linewidth=0.1)
	f1.suptitle('gwRVIS values across chromosome ' + chr, fontsize=12) 
	plt.xlabel('chr window index (genomic coordinate)', fontsize=10) 
	plt.ylabel('gwRVIS score', fontsize=10)
	#plt.show()

	rvis_score_quotients = np.array(rvis_score_quotients)
	f1a = plt.figure()
	plt.plot(rvis_score_quotients, linewidth=0.1)
	f1a.suptitle('common / all variants ratios across chromosome ' + chr, fontsize=12) 
	plt.xlabel('chr window index (genomic coordinate)', fontsize=12) 
	plt.ylabel('common / all variants ratio', fontsize=10)
	plt.ylim((-1.2, 1.2))

	binwidth = 0.001
	print(len(rvis_scores_arr))
	rvis_scores_arr = rvis_scores_arr[~np.isnan(rvis_scores_arr) ]
	print(len(rvis_scores_arr))
	f2 = plt.figure()
	plt.hist(rvis_scores_arr, bins=np.arange(min(rvis_scores_arr), max(rvis_scores_arr) + binwidth, binwidth))
	f2.suptitle('Histogram of gwRVIS values across chromosome ' + chr +'\n (excluding regions with no variation data)', fontsize=12) 
	plt.xlabel('chr window index (genomic coordinate)', fontsize=10) 
	plt.ylabel('Count', fontsize=10)


	# CDF plot 
	f3 = plt.figure()
	rvis_scores_arr = rvis_scores_arr[ rvis_scores_arr != 0 ]
	plt.hist(rvis_scores_arr, bins=np.arange(min(rvis_scores_arr), max(rvis_scores_arr) + binwidth, binwidth), cumulative=True, normed=True, histtype='step', alpha=0.55, color='purple')
	f3.suptitle('CDF of gwRVIS values across chromosome ' + chr + '\n (excluding regions with no variation data)', fontsize=12) 
	plt.xlabel('chr window index (genomic coordinate)', fontsize=10) 
	plt.ylabel('Count', fontsize=10)

	pp = PdfPages(plots_dir + "/rvis_chr" + chr + ".pdf")
	pp.savefig(f1)
	pp.savefig(f1a)
	pp.savefig(f2)
	pp.savefig(f3)
	pp.close()



print('Elapsed time (hh:mm:ss):', datetime.now() - startTime)
