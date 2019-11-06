import matplotlib 
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from datetime import datetime
import sys
import os
import re
import ntpath
from subprocess import call
from random import shuffle
import subprocess
from pathlib import Path
import multiprocessing 

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from custom_utils import create_out_dir, get_config_params


def pre_process_vcf_table(filtered_vcf, variant_filter=''):
	
	if config_params['filter_ccds_variants']:
		filtered_vcf_basename = ntpath.basename(filtered_vcf)
		print(filtered_vcf_basename)

		os.system("tail -n+2 " + filtered_vcf + ' | awk -v chrom=' + chrom + """ '{print "chr"chrom"\t"$1-1"\t"$1"\t"$2","$3","$4","$5","$6","$7","$8}' | subtractBed -a stdin  -b """ + genomic_classes_files['ccds'] + """ | awk -F"\t|," '{print $3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > """ + tmp_dir+"/"+filtered_vcf_basename+".tmp")

		os.system("cat " + filtered_vcf + " | head -1 > " + tmp_dir+"/"+filtered_vcf_basename+".header 2>/dev/null")
		print("[Broken pipe message may be safely ignored]")

		os.system("cat " + tmp_dir+"/"+filtered_vcf_basename+".header " + tmp_dir+"/"+filtered_vcf_basename+".tmp > " + tmp_dir+"/"+filtered_vcf_basename+".no_ccds")
		os.system("rm " +  tmp_dir+"/"+filtered_vcf_basename+".tmp " +  tmp_dir+"/"+filtered_vcf_basename+".header")
		print(tmp_dir+"/"+filtered_vcf_basename+'.no_ccds')
		df = pd.read_csv( tmp_dir+"/"+filtered_vcf_basename+".no_ccds", low_memory=False, sep='\t')
	else:
		df = pd.read_csv(filtered_vcf, low_memory=False, sep='\t')
	print(df.shape)
	print(df.head())


	# ===============  VCF Table Pre-Processing  ===============
	df['LEN_DIFF'] = df['REF'].str.len() - df['ALT'].str.len()
	# > Keep only SNVs - Filter out indels
	if variant_filter == 'snv':
		print('\n- Keeping SNVs only...')
		df = df.loc[ df['LEN_DIFF'] == 0, :]
	elif variant_filter == 'indels':
		print('\n- Keeping INDELs only...')
		df = df.loc[ df['LEN_DIFF'] != 0, :]
	print(df.shape)


	# Flag rows with missing AF values (if any)
	# TO-DO: why was I replacing it with -1? I should replace with 0 probably!
	df.loc[ df.AF.astype(str) == '.', 'AF'] = 0 
	df['AF'] = df['AF'].astype(float)

	# Clenaup rows with AF=0
	df = df.loc[ df.AF != 0]
	print(df.shape)
	print(df.info())


	# ------------ Window index arithmetic ----------------
	start_idx = df['POS'].iloc[0]
	end_idx = df['POS'].iloc[-1]

	chr_range = end_idx - start_idx + 1

	# >> Window index calculation (0-based) to allow comparison of respective windows from different datasets (e.g. gnomad vs topmed)
	chr_first_window_idx = int(start_idx / win_len)
	first_win_offset = start_idx - (chr_first_window_idx * win_len)


	# Calculate all full 'win_len' windows and add another 2 for the flanking windows at the start and end
	total_num_windows = (int(end_idx / win_len) * win_len - ((chr_first_window_idx + 1) * win_len)) / win_len + 2

	print('Start index:', start_idx, ' | End index:', end_idx)
	print('Chromosome range:', chr_range)
	print('Num. of rows:', str(df.shape[0]))
	print('First window index:', str(chr_first_window_idx))
	print('First window offset:', str(first_win_offset) + ' nt')
	print('Total genomic windows to scan:', total_num_windows)


	# Record start coordinate of first variant (exlcuding indels) at each chromosome
	chr_start_coords_file = out_dir +'/chr_start_coords.txt'
	tmp_fh = open(chr_start_coords_file, 'a')
	tmp_fh.write(chrom + '\t' + str(start_idx) + '\t' + str(chr_first_window_idx) + '\n')
	tmp_fh.close()



	## Essential df pre-processing to speed-up parsing for feature extraction
	df['WIN'] = (df['POS'] / win_len).astype(int)

	return df, chr_first_window_idx, total_num_windows




def get_mutability_rates(kmer=7):
	"""
		Read mutability rates by 7- or 3-nt and return them in a dictionary
	"""
	# k-mer = 3
	mut_3mer_matrix_file = '../other_datasets/mutability_matrices/mutation_rate_by_trinucleotide_matrix.txt'
	mut_3mer_matrix = pd.read_csv(mut_3mer_matrix_file, sep='\t', header=0, index_col=False)

	mut_3mer_matrix['sum'] = mut_3mer_matrix.loc[:, 'A'] + mut_3mer_matrix.loc[:, 'T'] + mut_3mer_matrix.loc[:, 'C'] + mut_3mer_matrix.loc[:, 'G']
	mut_3mer_matrix_dict = dict(zip(mut_3mer_matrix['trint'], mut_3mer_matrix['sum']))


	# k-mer = 7
	mut_7mer_matrix_file = '../other_datasets/mutability_matrices/heptamer_mutability_rates.processed_sums.txt'
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

	return mut_rate_dict



def get_regulatory_features_per_window(chrom, seq_start, seq_end):
	
	manager = multiprocessing.Manager()
	regul_features = manager.dict()

	tmp_bed_file = scratch_dir + '/' + str(chrom) + '.' + str(seq_start) + '_' + str(seq_end) + '.bed'
	with open(tmp_bed_file, 'w') as fh:
		fh.write(str(chrom) + '\t' + str(seq_start) + '\t' + str(seq_end) + '\n')


	def call_intersect_bed(ensembl_regul_file, regul_features):

		cmd = 'cat ' + tmp_bed_file + ' | intersectBed -a stdin -b ' + ensembl_regul_file
		#print(cmd)
		
		p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout, stderr = p.communicate()
		res = str(stdout, "utf-8").rstrip()

		# Also taking care of intersectBed results with > 1 overlaps in the results
		tmp_interval = 0
		for line in res.split('\n'):
			if line == '':
				continue
			#print(line)
			chrom, start, end = line.split('\t')
			tmp_interval += int(end) - int(start) + 1

		cur_feature = tmp_interval / win_len
		regul_elem = ntpath.basename(ensembl_regul_file).replace('.bed', '')
		regul_features[regul_elem] = cur_feature

	processes = []
	for ensembl_regul_file in ensembl_regulatory_files:
		p = multiprocessing.Process(target=call_intersect_bed, args=(ensembl_regul_file, regul_features))
		processes.append(p)
		p.start()	
	
	for p in processes:
		p.join()

	return dict(regul_features)



def extract_additional_features_per_window(df, chr_first_window_idx, total_num_windows, mut_rate_dict, placeholder_val=-1):
	"""
		Calculate aggreagate or average mutability rate, GC-content and CpG dinucleotides for each window
		and return a series with all window indexes and the calculated mut. rate, GC and GpG counts.
	"""
	
	# Skip mutability rate calculations per window if they have already been calculated
	additional_features_df_file = tmp_dir + '/chr' + chrom + '.additional_features_df.csv'

	# TEMP: comment next 4 lines for DEBUG
	if Path(additional_features_df_file).exists():
		print(">> additional_features_df_file already exists. Reading file...")
		additional_features_df = pd.read_csv(additional_features_df_file, index_col=0)
		return additional_features_df


	agg_mut_rates_per_window = dict()
	gc_content_per_window = dict()
	cpg_per_window = dict()
	cpg_islands_per_window = dict()
	regul_features_per_window = dict()

	print('win_len:', win_len)

	valid_window_indexes = df.loc[:, 'WIN'].unique()
	print('Num of valid windows:', len(valid_window_indexes))	



	cc = 0
	for win in valid_window_indexes:
		seq_start = win * win_len
		seq_end = seq_start + win_len - 1

		#print(chrom + ':' + str(seq_start) + '-' + str(seq_end))

		# Get regulatory features (open chromatin, CTCF, TG binding, etc.) for current window
		# DEBUG - move to the end of pipeline
		#regul_features = get_regulatory_features_per_window(chrom, seq_start, seq_end)
		#print(regul_features)

		cmd = '../utils/twoBitToFa ' + human_ref_genome_2bit + ':' + chrom + ':' + str(seq_start) + '-' + str(seq_end) + " /dev/stdout | grep -v '>' | tr '\n' ' ' | sed 's/ //g'"
		#print(cmd)
		p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

		stdout, stderr = p.communicate()
		seq = str(stdout, "utf-8").rstrip()
		#print('Seq:', seq)
		stderr_str = str(stderr, "utf-8").rstrip()
		if stderr_str != '':
			print(stderr_str)


		# Mutability rate
		aggregate_rate = 0
		for k in range(0, kmer):
		    sub_seqs = split_seq_into_same_size_susbseq(seq[k:], kmer)
		    aggregate_rate += sum([mut_rate_dict[s] for s in sub_seqs if 'N' not in s])
		#print('Mut. rate (aggregate):', aggregate_rate)


		# CpG dinucleotides
		cpg = seq.count('CG')
		#print('CpG:', cpg)

		# CpG islands (>=2 CpG dinucleotides, e.g.: CGCG, CGCGCGCG, etc.)
		cpg_islands = len(re.findall(re.compile('CG(CG)+'), seq))
		#print('CpG islands (CGCG, CGCGCG, ...):', cpg_islands)
		
		
		# GC content (ratio over the entire window size)
		gc_content = (seq.count('G') + seq.count('C')) / win_len
		#print('GC ratio:', gc_content)


		win_idx = win - chr_first_window_idx
		agg_mut_rates_per_window[ win_idx ] = aggregate_rate
		cpg_per_window[ win_idx ] = cpg
		cpg_islands_per_window[ win_idx ] = cpg_islands
		gc_content_per_window[ win_idx ] = gc_content
		# DEBUG - move to the end of pipeline
		#regul_features_per_window[ win_idx ] = regul_features

		cc += 1
		#if cc % 100 == 0:
		if cc % 10 == 0:
			#break
			print('Chr' + str(chrom) + ' - Current window:', cc , 'out of', len(valid_window_indexes))
		

	# Mutability-rate
	agg_mut_rates_per_window_df = pd.DataFrame.from_dict(agg_mut_rates_per_window, orient='index')
	agg_mut_rates_per_window_df.columns = ['mut_rate']
	print(agg_mut_rates_per_window_df.head(10))
	
	# CpG di-nucleotides
	cpg_per_window_df = pd.DataFrame.from_dict(cpg_per_window, orient='index')
	cpg_per_window_df.columns = ['cpg']
	print(cpg_per_window_df.head(10))

	# CpG islands
	cpg_islands_per_window_df = pd.DataFrame.from_dict(cpg_islands_per_window, orient='index')
	cpg_islands_per_window_df.columns = ['cpg_islands']
	print(cpg_islands_per_window_df.head(10))

	# GC content
	gc_content_per_window_df = pd.DataFrame.from_dict(gc_content_per_window, orient='index')
	gc_content_per_window_df.columns = ['gc_content']
	print(gc_content_per_window_df.head(10))

	# Regulatory features
	# DEBUG - move to the end of pipeline
	#regul_features_df = pd.DataFrame(regul_features_per_window).T
	#print(regul_features_df.head(10))


	# DEBUG - move to the end of pipeline
	#additional_features_df = pd.concat([agg_mut_rates_per_window_df, gc_content_per_window_df, cpg_per_window_df, cpg_islands_per_window_df, regul_features_df], axis=1)
	additional_features_df = pd.concat([agg_mut_rates_per_window_df, gc_content_per_window_df, cpg_per_window_df, cpg_islands_per_window_df], axis=1)
	print(additional_features_df.head())
	print(additional_features_df.tail())

	# make start-window-index: 0
	offset_valid_window_indexes = np.array(valid_window_indexes) - chr_first_window_idx


	all_window_indexes = np.arange(total_num_windows)
	#print('all_window_indexes:', len(all_window_indexes))
	zero_window_indexes = np.setdiff1d(all_window_indexes, offset_valid_window_indexes)
	#print('zero_window_indexes:', len(zero_window_indexes))
	#print('zero_window_indexes:', zero_window_indexes[:10], '...')
	#print('valid_window_indexes:', offset_valid_window_indexes[:10], '...')	



	zero_variants_df = pd.DataFrame(placeholder_val, index=zero_window_indexes, columns=additional_features_df.columns.values)
	#print(zero_variants_df.head())
	#print(zero_variants_df.tail())


	additional_features_df = pd.concat([additional_features_df, zero_variants_df])
	additional_features_df = additional_features_df.sort_index()
	print(additional_features_df.head(10))
	print(additional_features_df.tail(10))

	additional_features_df.to_csv(additional_features_df_file)
	

		
	return additional_features_df




def split_seq_into_same_size_susbseq(seq, size):
	return [''.join(x) for x in zip(*[list(seq[z::size]) for z in range(size)])]




def get_collapsed_counts(df, placeholder_val=-1):

	print(df.head())
	print(df.tail())
	

	## fill in windows with no variants in the original VCF file with a placeholder value
	# ??

	mean_ac_collapsed_df = df.groupby(['WIN'])['AC'].agg('mean')
	mean_af_collapsed_df = df.groupby(['WIN'])['AF'].agg('mean')

	all_var_collapsed_df = df.groupby(['WIN']).agg(['count'])
	#all_variants_df = pd.DataFrame(pd.Series(all_var_collapsed_df.iloc[:, 0]), pd.Series(ac_collapsed_df), pd.Series(af_collapsed_df))
	all_variants_df = pd.concat([pd.Series(all_var_collapsed_df.iloc[:, 0]), pd.Series(mean_ac_collapsed_df), pd.Series(mean_af_collapsed_df)], axis=1)
	all_variants_df.columns = ['count', 'mean_ac', 'mean_af']
	print(all_variants_df.shape)
	print(all_variants_df.head())
	print(all_variants_df.tail())


	bins = [[1,1], [2,5], [6,10], [11,50], [51,200], [201,10000000]]
	#bins = [[1,5], [6,10], [11,50], [51,200], [201, 500], [501,10000000]]

	bin_cnt = 1
	for b in bins:
		low, high = b[0], b[1]

		print('--------------')
		print('\n> bin:', b)
		subset = df.loc[(df['AC'] >= low) & (df['AC'] <= high)]
		print(subset.head())
		print(subset.shape)

		grouped_subset = subset.groupby(['WIN'])['AC'].agg('count')
		print(grouped_subset.head())
		print(grouped_subset.tail())
		print(grouped_subset.shape)

		all_variants_df = pd.concat([all_variants_df, grouped_subset.rename('bin_' + str(bin_cnt))], axis=1, )
		all_variants_df.fillna(0, inplace=True)

		bin_cnt += 1

	all_variants_df.index.name = None
	all_variants_df.index -= chr_first_window_idx
	print(all_variants_df.head())
	

	## (Deprecated) ===== Scale Allele Counts (either AC or plain numbers of variants) =====
	# > scale by AF
	# all_variants_df['count'] = all_variants_df['count'] * pd.Series(af_collapsed_df)

	# > scale by AC
	# all_variants_df['count'] = all_variants_df['count'] * pd.Series(ac_collapsed_df)


	all_variant_indexes = all_variants_df.index
	print('all_variant_indexes:', len(all_variant_indexes))

	all_window_indexes = np.arange(total_num_windows)
	print('all_window_indexes:', len(all_window_indexes))
	zero_window_indexes = np.setdiff1d(all_window_indexes, all_variant_indexes)
	print('zero_window_indexes:', len(zero_window_indexes))


	zero_variants_df = pd.DataFrame(placeholder_val, index=zero_window_indexes, columns=all_variants_df.columns.values)


	concat_variant_cnt_df = pd.concat([all_variants_df, zero_variants_df])
	concat_variant_cnt_df = concat_variant_cnt_df.sort_index()

	return concat_variant_cnt_df


def get_common_and_all_variants(df):

	# > Get all variants across all windows
	print('> Getting all variants ...')
	all_variants_df = get_collapsed_counts(df)
	print(all_variants_df.shape)
	print(all_variants_df.head())


	# > Get common variants across all windows
	print('> Getting common variants ...')
	tmp_com_var_df = df.loc[ df['AF'] >= MAF_thres, ]
	print(tmp_com_var_df.shape)
	print(tmp_com_var_df.head())

	## TO-DO: why is this not -1 as in the all_variants_df calculations [DONE: signal vanishes]
	common_variants_df = get_collapsed_counts(tmp_com_var_df, placeholder_val=0)   
	print(common_variants_df.shape)
	print(common_variants_df.head())


	# Get common / all variant ratios
	gwrvis_score_quotients = common_variants_df['count'] / all_variants_df['count']
	gwrvis_score_quotients = gwrvis_score_quotients.fillna(-1)
	gwrvis_score_quotients = gwrvis_score_quotients.to_frame()
	gwrvis_score_quotients.columns = ['common_vs_all_variants_ratio']
	print(gwrvis_score_quotients.head())
	print(type(gwrvis_score_quotients))
	gwrvis_score_quotients.to_csv(var_ratios_dir + '/common_vs_all_variant_ratios.chr' + chrom + '.csv')


	return common_variants_df, all_variants_df, gwrvis_score_quotients


def prepare_data_for_regression(common_variants_df, all_variants_df, gwrvis_score_quotients, additional_features_df):

	print(common_variants_df.head())
	print(all_variants_df.head())

	# Create linear regression object 
	#all_variants = list(all_variants_df['count'])  # X in regression: passed on with the entire all_variants_df data frame

	allv_mean_ac = list(all_variants_df['mean_ac'])
	allv_mean_af = list(all_variants_df['mean_af'])


	idx = list(all_variants_df.index)
	common_variants = list(common_variants_df['count'])  # y in regression
	gwrvis_ratios = list(gwrvis_score_quotients['common_vs_all_variants_ratio'])


	tmp_df = pd.DataFrame({'idx': idx, 'y': common_variants, 'common_vs_all_variants_ratio': gwrvis_ratios})
	tmp_df = pd.concat([tmp_df, all_variants_df], axis=1)
	tmp_df = tmp_df.rename(columns = {'count': 'all_variants'})

	tmp_df = pd.concat([tmp_df, additional_features_df], axis=1)
	print(tmp_df.head())
	print(tmp_df.tail())

	xy_file = tmp_dir + '/Xy.chr' + chrom + '.txt'
	tmp_df.to_csv(xy_file, index=False, line_terminator='\r\n')





if __name__ == '__main__':

	startTime = datetime.now()

	args = sys.argv
	chrom = args[1]
	config_file = args[2] #'config.yaml'



	# Read run parameters from config file and store into a dictionary
	config_params = get_config_params(config_file)
	print(config_params)
	hg_version = config_params['hg_version']
	grch = {'hg19': '37', 'hg38': '38'}

	genomic_classes_files = {}
	print('cwd:', os.getcwd())

	with open(config_params['genomic_classes']) as fh:
		for line in fh:
			line = line.rstrip()
			genomic_class, cur_path, _ = line.split('\t')
			genomic_classes_files[genomic_class] = cur_path


	# ==================== Initialisation ====================
	dataset = config_params['dataset']		# e.g. 'gnomad'
	population = config_params['population']	# e.g. 'all', 'FIN', etc.
	win_len = config_params['win_len']		# e.g. 250 (window length in nt)
	variant_filter = config_params['variant_filter']	# 'snv' or 'cnv', anything else retains all
	kmer = config_params['kmer']			# 7 or 3
	all_variants_upper_thres = config_params['all_variants_upper_thres']	# e.g. 200 (filter out windows with more than 200 variants before fitting regression)
	MAF_thres = config_params['MAF_thres']	        # e.g. 0.0001 (Minor Allele Frequency)
	variants_table_dir = config_params['variants_table_dir']	# gnomad-filtered_variant_tables-all-PASS_ONLY-NO_SEGDUP-NO_LCR-high_conf_regions
	# ----------------------


	human_ref_genome_2bit = '../' + hg_version + '/homo_sapiens_GRCh' + grch[hg_version] + '_FASTA/hsa' + grch[hg_version] +'.2bit'
	data_dir = '../' + dataset + '/out/' + variants_table_dir
	print('> data_dir: ' + data_dir)
	filtered_vcf = data_dir + '/chr' + chrom + '_' + dataset + '_table.' + population + '.txt.filtered'
	# ----------------------


	# Create out_MAF{threshold} dir to store output results (plots and gwRVIS csv files)
	out_dir = create_out_dir(config_file)
	print('> out_dir: ' + out_dir)

	var_ratios_dir = out_dir + '/var_ratios'
	if not os.path.exists(var_ratios_dir):     
		os.makedirs(var_ratios_dir, exist_ok=True)

	gwrvis_dir = out_dir + '/gwrvis_scores'
	if not os.path.exists(gwrvis_dir):     
		os.makedirs(gwrvis_dir, exist_ok=True)

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

	# Ensembl regulatory elements for additional features extraction
	scratch_dir = '../scratch'
	if not os.path.exists(scratch_dir):
		os.makedirs(scratch_dir)
	
	regulatory_elements = ['CTCF_binding_sites', 'Enhancers',
				'Open_chromatin', 'TF_binding_sites', 'H3K27ac',
				'H3K27me3', 'H4K20me1', 'H3K9ac', 'H3K4me1',	
				'H3K4me2', 'H3K4me3', 'H3K36me3']
	cell_line = 'Monocytes_CD14plus'

	ensembl_regulatory_files = []
	for elem in regulatory_elements:
		ensembl_regulatory_files.append('../ensembl/GRCh' + grch[hg_version] + '-Regulatory_Features/' + cell_line + '.' + elem + '.sorted.bed')
	# -------------------------------------------------------



	# ============================= Main Analysis ================================
	df, chr_first_window_idx, total_num_windows = pre_process_vcf_table(filtered_vcf, variant_filter=variant_filter)
	print(df.head())

	mut_rate_dict = get_mutability_rates(kmer=kmer)
	additional_features_df = extract_additional_features_per_window(df, chr_first_window_idx, total_num_windows, mut_rate_dict)

	common_variants_df, all_variants_df, gwrvis_score_quotients = get_common_and_all_variants(df)

	prepare_data_for_regression(common_variants_df, all_variants_df, gwrvis_score_quotients, additional_features_df)

	print('Elapsed time (hh:mm:ss):', datetime.now() - startTime)
