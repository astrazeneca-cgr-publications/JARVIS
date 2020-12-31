import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import subprocess
from sklearn.preprocessing import LabelEncoder, OneHotEncoder, MinMaxScaler
from sklearn.model_selection import train_test_split
import tensorflow as tf
import pickle
from scipy.io import savemat
import sys, os
import re
from multiprocessing import Process


from for_prediction_train_nn_model import train_with_cv

sys.path.insert(1, os.path.join(sys.path[0], '../..'))
import custom_utils



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
	



class JarvisDataPreprocessing:

	def __init__(self, config_file, input_features, chrom, NTHREADS=20):

		print("Initialising new JarvisDataPreprocessing object...")

		self.input_features = input_features
		self.chrom = chrom
		self.NTHREADS = NTHREADS
		
		
		
		# ==== Read config parameters ====
		config_params = custom_utils.get_config_params(config_file)
		self.hg_version = config_params['hg_version']
		print('\n\nhg_version:', self.hg_version)
		self.grch = {'hg19': '37', 'hg38': '38'}

		pathogenic_set = config_params['pathogenic_set']
		benign_set = config_params['benign_set']

		self.patho_benign_sets = pathogenic_set + '_' + benign_set
		self.win_len = config_params['win_len']
		#self.win_len = int(config_params['win_len'] / 2)
		self.Y_label = config_params['Y_label']

		# ==== Define dir structure ====
		out_dir = custom_utils.create_out_dir(config_file)

		self.ml_data_dir = out_dir + '/ml_data'
		if not os.path.exists(self.ml_data_dir):
			os.makedirs(self.ml_data_dir)
		self.seq_out_dir = self.ml_data_dir + '/raw_seq'
		if not os.path.exists(self.seq_out_dir):
			os.makedirs(self.seq_out_dir)
		self.feature_tables_dir = self.ml_data_dir + '/clinvar_feature_tables'
		if not os.path.exists(self.feature_tables_dir):
			os.makedirs(self.feature_tables_dir)
			
			
		self.jarvis_predictions_dir = self.ml_data_dir + '/jarvis_predictions'
		if not os.path.exists(self.jarvis_predictions_dir):
			os.makedirs(self.jarvis_predictions_dir)

		
		self.jarvis_predictions_per_chr_dir = self.jarvis_predictions_dir + '/chr' + str(self.chrom)
		if not os.path.exists(self.jarvis_predictions_per_chr_dir):
			os.makedirs(self.jarvis_predictions_per_chr_dir)
			

		# Specificy input (static) files
		self.human_ref_genome_2bit = '../' + self.hg_version +  '/homo_sapiens_GRCh' + self.grch[self.hg_version] + '_FASTA/hsa' + self.grch[self.hg_version] + '.2bit'







	def onehot_encoded_seq(self, seq):

		# remove any trailing new line characters
		seq = seq.rstrip()
		seq_array = np.array(list(seq))
		

		label_encoder = LabelEncoder()
		integer_encoded_seq = label_encoder.fit_transform(seq_array)
		integer_encoded_seq = integer_encoded_seq.reshape(len(integer_encoded_seq), 1)

		onehot_encoder = OneHotEncoder(sparse=False, categories='auto')
		onehot_seq = onehot_encoder.fit_transform(integer_encoded_seq) #.tolist()
		onehot_seq = np.transpose(onehot_seq)


		# Fill in rows in onehot encoded array with nts not present in the sequence
		nt_onehot_index = {'A': 0, 'T': 3, 'G': 2, 'C': 1}
		
		all_nts = set(nt_onehot_index.keys())
		present_nts = set(seq)
		absent_nts = list(all_nts - present_nts)

		#print('absent_nts:', absent_nts, 'present_nts:', present_nts)
		
		for nt, _ in sorted(nt_onehot_index.items(), key=lambda x: x[1]):
			if nt in absent_nts:
				#print('nt:', nt, 'nt_onehot_index[nt]:', nt_onehot_index[nt])
				onehot_seq = np.insert(onehot_seq, nt_onehot_index[nt], [0], axis=0)


		return onehot_seq



	def one_hot_encode_genomic_data(self, seq_file):

		print("One-hot encoding of fixed-length genomic windows...")

		# These will just be the indexes in the entire feature/gwrvis_index data frames 
		valid_bed_coords = []

		num_lines = sum(1 for line in open(seq_file))
		all_onehot_seqs = np.empty(shape=(num_lines, 4, self.win_len))

		seqs_with_n = {}

		tmp_win_id = -1
		with open(seq_file) as f:
			for seq in f.readlines():
			
				tmp_win_id += 1
			
				# discard sequences shorter than win_len
				if len(seq) != (self.win_len + 1):
					continue

				#if seq.count('N') > 10:
				if 'N' in seq:
					seqs_with_n[seq] = seqs_with_n.get(seq, 0) + 1
					continue
				else:
					valid_bed_coords.append(tmp_win_id)
					
				# A=[1, 0, 0, 0]
				# T=[0, 0, 0, 1]
				# G=[0, 0, 1, 0]
				# C=[0, 1, 0, 0]
				onehot_seq = self.onehot_encoded_seq(seq)

				try:
					all_onehot_seqs[tmp_win_id] = onehot_seq
				except:
					print('all_onehot_seqs:', all_onehot_seqs.shape)
					print('onehot_seq:', onehot_seq.shape)
					print(onehot_seq[:, :10])
					print(seq)
					sys.exit()
				# DEBUG
				#all_onehot_seqs[tmp_win_id] = None

		# Keep only sequences that didn't contain 'N's
		all_onehot_seqs = all_onehot_seqs[valid_bed_coords]
		print(all_onehot_seqs)
		print(all_onehot_seqs.shape)

		print('Total num. of unique seqs with Ns:', len(seqs_with_n.keys()))
		print('Total num. of variant entris to filter out:', sum(seqs_with_n.values()))

		return all_onehot_seqs, valid_bed_coords


		




	def filter_invalid_entries(self, subset_feature_df, valid_singlent_coords):

		print("\nFiltering invalid entries from 'subset_feature_df' ...")
		print('(Original) subset_feature_df:', subset_feature_df.shape)

		subset_feature_df = subset_feature_df.iloc[ valid_singlent_coords, :]
		
		print('(After filtering of Ns) subset_feature_df:', subset_feature_df.shape)
		
		
		return subset_feature_df







	def prepare_input_table_and_predict(self, subset_feature_df, cur_start_index, cur_end_index, pred_out_file):
		
		
		mut_rate_dict = get_mutability_rates(kmer=7)
		
		
		# Create seq_coords_file with windows centered around each variant
		# gwrvis_index_df will have just the row number of the original entries, and those having seqs with 'N' will be filtered out eventually
		seq_coords_file = self.seq_out_dir + '/' + self.patho_benign_sets + '.seq_coords.txt'
		seq_coords_df = subset_feature_df[['chr', 'start', 'end']].copy()

	


		seq_coords_df['chr'] = seq_coords_df['chr'].str.replace('chr', '').astype(int)
		seq_coords_df['start'] = seq_coords_df['start'] - int(self.win_len / 2)
		seq_coords_df['end'] = seq_coords_df['end'] + int(self.win_len / 2) - 1
		
		# Make sure genomic coordinates are non-negative - Discard those sequence (shorter than win_len) in the one_hot_encode_genomic_data() function
		seq_coords_df.loc[ seq_coords_df['start'] < 0, 'start'] = 0


		
		seq_coords_df.sort_values(by=['chr', 'start'], inplace=True)
		

		print(seq_coords_df.head())
		print('seq_coords_df:', seq_coords_df.shape)


		seq_list = seq_coords_df['chr'].map(str) + ':' + seq_coords_df['start'].astype(str) + '-' + seq_coords_df['end'].astype(str)

		
		# Write seq coords to file
		cur_file_identifier = 'chr' + str(self.chrom) + '_' + str(cur_start_index) + '_' + str(cur_end_index)
		
		chr_windows_dir = self.seq_out_dir + '/chr' + str(self.chrom)
		if not os.path.exists(chr_windows_dir):
			os.makedirs(chr_windows_dir)
			
		windows_file = chr_windows_dir + '/' + cur_file_identifier + '.windows.txt'
		print(windows_file)
		
		
		cur_file_identifier
		
		out_fh = open(windows_file, 'w')	
		out_fh.write("\n".join(seq_list))
		out_fh.close()	

		
		

		print("Extracting raw sequences from reference genome...")		
		raw_seq_out_file = chr_windows_dir + '/' + cur_file_identifier + '.raw_seqs.out' 

		awk_fasta_collapser = "awk '/^>/ {printf(\"\\n%s\\n\",$0);next; } { printf(\"%s\",$0);}  END {printf(\"\\n\");}' | tail -n +2"
		cmd = '../utils/twoBitToFa ' + self.human_ref_genome_2bit + ' -seqList=' + windows_file + " /dev/stdout | " + awk_fasta_collapser + " | grep -v '>' > " + raw_seq_out_file 
		print('\n' + cmd)

		p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout, stderr = p.communicate()
		print(str(stderr, "utf-8").rstrip())
		

		print(raw_seq_out_file)
		print(subset_feature_df.head())
		# NEW - Get updated features for mut_rate, gc_content, cpg and cpg_islands
		

		mut_rate_list = []
		cpg_list = []
		cpg_islands_list = []
		gc_content_list = []
			

		with open(raw_seq_out_file) as fh:
			for seq in fh:
				#print(seq)

				
				# 7-mer mutability rate
				try:
					center_index = int(self.win_len / 2)
					center_nt = seq[center_index]
				
					heptamer_seq = seq[center_index-3: center_index+4]
					
				except:
					print('\nException: center_index:', center_index, '; seq:', heptamer_seq, '; len:', len(seq))
					mut_rate = 'NA'
					
				try:
					mut_rate = mut_rate_dict[heptamer_seq]
				except:
					print('\nException: no mut_rate for current 7-mer (' + heptamer_seq + '); Indexes:', cur_start_index, '-', cur_end_index)
					mut_rate = 'NA'
					
					
				mut_rate_list.append(mut_rate)

				# CpG dinucleotides
				cpg = seq.count('CG')
				cpg_list.append(cpg)
			

				# CpG islands (>=2 CpG dinucleotides, e.g.: CGCG, CGCGCGCG, etc.)
				cpg_islands = len(re.findall(re.compile('CG(CG)+'), seq))
				cpg_islands_list.append(cpg_islands)

				# GC content (ratio over the entire window size)
				gc_content = (seq.count('G') + seq.count('C')) / self.win_len
				gc_content_list.append(gc_content)
	

		print('Indexes:', cur_start_index, '-', cur_end_index, '; subset_feature_df:', subset_feature_df.shape, '; mut_rate_list:', len(mut_rate_list))
	
		# Update full feature table and additional features files	
		subset_feature_df['mut_rate'] = mut_rate_list
		subset_feature_df['cpg'] = cpg_list
		subset_feature_df['cpg_islands'] = cpg_islands_list
		subset_feature_df['gc_content'] = gc_content_list
		print(subset_feature_df.head())

		


		# > One-hot encode raw genomic sequences and filter out sequences with 'N's
		filtered_onehot_seqs, valid_singlent_coords = self.one_hot_encode_genomic_data(raw_seq_out_file)
		print('valid_singlent_coords:', len(valid_singlent_coords))


		# Filter out entries that had sequences with Ns
		subset_feature_df = self.filter_invalid_entries(subset_feature_df, valid_singlent_coords)


		
		ordered_feature_cols = ['gwrvis', 'mut_rate', 'gc_content', 'cpg', 'cpg_islands', 'CTCF_Binding_Site', 'Enhancer', 'Open_chromatin', 'TF_binding_site', 'H3K27ac', 'H3K27me3', 'H4K20me1', 'H3K9ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K36me3']

		print(subset_feature_df.columns)
		print(subset_feature_df.shape)
		
		
		
		feature_df = subset_feature_df[ ordered_feature_cols ]
		print(feature_df.head())
		print(feature_df.columns)
		print(feature_df.shape)		


		# Get X, y, seqs arrays and run through pre-trained model		
		X = feature_df.values
		print('X:', X.shape)

		
		
		

		# 'seqs': reshape to (num_of_seqs, seq_length, 4)
		filtered_onehot_seqs = np.transpose(filtered_onehot_seqs, axes=(0,2,1))
		print('filtered_onehot_seqs:', filtered_onehot_seqs.shape)

		
		pred_scores = train_with_cv(X, input_features=self.input_features, seqs=filtered_onehot_seqs)
		
		
		# force garbage collector to free up memory
		X = None
		filtered_onehot_seqs = None
		
		

		prediction_df = subset_feature_df[['chr', 'start', 'end']].copy()
		del subset_feature_df
		
		prediction_df['jarvis'] = pred_scores
		
		print(prediction_df.head())
		print(prediction_df.tail())
		print(prediction_df.shape)
		
		
		# TODO: change output filename to reflect start-end indexes
		prediction_df.to_csv(pred_out_file, index=False, header=False, sep='\t')
		
		print('\nSaved jarvis predicted scores in:', pred_out_file)
		
		
		
		


	def compile_feature_table_incl_raw_seqs(self):
		
		self.input_table_name = 'full_feature_table.' + self.patho_benign_sets + '.bed.chr' + str(chrom) + '.single_nt_gwrvis.bed'
		
		full_feature_table_file = self.feature_tables_dir + '/' + self.input_table_name
		print('\nPredicting on\n', full_feature_table_file)
	
		# DEBUG
		#full_feature_table_file = self.feature_tables_dir + '/sample.chr21.single_nt.bed'


		print('\nReading full table for chromosome:', self.chrom)
		full_feature_df = pd.read_csv(full_feature_table_file, sep='\t', low_memory=False)
		print('full_feature_df:', full_feature_df.shape)
		print(full_feature_df.head())
		print(full_feature_df.shape)

		
		data_batch_size = 50000
		global_start_index = full_feature_df.index[0]
		global_end_index = full_feature_df.index[-1]


		# ======= TEMP - DEBUG =======
		#data_batch_size = 5
		#global_start_index = 50656165
		#global_end_index = 50656170
		

		
		print('global_start_index:', global_start_index)
		print('global_end_index:', global_end_index)
		
		
		
		# Multi-processing run across different start-end indexes
		process_jobs = []
		
		for cur_start_index in range(global_start_index, global_end_index, data_batch_size):
					
			cur_end_index = cur_start_index + data_batch_size
			
			#print('\nCurrent range:', cur_start_index, '-', cur_end_index)
			subset_feature_df = full_feature_df.iloc[cur_start_index : cur_end_index].copy()
			#print(subset_feature_df.head())
			#print(subset_feature_df.tail())
			#print(subset_feature_df.shape)

			pred_out_file = self.jarvis_predictions_per_chr_dir + '/chr' + self.chrom + '.' + self.input_features + '-features.' + str(cur_start_index) + '_' + str(cur_end_index) + '.jarvis'
			
			# DEBUG
			#self.prepare_input_table_and_predict(subset_feature_df, cur_start_index, cur_end_index, pred_out_file)


			if os.path.exists(pred_out_file):
				continue
			else:
				print('> Missing jarvis file:', pred_out_file)


			
			# Prepare input table (including raw sequences) and predict jarvis scores
			p = Process(target=self.prepare_input_table_and_predict, args=(subset_feature_df, cur_start_index, cur_end_index, pred_out_file))

			process_jobs.append(p)
			p.start()
			
			
			if len(process_jobs) >= self.NTHREADS:
				for p in process_jobs:
					p.join()
					process_jobs = []

			
		# Wait for any running process left
		for p in process_jobs:
			p.join()


	

if __name__ == '__main__':

	config_file = sys.argv[1]
	input_features = sys.argv[2]
	chrom = sys.argv[3]
	
	
	data_preprocessor = JarvisDataPreprocessing(config_file, input_features, chrom)

	
	# Extract raw sequences from input variant windows and combine with original feature set
	data_preprocessor.compile_feature_table_incl_raw_seqs()
	
