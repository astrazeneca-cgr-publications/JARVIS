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

	def __init__(self, config_file, predict_on_test_set=False, test_indexes=[]):

		print("Initialising new JarvisDataPreprocessing object...")

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
		gwrvis_scores_dir = out_dir + '/gwrvis_scores'

		self.ml_data_dir = out_dir + '/ml_data'
		if not os.path.exists(self.ml_data_dir):
			os.makedirs(self.ml_data_dir)
		self.seq_out_dir = self.ml_data_dir + '/raw_seq'
		if not os.path.exists(self.seq_out_dir):
			os.makedirs(self.seq_out_dir)
		self.feature_tables_dir = self.ml_data_dir + '/clinvar_feature_tables'
		if not os.path.exists(self.feature_tables_dir):
			os.makedirs(self.feature_tables_dir)


		# Specificy input (static) files
		self.human_ref_genome_2bit = '../' + self.hg_version +  '/homo_sapiens_GRCh' + self.grch[self.hg_version] + '_FASTA/hsa' + self.grch[self.hg_version] + '.2bit'
		self.all_gwrvis_bed_file = gwrvis_scores_dir + '/full_genome.all_gwrvis.bed'
		self.full_feature_table_file = self.feature_tables_dir + '/full_feature_table.' + self.patho_benign_sets + '.bed'



		# =============  Ad-hoc analysis: match TSS distance of control variants with pathogenic ones  =============
		"""
		self.full_feature_table_file = self.feature_tables_dir + '/full_feature_table.' + self.patho_benign_sets + '.matched.bed'
		"""




		self.predict_on_test_set = predict_on_test_set
		self.test_indexes = test_indexes
		if self.predict_on_test_set:

			# -- Reading full_gwrvis_and_regulatory_features.All_genomic_classes.tsv [Deprecated]
			# Pre-process full table: add 'clinvar_annot' column with placeholder 'Benign' values for consistent processing
			#tmp_full_feature_table_file = self.ml_data_dir + '/feature_tables/full_gwrvis_and_regulatory_features.All_genomic_classes.tsv'
			#tmp_df = pd.read_csv(tmp_full_feature_table_file, sep='\t')
			
			#tmp_df.insert(3, 'clinvar_annot', 'Benign')
			#tmp_df.drop( tmp_df.columns[len(tmp_df.columns)-1], axis=1, inplace=True)

			# -- NEW (includes clinvar_annot already + phastCons_primate)
			tmp_df = pd.read_csv(self.full_feature_table_file, sep='\t')

			print(tmp_df.head())
			print(tmp_df.columns)
			print(tmp_df.shape)

			
			tmp_df = tmp_df.iloc[ self.test_indexes[0]:self.test_indexes[1], : ]
			
			self.full_feature_table_file = self.ml_data_dir + '/feature_tables/full_gwrvis_and_regulatory_features.All_genomic_classes.tsv.for_prediction.' + str(self.test_indexes[0]) + '_' + str(self.test_indexes[1])
			tmp_df.to_csv(self.full_feature_table_file, sep='\t', header=True, index=False)
			
			print("Succesfully created new JarvisDataPreprocessing object!")



	def read_all_gwrvis_scores(self):

		all_gwrvis_bed_df = pd.read_csv(self.all_gwrvis_bed_file, header=0, sep='\t')
		print('all_gwrvis_bed_df:', all_gwrvis_bed_df.shape)

		# remove NAs
		all_gwrvis_bed_df.dropna(axis=0, inplace=True)
		print('all_gwrvis_bed_df (no NAs):',all_gwrvis_bed_df.shape)

		return all_gwrvis_bed_df




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
		valid_global_win_ids = []

		num_lines = sum(1 for line in open(seq_file))
		all_onehot_seqs = np.empty(shape=(num_lines, 4, self.win_len))

		seqs_with_n = {}

		tmp_win_id = -1
		with open(seq_file) as f:
			for seq in f.readlines():
			
				tmp_win_id += 1
				#if seq.count('N') > 10:
				if 'N' in seq:
					seqs_with_n[seq] = seqs_with_n.get(seq, 0) + 1
					continue
				else:
					valid_global_win_ids.append(tmp_win_id)
					
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
		all_onehot_seqs = all_onehot_seqs[valid_global_win_ids]
		print(all_onehot_seqs)
		print(all_onehot_seqs.shape)

		print('Total num. of unique seqs with Ns:', len(seqs_with_n.keys()))
		print('Total num. of variant entris to filter out:', sum(seqs_with_n.values()))

		return all_onehot_seqs, valid_global_win_ids


		

	def get_additional_feature_df(self, additional_features_file, additional_features_columns):

		print('\n', additional_features_file)
		additional_features_df = pd.read_csv(additional_features_file, sep='\t', header=None, index_col=None)
		print(additional_features_df.head())
		
		print(additional_features_columns)
		
		additional_features_columns = ['win_index'] + additional_features_columns
		additional_features_df.columns = additional_features_columns

		return additional_features_df




	# @anchor - use_gwrvis_windows: True or False
	def get_raw_seqs_from_variant_windows(self, use_gwrvis_windows=False):
		
		
		mut_rate_dict = get_mutability_rates(kmer=7)
		
		raw_seq_out_file = self.seq_out_dir + '/' + self.patho_benign_sets + '.raw_seqs.out' 


		gwrvis_seq_coords_file = self.seq_out_dir + '/' + self.patho_benign_sets + '.gwrvis_seq_coords.txt'



		os.system("tail -n+2 " + self.all_gwrvis_bed_file  + """ | awk '{print $2"\t"$3"\t"$4"\t"$1 }' """ + " | sed 's/^chr//' > " + self.all_gwrvis_bed_file+".tmp")
		print(self.all_gwrvis_bed_file+".tmp")
		os.system("tail -n+2 " + self.full_feature_table_file + " | sed 's/^chr//' > " + self.full_feature_table_file+".tmp")
		print(self.full_feature_table_file+".tmp")
		os.system("intersectBed -wo -a " + self.all_gwrvis_bed_file+".tmp" + " -b " + self.full_feature_table_file+".tmp" + " | cut -f1,2,3,4 > " + gwrvis_seq_coords_file)
		

		


		gwrvis_seq_coords_df = pd.read_csv(gwrvis_seq_coords_file, sep='\t', header=None)
		gwrvis_seq_coords_df.columns = ['chr', 'start', 'end', 'win_index']

		# record window indexes that correspond to each variant entry
		gwrvis_and_index_df = gwrvis_seq_coords_df[['chr', 'win_index']].copy()
		gwrvis_and_index_df.reset_index(drop=True, inplace=True)


		full_feature_df = pd.read_csv(self.full_feature_table_file, sep='\t')
		

		
		

		if not use_gwrvis_windows:
			# Create seq_coords_file with windows centered around each variant
			# gwrvis_index_df will have just the row number of the original entries, and those having seqs with 'N' will be filtered out eventually
			seq_coords_file = self.seq_out_dir + '/' + self.patho_benign_sets + '.seq_coords.txt'

			seq_coords_df = full_feature_df[['chr', 'start', 'end']].copy()
		


			seq_coords_df['chr'] = seq_coords_df['chr'].str.replace('chr', '').astype(int)
			seq_coords_df['start'] = seq_coords_df['start'] - int(self.win_len / 2)
			seq_coords_df['end'] = seq_coords_df['end'] + int(self.win_len / 2) - 1
			
			seq_coords_df.sort_values(by=['chr', 'start'], inplace=True)
			
	
			print(seq_coords_df.head())
			print('seq_coords_df:', seq_coords_df.shape)

		else:
			seq_coords_df = gwrvis_seq_coords_df.copy()
			seq_coords_df['end'] += 1 # include right-most part of the interval  -- redundant


		seq_list = seq_coords_df['chr'].map(str) + ':' + seq_coords_df['start'].astype(str) + '-' + seq_coords_df['end'].astype(str)

		

		windows_file = self.seq_out_dir + '/' + self.patho_benign_sets + '.windows.txt'
		out_fh = open(windows_file, 'w')	
		out_fh.write("\n".join(seq_list))
		out_fh.close()	

		print(windows_file)
		

		print("Extracting raw sequences from reference genome...")
		awk_fasta_collapser = "awk '/^>/ {printf(\"\\n%s\\n\",$0);next; } { printf(\"%s\",$0);}  END {printf(\"\\n\");}' | tail -n +2"
		cmd = '../utils/twoBitToFa ' + self.human_ref_genome_2bit + ' -seqList=' + windows_file + " /dev/stdout | " + awk_fasta_collapser + " | grep -v '>' > " + raw_seq_out_file 
		print('\n' + cmd)

		p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout, stderr = p.communicate()
		print(str(stderr, "utf-8").rstrip())
		
		
		print(full_feature_df.head())
		# NEW - Get updated features for mut_rate, gc_content, cpg and cpg_islands
		

		mut_rate_list = []
		cpg_list = []
		cpg_islands_list = []
		gc_content_list = []
			

		with open(raw_seq_out_file) as fh:
			for seq in fh:
				#print(seq)

				center_index = int(self.win_len / 2)

				center_nt = seq[center_index]
				
				# 7-mer mutability rate
				heptamer_seq = seq[center_index-3: center_index+4]
				mut_rate = mut_rate_dict[heptamer_seq]
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
	


		# Update full feature table and additional features files	
		full_feature_df['mut_rate'] = mut_rate_list
		full_feature_df['cpg'] = cpg_list
		full_feature_df['cpg_islands'] = cpg_islands_list
		full_feature_df['gc_content'] = gc_content_list
		full_feature_df.to_csv(self.full_feature_table_file, index=False, sep='\t')




		# Save into files
		os.system("tail -n+2 " + self.full_feature_table_file + " | sed 's/^chr//' > " + self.full_feature_table_file+".tmp")

		additional_features_file = self.seq_out_dir + '/' + self.patho_benign_sets + '.additional_features.tsv'
		os.system("intersectBed -wo -a " + self.all_gwrvis_bed_file+".tmp" + " -b " + self.full_feature_table_file+".tmp" + " | cut --complement -f1,2,3 | rev | cut --complement -f1 | rev > " + additional_features_file)




		return raw_seq_out_file, gwrvis_and_index_df, additional_features_file


		

	def filter_invalid_entries(self, gwrvis_and_index_df, additional_features_df, valid_global_win_ids):

		print("Filtering invalid entries from 'gwrvis_and_index_df' and 'additional_features_df' ...")
		print('(original) gwrvis_and_index_df:', gwrvis_and_index_df.shape)
		print(gwrvis_and_index_df.head())
		print(gwrvis_and_index_df.tail())

		print('(original) additional_features_df:', additional_features_df.shape)
		print(additional_features_df.head())
		print(additional_features_df.tail())


		gwrvis_and_index_df = gwrvis_and_index_df.iloc[ valid_global_win_ids, :]
		additional_features_df = additional_features_df.iloc[ valid_global_win_ids, :]
		
		
		if self.predict_on_test_set:
			# Save full feature table with valid global win indexes only
		
			tmp_df = pd.read_csv(self.full_feature_table_file, sep='\t')
			tmp_df = tmp_df.iloc[valid_global_win_ids, :]
		
		
			valid_full_feature_table_file = self.ml_data_dir + '/feature_tables/full_gwrvis_and_regulatory_features.All_genomic_classes.tsv.valid_windows.' + str(self.test_indexes[0]) + '_' + str(self.test_indexes[1])
			tmp_df.to_csv(valid_full_feature_table_file, sep='\t', header=True, index=False)
		
		
		
		
		return gwrvis_and_index_df, additional_features_df




	def compile_feature_table_incl_raw_seqs(self):

		# Read feature table with variant positions, gwRVIS, GC-related features and external annotations
		full_feature_table = pd.read_csv(self.full_feature_table_file, sep='\t', low_memory=False)
		full_feature_table.dropna(inplace=True)
		full_feature_table.to_csv(self.full_feature_table_file, sep='\t', header=True, index=False)

		print('full_feature_table:', full_feature_table.shape)
		print(full_feature_table.head())


		# Read all gwRVIS scores with BED-style genomic coordinates into a data frame
		all_gwrvis_bed_df = self.read_all_gwrvis_scores()	
		print(all_gwrvis_bed_df.head())
		print(all_gwrvis_bed_df.tail())
		print(max(all_gwrvis_bed_df['win_index']))
		print(len(all_gwrvis_bed_df['win_index']))
		print(len(all_gwrvis_bed_df['win_index'].unique()))


		# > Get raw sequences and windows indexes of most intolerant/tolerant windows
		seq_out_file, gwrvis_and_index_df, additional_features_file = self.get_raw_seqs_from_variant_windows()
		# The pair (chr, win_index) define uniquely each genomic window
		print('gwrvis_and_index_df:', gwrvis_and_index_df.shape)
		print(gwrvis_and_index_df.head())


		additional_features_df = self.get_additional_feature_df(additional_features_file, list(full_feature_table.columns.values))
		print('additional_features_df:', additional_features_df.shape)
		print(additional_features_df.head())
		

		# > One-hot encode raw genomic sequences and filter out sequences with 'N's
		filtered_onehot_seqs, valid_global_win_ids = self.one_hot_encode_genomic_data(seq_out_file)
		print('valid_global_win_ids:', len(valid_global_win_ids))



		# > Filter out windows with 'N's in their sequence
		gwrvis_and_index_df, additional_features_df = self.filter_invalid_entries(gwrvis_and_index_df, additional_features_df, valid_global_win_ids)
		print(gwrvis_and_index_df.head())
		print(additional_features_df.head())
		
		print('gwrvis_and_index_df:', gwrvis_and_index_df.shape)
		print('additional_features_df:', additional_features_df.shape)
		print('filtered_onehot_seqs:', filtered_onehot_seqs.shape)


		
		if self.Y_label == 'conservation_annot':
			# - For conservation labels
			additional_features_df[self.Y_label].replace(regex={r'.*Conserved.*': 1, r'.*Non_conserved.*': 0}, inplace=True)
		else:
			# - For pathogenicity labels
			additional_features_df[self.Y_label].replace(regex={r'.*Pathogenic.*': 1, r'.*Benign.*': 0}, inplace=True)
		print(additional_features_df.head())
		print(additional_features_df.tail())

		return additional_features_df, filtered_onehot_seqs




	def transform_and_save_data(self, additional_features_df, filtered_onehot_seqs):

		additional_features_df.reset_index(drop=True, inplace=True)
		additional_features_df['global_index'] = additional_features_df.index.values


		# Impute NA data with medians
		cols_to_impute = ['phastCons_primate']
		for col in cols_to_impute:
			tmp_median = additional_features_df.loc[ additional_features_df[col] != '.', col].astype(float).tolist()
			additional_features_df[col] = additional_features_df[col].apply(lambda x: np.median(tmp_median)  if x == '.' else float(x))

		
		print(additional_features_df.head())
		print(additional_features_df.tail())
		print(additional_features_df.info())
		print(filtered_onehot_seqs.shape)
		additional_features_df.to_csv('prepared_data.ncER-mendelian.tsv', sep='\t', index=False)
		sys.exit()



		
		# drop non-informative classes
		cols_to_drop = ['win_index', 'chr', 'start', 'end']
		coords_df = additional_features_df[['chr', 'start']].values
		print('coords_df:', coords_df)
		additional_features_df.drop(cols_to_drop, inplace=True, axis=1)
		
		# Get genomic classes
		genomic_classes = additional_features_df['genomic_class'].values
		additional_features_df.drop(['genomic_class'], inplace=True, axis=1)
		print('genomic_classes:', genomic_classes)

		# Get VCF-based features
		vcf_based_features = ['common_variants', 'common_vs_all_variants_ratio', 'all_variants', 'mean_ac', 'mean_af', 'bin_1', 'bin_2', 'bin_3', 'bin_4', 'bin_5', 'bin_6']
		vcf_features_df = additional_features_df[vcf_based_features].values
		additional_features_df.drop(vcf_based_features, inplace=True, axis=1)
		print('vcf_features_df:', vcf_features_df)

		# Get global index (to be used in other module for CV train/test splits)
		global_index = additional_features_df['global_index'].values
		additional_features_df.drop(['global_index'], inplace=True, axis=1)
		print('global_index:', global_index)

		# Get y label
		y = additional_features_df[self.Y_label].values
		y = np.reshape(y, (y.shape[0], 1))
		y = tf.keras.utils.to_categorical(y)
		additional_features_df.drop([self.Y_label], inplace=True, axis=1)
		print('y:', y)




		# Get main features (X)
		X = additional_features_df.copy().values
		print('X:', X)


		# [Deprecated] - Normalise features
		#X = MinMaxScaler().fit_transform(X)



		# 'seqs': reshape to (num_of_seqs, seq_length, 4)
		filtered_onehot_seqs = np.transpose(filtered_onehot_seqs, axes=(0,2,1))
		print('filtered_onehot_seqs:', filtered_onehot_seqs)


		# Integrate data into a dictionary
		data_dict = {'X': X, 'y': y, 'seqs': filtered_onehot_seqs, 
			     'vcf_features': vcf_features_df, 'genomic_classes': genomic_classes, 'global_index': global_index,
			     'X_cols': additional_features_df.columns.values, 'vcf_features_cols': vcf_based_features, 'coords_df': coords_df }
		

		out_file = self.save_data_to_files(data_dict)

		return out_file



			

	def save_data_to_files(self, data_dict):

		if not self.predict_on_test_set:
			out_file = self.ml_data_dir + '/jarvis_data.pkl'
		else:
			out_file = self.ml_data_dir + '/jarvis_data.' + str(self.test_indexes[0]) + '_' + str(self.test_indexes[1]) + '.pkl'
		
		# Save data dict to a pickle file
		pkl_out = open(out_file, 'wb')
		pickle.dump(data_dict, pkl_out, protocol=4)
		pkl_out.close()
		print('Saved all data to ' + out_file)

		return out_file
		
		
		

if __name__ == '__main__':

	config_file = sys.argv[1]
	run_params = custom_utils.get_config_params(config_file)

	data_preprocessor = JarvisDataPreprocessing(config_file)

	
	# Extract raw sequences from input variant windows and combine with original feature set
	additional_features_df, filtered_onehot_seqs = data_preprocessor.compile_feature_table_incl_raw_seqs()
	
	# Merge data, transform into form appropriate for DNN training and save into file
	data_dict_file = data_preprocessor.transform_and_save_data(additional_features_df, filtered_onehot_seqs)
	print(data_dict_file)
