import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import subprocess
from sklearn.preprocessing import LabelEncoder 
from sklearn.preprocessing import OneHotEncoder
from sklearn.model_selection import train_test_split
import tensorflow as tf
import pickle
from scipy.io import savemat
import sys, os

sys.path.insert(1, os.path.join(sys.path[0], '../..'))
import custom_utils


def read_all_gwrvis_scores(all_gwrvis_bed_file):

	all_gwrvis_bed_df = pd.read_csv(all_gwrvis_bed_file, header=0, sep='\t')
	print(all_gwrvis_bed_df.head())
	print(all_gwrvis_bed_df.shape)

	# remove NAs
	all_gwrvis_bed_df.dropna(axis=0, inplace=True)
	print(all_gwrvis_bed_df.shape)

	return all_gwrvis_bed_df




def onehot_encoded_seq(seq):

	# remove any trailing new line characters
	seq = seq.rstrip()
	seq_array = np.array(list(seq))
	

	label_encoder = LabelEncoder()
	integer_encoded_seq = label_encoder.fit_transform(seq_array)
	integer_encoded_seq = integer_encoded_seq.reshape(len(integer_encoded_seq), 1)

	onehot_encoder = OneHotEncoder(sparse=False, categories='auto')
	onehot_encoded_seq = onehot_encoder.fit_transform(integer_encoded_seq) #.tolist()
	onehot_encoded_seq = np.transpose(onehot_encoded_seq)

	return onehot_encoded_seq



def one_hot_encode_genomic_data(seq_file):

	print("One-hot encoding of fixed-length genomic windows...")

	valid_win_ids = []

	num_lines = sum(1 for line in open(seq_file))
	all_onehot_seqs = np.empty(shape=(num_lines, 4, win_len))

	seqs_with_n = {}

	with open(seq_file) as f:
		tmp_win_id = -1
		for seq in f.readlines():
			tmp_win_id += 1

			#if seq.count('N') > 10:
			if 'N' in seq:
				seqs_with_n[seq] = seqs_with_n.get(seq, 0) + 1
				continue
			else:
				valid_win_ids.append(tmp_win_id)
				
			# A = [1, 0, 0, 0]
			# T = [0, 0, 0, 1]
			# G = [0, 0, 1, 0]
			# C = [0, 1, 0, 0]
			onehot_seq = onehot_encoded_seq(seq)
			all_onehot_seqs[tmp_win_id] = onehot_seq

	# Keep only sequences that didn't contain 'N's
	all_onehot_seqs = all_onehot_seqs[valid_win_ids]
	print(all_onehot_seqs)
	print(all_onehot_seqs.shape)

	print('Total num. of unique seqs with Ns:', len(seqs_with_n.keys()))
	print('Total num. of variant entris to filter out:', sum(seqs_with_n.values()))

	return all_onehot_seqs, valid_win_ids


		
def filter_gwrvis_index_df(gwrvis_and_index_df, valid_win_ids):

	print(gwrvis_and_index_df.shape)
	gwrvis_and_index_df = gwrvis_and_index_df.loc[ gwrvis_and_index_df.win_index.isin(valid_win_ids), :]
	print(gwrvis_and_index_df.shape)
	
	return gwrvis_and_index_df



def integrate_additional_data(gwrvis_and_index_df):

	# Aggregate additional data
	agg_additional_data_df = pd.DataFrame()

	for chrom in range(1,23):
		cur_addit_data_file = additional_data_dir + '/chr' + str(chrom) + '.additional_features_df.csv'
		cur_addit_data_df = pd.read_csv(cur_addit_data_file, index_col=0)
		cur_addit_data_df['chr'] = 'chr' + str(chrom)
		cur_addit_data_df.reset_index(inplace=True)
		cur_addit_data_df.columns.values[0] = 'win_index'
		cur_addit_data_df['win_index'] = cur_addit_data_df['win_index'].astype(int)

		#print(cur_addit_data_df.head())
		#print(gwrvis_and_index_df.head())
		if agg_additional_data_df.shape[0] > 0:
			agg_additional_data_df = pd.concat([agg_additional_data_df, cur_addit_data_df], axis=0)
			print(agg_additional_data_df.shape)
		else:
			agg_additional_data_df = cur_addit_data_df
			print(agg_additional_data_df.shape)


	merged_features_df = gwrvis_and_index_df.merge(agg_additional_data_df, on=['chr', 'win_index'])
	print(merged_features_df.head())
	print(merged_features_df.tail())
	print(merged_features_df.shape)


	return merged_features_df
	





def get_raw_seqs_from_variant_windows(all_gwrvis_bed_file, full_feature_table_file, use_gwrvis_windows=True):
	
	raw_seq_out_file = seq_out_dir + '/' + patho_benign_sets + '.raw_seqs.out' 
	windows_file = seq_out_dir + '/' + patho_benign_sets + '.windows.txt'
	out_fh = open(windows_file, 'w')	


	if use_gwrvis_windows:
		seq_coords_file = seq_out_dir + '/' + patho_benign_sets + '.seq_coords.txt'
		os.system("tail -n+2 " + all_gwrvis_bed_file  + """ | awk '{print $2"\t"$3"\t"$4"\t"$1 }' """ + " | sed 's/^chr//' > " + all_gwrvis_bed_file+".tmp")
		os.system("tail -n+2 " + full_feature_table_file + " | sed 's/^chr//' > " + full_feature_table_file+".tmp")

		os.system("intersectBed -wo -a " + all_gwrvis_bed_file+".tmp" + " -b " + full_feature_table_file+".tmp" + " | cut -f1,2,3,4 > " + seq_coords_file)
	else:
		# Create seq_coords_file with windows centered around each variant
		pass
	print(seq_coords_file)


	seq_coords_df = pd.read_csv(seq_coords_file, sep='\t', header=None)
	seq_coords_df.columns = ['chr', 'start', 'end', 'win_index']
	print(seq_coords_df.head())

	# record window indexes that correspond to each variant entry
	gwrvis_and_index_df = seq_coords_df[['win_index']].copy()
	gwrvis_and_index_df.reset_index(drop=True, inplace=True)
	print(gwrvis_and_index_df.head())


	seq_coords_df['end'] += 1 # include right-most part of the interval
	print(seq_coords_df.head())
	seq_list = seq_coords_df['chr'].map(str) + ':' + seq_coords_df['start'].astype(str) + '-' + seq_coords_df['end'].astype(str)

	out_fh.write("\n".join(seq_list))
	out_fh.close()	

	# Extract raw sequences from reference genome
	awk_fasta_collapser = "awk '/^>/ {printf(\"\\n%s\\n\",$0);next; } { printf(\"%s\",$0);}  END {printf(\"\\n\");}' | tail -n +2"
	cmd = '../utils/twoBitToFa ' + human_ref_genome_2bit + ' -seqList=' + windows_file + " /dev/stdout | " + awk_fasta_collapser + " | grep -v '>' > " + raw_seq_out_file 
	print('\n' + cmd)

	p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = p.communicate()
	print(str(stderr, "utf-8").rstrip())

	return raw_seq_out_file, gwrvis_and_index_df


	


def compile_feature_table_incl_raw_seqs(all_gwrvis_bed_file, full_feature_table_file):

	# Read feature table with variant positions, gwRVIS, GC-related features and external annotations
	full_feature_table = pd.read_csv(full_feature_table_file, sep='\t', low_memory=False)
	print(full_feature_table.head())

	# Read all gwRVIS scores with BED-style genomic coordinates into a data frame
	all_gwrvis_bed_df = read_all_gwrvis_scores(all_gwrvis_bed_file)	
	print(all_gwrvis_bed_df.head())


	# > Get raw sequences and windows indexes of most intolerant/tolerant windows
	seq_out_file, gwrvis_and_index_df = get_raw_seqs_from_variant_windows(all_gwrvis_bed_file, full_feature_table_file)



	# > One-hot encode raw genomic sequences and filter out sequences with 'N's
	filtered_onehot_seqs, valid_win_ids = one_hot_encode_genomic_data(seq_out_file)


	# > Filter out windows with 'N's in their sequence
	gwrvis_and_index_df = filter_gwrvis_index_df(gwrvis_and_index_df, valid_win_ids)
	sys.exit()	

	# Add additional features in the compiled table	
	merged_features_df = integrate_additional_data(gwrvis_and_index_df)

	if tol_type == 'intolerant':
		merged_features_df['y'] = 1
	else:
		merged_features_df['y'] = 0

	return all_merged_df, variant_win_onehot_seqs





def check_gwrvis_extremes_distribution(all_merged_df):
	
	intol_gwrvis = all_merged_df.loc[ all_merged_df.y ==1, 'gwrvis'].values
	tol_gwrvis = all_merged_df.loc[ all_merged_df.y == 0, 'gwrvis'].values

	intol_median = np.median(intol_gwrvis)
	tol_median = np.median(tol_gwrvis)


	fig, ax = plt.subplots(figsize=(15, 15))
	p1 = sns.kdeplot(intol_gwrvis, shade=True, color='r', label='intolerant')
	p2 = sns.kdeplot(tol_gwrvis, shade=True, color='b', label='tolerant')

	plt.axvline(x=intol_median, linestyle='--', linewidth=1.5, color='r', label='median intol.')
	plt.axvline(x=tol_median, linestyle='--', linewidth=1.5, color='b', label='median tol.')
	
	plt.title('gwRVIS most intolerant/tolerant distribution')
	plt.xlabel('gwRVIS')
	plt.ylabel('Density')
	plt.legend(fontsize=18, markerscale=2)

	fig.savefig(out_dir + '/gwrvis_extremes_distribution.pdf', bbox_inches='tight')



def add_regulatory_features(all_merged_df):

	print(all_merged_df.head())
	



def split_data_into_train_val_test_sets(all_merged_df, all_filtered_onehot_seqs, test_size=0.2):

	validation_size = int(all_merged_df.shape[0] * test_size * 0.4)

	print(all_merged_df.head())
	print(all_merged_df.tail())
	print(all_filtered_onehot_seqs.shape)
	
	all_merged_df['global_index'] = all_merged_df.index.values

	# Get a stratified random validation set for additional features and raw sequences
	validation_set_positive = all_merged_df.loc[ all_merged_df.y == 1, :].sample(n=int(validation_size/2), replace=False)
	validation_set_negative = all_merged_df.loc[ all_merged_df.y == 0, :].sample(n=int(validation_size/2), replace=False)
	validation_set = pd.concat([validation_set_positive, validation_set_negative], axis=0)
	validation_indexes = validation_set.global_index.values
	X_validation = validation_set.loc[:, validation_set.columns != 'y']
	y_validation = validation_set['y']

	validation_seq_set = all_filtered_onehot_seqs[validation_indexes]
	

	# Remove validation set from the whole dataset
	all_merged_df.drop(validation_indexes, axis=0, inplace=True)
	#print(all_merged_df.head())
	X = all_merged_df.loc[:, all_merged_df.columns != 'y']
	y = all_merged_df['y']


	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size)
	train_indexes = X_train.index.values
	test_indexes = X_test.index.values

	train_seq_set = all_filtered_onehot_seqs[train_indexes]
	test_seq_set = all_filtered_onehot_seqs[test_indexes]


	# Integrate data into dictionaries
	train_dict = {'X': X_train, 'y': y_train, 'seqs': train_seq_set}
	validation_dict = {'X': X_validation, 'y': y_validation, 'seqs': validation_seq_set}
	test_dict = {'X': X_test, 'y': y_test, 'seqs': test_seq_set}
	
	return train_dict, validation_dict, test_dict
	
	

def transform_data(train_dict, validation_dict, test_dict):

	# Store 'gwrvis' in separate field in the dictionary
	train_gwrvis = train_dict['X']['gwrvis'].values
	train_gwrvis = np.reshape(train_gwrvis, (train_gwrvis.shape[0], 1))
	train_dict['gwrvis'] = train_gwrvis

	test_gwrvis = test_dict['X']['gwrvis'].values
	test_gwrvis = np.reshape(test_gwrvis, (test_gwrvis.shape[0], 1))
	test_dict['gwrvis'] = test_gwrvis

	validation_gwrvis = validation_dict['X']['gwrvis'].values
	validation_gwrvis = np.reshape(validation_gwrvis, (validation_gwrvis.shape[0], 1))
	validation_dict['gwrvis'] = validation_gwrvis

	# Store other/additional features in a separate field
	other_features = ['mut_rate', 'gc_content', 'cpg', 'cpg_islands']
	cols_to_drop = ['gwrvis'] + other_features

	train_dict['X_other'] = train_dict['X'][other_features].values
	train_dict['X'].drop(cols_to_drop, axis=1, inplace=True)

	test_dict['X_other'] = test_dict['X'][other_features].values
	test_dict['X'].drop(cols_to_drop, axis=1, inplace=True)

	validation_dict['X_other'] = validation_dict['X'][other_features].values
	validation_dict['X'].drop(cols_to_drop, axis=1, inplace=True)


	# 'seqs': reshape to (num_of_seqs, seq_length, 4)
	train_dict['seqs'] = np.transpose(train_dict['seqs'], axes=(0,2,1))
	test_dict['seqs'] = np.transpose(test_dict['seqs'], axes=(0,2,1))
	validation_dict['seqs'] = np.transpose(validation_dict['seqs'], axes=(0,2,1))


	# 'y': reshape and convert to_categorical
	train_y = train_dict['y'].values
	train_y = np.reshape(train_y, (train_y.shape[0], 1))
	train_dict['y'] = tf.keras.utils.to_categorical(train_y)

	validation_y = validation_dict['y'].values
	validation_y = np.reshape(validation_y, (validation_y.shape[0], 1))
	validation_dict['y'] = tf.keras.utils.to_categorical(validation_y)

	test_y = test_dict['y'].values
	test_y = np.reshape(test_y, (test_y.shape[0], 1))
	test_dict['y'] = tf.keras.utils.to_categorical(test_y)	

	print('\n> X_train:', train_dict['X'].shape)
	print('> y_train:', train_dict['y'].shape)
	print('> train_seqs:', train_dict['seqs'].shape)
	
	print('\n> X_test:', test_dict['X'].shape)
	print('> y_test:', test_dict['y'].shape)
	print('> test_seqs:', test_dict['seqs'].shape)

	print('\n> X_validation:', validation_dict['X'].shape)
	print('> y_validation:', validation_dict['y'].shape)
	print('> validation_seqs:', validation_dict['seqs'].shape)

	print(test_dict['X'])
	
	return train_dict, validation_dict, test_dict

		

def save_data_to_files(train_dict, validation_dict, test_dict, random_seqs=False):

	top_ratio_str = '.top_' + str(top_ratio)

	if random_seqs:
		top_ratio_str += '.random'
	
	# Save data to files
	pkl_out = open(ml_data_dir + '/train' + top_ratio_str + '.pkl', 'wb')
	pickle.dump(train_dict, pkl_out, protocol=4)
	pkl_out.close()
	
	pkl_out = open(ml_data_dir + '/validation' + top_ratio_str + '.pkl', 'wb')
	pickle.dump(validation_dict, pkl_out, protocol=4)
	pkl_out.close()

	pkl_out = open(ml_data_dir + '/test' + top_ratio_str + '.pkl', 'wb')
	pickle.dump(test_dict, pkl_out, protocol=4)
	pkl_out.close()

	
	
	

if __name__ == '__main__':

	config_file = sys.argv[1]


	# ==== Read config parameters ====
	config_params = custom_utils.get_config_params(config_file)
	win_len = config_params['win_len']
	pathogenic_set = config_params['pathogenic_set']
	benign_set = config_params['benign_set']
	patho_benign_sets = pathogenic_set + '_' + benign_set


	# ==== Define dir structure ====
	out_dir = custom_utils.create_out_dir(config_file)
	additional_data_dir = out_dir + '/tmp'

	gwrvis_scores_dir = out_dir + '/gwrvis_scores'
	ml_data_dir = out_dir + '/ml_data'
	if not os.path.exists(ml_data_dir):
		os.makedirs(ml_data_dir)
	seq_out_dir = ml_data_dir + '/raw_seq'
	if not os.path.exists(seq_out_dir):
		os.makedirs(seq_out_dir)
	feature_tables_dir = ml_data_dir + '/clinvar_feature_tables'
	if not os.path.exists(feature_tables_dir):
		os.makedirs(feature_tables_dir)


	# Specificy input (static) files
	human_ref_genome_2bit = '../hg19/homo_sapiens_GRCh37_FASTA/hsa37.2bit'
	all_gwrvis_bed_file = gwrvis_scores_dir + '/full_genome.all_gwrvis.bed'
	full_feature_table_file = feature_tables_dir + '/full_feature_table.' + patho_benign_sets + '.bed'
	print(all_gwrvis_bed_file)
	print(full_feature_table_file)



	
	# Extract raw sequences from input variant windows and combine with original feature set
	all_merged_df, variant_win_onehot_seqs = compile_feature_table_incl_raw_seqs(all_gwrvis_bed_file, full_feature_table_file)
	sys.exit()
	
	
	print('===========================================')


	all_merged_df = pd.concat([intol_merged_df, tol_merged_df], axis=0)
	all_merged_df.reset_index(inplace=True, drop=True)
	all_merged_df.drop(['tolerance_rank'], axis=1, inplace=True)
	print(all_merged_df.shape)
	
	check_gwrvis_extremes_distribution(all_merged_df)


	#all_merged_df = add_regulatory_features(all_merged_df)
	

	all_filtered_onehot_seqs = np.concatenate((intol_filtered_onehot_seqs, tol_filtered_onehot_seqs), axis=0)
	print(all_filtered_onehot_seqs.shape)


	train_dict, validation_dict, test_dict = split_data_into_train_val_test_sets(all_merged_df, all_filtered_onehot_seqs)

	
	train_dict, validation_dict, test_dict = transform_data(train_dict, validation_dict, test_dict)
	print(test_dict)
	
	save_data_to_files(train_dict, validation_dict, test_dict, random_seqs=random_seqs)
