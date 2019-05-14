import pandas as pd
import numpy as np
import sys, os
import subprocess
from sklearn.preprocessing import LabelEncoder 
from sklearn.preprocessing import OneHotEncoder
from sklearn.model_selection import train_test_split
sys.path.insert(1, os.path.join(sys.path[0], '..'))
import custom_utils


def read_all_gwrvis_scores(all_gwrvis_bed_file):

	all_gwrvis_bed_df = pd.read_csv(all_gwrvis_bed_file, header=0, sep='\t')
	print(all_gwrvis_bed_df.head())
	print(all_gwrvis_bed_df.shape)

	# remove NAs
	all_gwrvis_bed_df.dropna(axis=0, inplace=True)
	print(all_gwrvis_bed_df.shape)

	return all_gwrvis_bed_df



def get_most_intol_or_toler_sequences(all_gwrvis_bed_df, tol_type='intolerant', top_ratio=0.01):
	
	# Sort by gwRVIS value in ascending (for intolerant) or descending (for tolerant) order
	ascending = True
	if tol_type == 'tolerant':
		ascending = False
	all_gwrvis_bed_df.sort_values(by='gwrvis', inplace=True, ascending=ascending)
	all_gwrvis_bed_df.reset_index(drop=True, inplace=True)
	print(all_gwrvis_bed_df.head())
	print(all_gwrvis_bed_df.tail())


	top_seqs_set_size = int(all_gwrvis_bed_df.shape[0] * top_ratio)
	print(top_seqs_set_size)


	raw_seq_out_file = seq_out_dir + '/most_' + tol_type + '.raw_seqs.out' 
	windows_file = seq_out_dir + '/most_' + tol_type + '.windows.txt'
	out_fh = open(windows_file, 'w')	

	
	seq_list_df = all_gwrvis_bed_df[:top_seqs_set_size]
	gwrvis_and_index_df = seq_list_df.copy()
	gwrvis_and_index_df.reset_index(inplace=True)	
	gwrvis_and_index_df.columns.values[0] = 'tolerance_rank'
	#print(gwrvis_and_index_df.head())


	seq_list_df = seq_list_df[['chr', 'start', 'end']]
	seq_list_df['chr'] = seq_list_df['chr'].str.replace('chr', '')
	seq_list_df['end'] += 1 # include right-most part of the interval
	seq_list = seq_list_df['chr'].map(str) + ':' + seq_list_df['start'].astype(str) + '-' + seq_list_df['end'].astype(str)

	out_fh.write("\n".join(seq_list))
	out_fh.close()	

	# Extract raw sequences from reference genome
	awk_fasta_collapser = "awk '/^>/ {printf(\"\\n%s\\n\",$0);next; } { printf(\"%s\",$0);}  END {printf(\"\\n\");}' | tail -n +2"
	cmd = '../util/twoBitToFa ' + human_ref_genome_2bit + ' -seqList=' + windows_file + " /dev/stdout | " + awk_fasta_collapser + " | grep -v '>' > " + raw_seq_out_file # | grep -v '>' | tr '\n' ' ' | sed 's/ //g'"
	print('\n' + cmd)

	p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = p.communicate()
	#print(str(stderr, "utf-8").rstrip())

	return raw_seq_out_file, gwrvis_and_index_df



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

	valid_rank_ids = []

	num_lines = sum(1 for line in open(seq_file))
	all_onehot_seqs = np.empty(shape=(num_lines, 4, win_len))

	with open(seq_file) as f:
		rank_id = -1
		for seq in f.readlines():
			rank_id += 1

			if 'N' in seq:
				continue
			else:
				valid_rank_ids.append(rank_id)
				
			onehot_seq = onehot_encoded_seq(seq)
			all_onehot_seqs[rank_id] = onehot_seq

	# Keep only sequences that didn't contain 'N's
	all_onehot_seqs = all_onehot_seqs[valid_rank_ids]
	print(all_onehot_seqs)
	print(all_onehot_seqs.shape)

	return all_onehot_seqs, valid_rank_ids


		
def filter_gwrvis_index_df(gwrvis_and_index_df, valid_rank_ids):

	print(gwrvis_and_index_df.shape)
	gwrvis_and_index_df = gwrvis_and_index_df.loc[ gwrvis_and_index_df.tolerance_rank.isin(valid_rank_ids), :]
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
	
	

def compile_feature_table_per_tolerance_class(all_gwrvis_bed_df, tol_type='intolerant'):

	# Get raw sequences and windows indexes of most intolerant/tolerant windows
	seq_out_file, gwrvis_and_index_df = get_most_intol_or_toler_sequences(all_gwrvis_bed_df, tol_type=tol_type)

	# One-hot encode raw genomic sequences and filter out sequences with 'N's
	filtered_onehot_seqs, valid_rank_ids = one_hot_encode_genomic_data(seq_out_file)

	# Filter out windows with 'N's in their sequence
	gwrvis_and_index_df = filter_gwrvis_index_df(gwrvis_and_index_df, valid_rank_ids)
	
	# Add additional features in the compiled table	
	merged_features_df = integrate_additional_data(gwrvis_and_index_df)

	if tol_type == 'intolerant':
		merged_features_df['y'] = 1
	else:
		merged_features_df['y'] = 0

	return merged_features_df, filtered_onehot_seqs



def split_data_into_train_val_test_sets(all_merged_df, all_filtered_onehot_seqs, validation_size=1000, test_size=0.2):

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
	print(X_validation.head())
	print(y_validation.head())

	validation_seq_set = all_filtered_onehot_seqs[validation_indexes]
	print(validation_seq_set.shape)
	print(validation_seq_set)
	

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

	print('\n> X_train:', X_train.shape)
	print(X_train.head())
	print('\n> y_train:', y_train.shape)
	print(y_train.head())
	print(train_seq_set.shape)

	print('\n> X_test:', X_test.shape)
	print(X_test.head())
	print('\n> y_test:', y_test.shape)
	print(y_test.head())
	print(test_seq_set.shape)


	# Save data to files
	X_validation.to_csv(ml_data_dir + '/X_validation.tsv', sep='\t')
	y_validation.to_csv(ml_data_dir + '/y_validation.tsv', sep='\t')
	np.savetxt(ml_data_dir + "/validation_seq_set.txt", validation_seq_set)

	X_train.to_csv(ml_data_dir + '/X_train.tsv', sep='\t')
	y_train.to_csv(ml_data_dir + '/y_train.tsv', sep='\t')
	np.savetxt(ml_data_dir + "/train_seq_set.txt", train_seq_set)
	
	X_test.to_csv(ml_data_dir + '/X_test.tsv', sep='\t')
	y_test.to_csv(ml_data_dir + '/y_test.tsv', sep='\t')
	np.savetxt(ml_data_dir + "/test_seq_set.txt", test_seq_set)




if __name__ == '__main__':

	config_file = sys.argv[1]
	config_params = custom_utils.get_config_params(config_file)
	win_len = config_params['win_len']

	human_ref_genome_2bit = '../../hg19/homo_sapiens_GRCh37_FASTA/hsa37.2bit'

	out_dir = custom_utils.create_out_dir(config_file)
	out_dir = '../' + out_dir
	additional_data_dir = out_dir + '/tmp'

	gwrvis_scores_dir = out_dir + '/gwrvis_scores'
	ml_data_dir = out_dir + '/ml_data'
	if not os.path.exists(ml_data_dir):
		os.makedirs(ml_data_dir)
	seq_out_dir = ml_data_dir + '/raw_seq'
	if not os.path.exists(seq_out_dir):
		os.makedirs(seq_out_dir)

	all_gwrvis_bed_file = gwrvis_scores_dir + '/full_genome.all_gwrvis.bed'

	# Read all gwRVIS scores with BED-style genomic coordinates into a data frame
	all_gwrvis_bed_df = read_all_gwrvis_scores(all_gwrvis_bed_file)	

	
	# Get features (including raw sequences) for most *intolerant* windows
	intol_merged_df, intol_filtered_onehot_seqs = compile_feature_table_per_tolerance_class(all_gwrvis_bed_df, tol_type='intolerant')
	# Get features (including raw sequences) for most *tolerant* windows
	tol_merged_df, tol_filtered_onehot_seqs = compile_feature_table_per_tolerance_class(all_gwrvis_bed_df, tol_type='tolerant')
	

	print('===========================================')


	all_merged_df = pd.concat([intol_merged_df, tol_merged_df], axis=0)
	all_merged_df.reset_index(inplace=True, drop=True)
	all_merged_df.drop(['tolerance_rank'], axis=1, inplace=True)
	print(all_merged_df.shape)

	all_filtered_onehot_seqs = np.concatenate((intol_filtered_onehot_seqs, tol_filtered_onehot_seqs), axis=0)
	print(all_filtered_onehot_seqs.shape)

	split_data_into_train_val_test_sets(all_merged_df, all_filtered_onehot_seqs)
