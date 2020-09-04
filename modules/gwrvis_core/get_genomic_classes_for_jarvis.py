import pandas as pd
import sys, os

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from custom_utils import create_out_dir, get_config_params


def read_input_classes(input_classes_file):
	
	input_classes = []
	with open(input_classes_file) as fh:
		for l in fh:
			if not l.startswith('#'):
				input_classes.append( l.rstrip() )

	return input_classes



if __name__ == '__main__':

	config_file = sys.argv[1]
	input_classes_file = sys.argv[2]



	out_dir = create_out_dir(config_file)
	run_params = get_config_params(config_file)


	# read input classes
	input_classes = read_input_classes(input_classes_file)


	pathogenic_set = run_params['pathogenic_set']
	benign_set = run_params['benign_set']
	patho_benign_sets = pathogenic_set + '_' + benign_set


	clinvar_tables_dir = out_dir + '/ml_data/clinvar_feature_tables/'


	full_clinvar_feature_table = clinvar_tables_dir + './full_feature_table.' + patho_benign_sets + '.bed'
	full_df = pd.read_csv(full_clinvar_feature_table, sep='\t')

	full_df['aux'] = full_df['chr'] + '-' + full_df['start'].astype(str) + '-' + full_df['end'].astype(str)
	print(full_df.head())
	print(full_df.columns)




	for genomic_class in reversed(input_classes):

		print('>> ', genomic_class)

		cur_class_bed = clinvar_tables_dir + './full_feature_table.' + patho_benign_sets + '.bed.tmp.' + genomic_class + '_tmp'

		#print(cur_class_bed)

		cur_df = pd.read_csv(cur_class_bed, sep='\t', header=None)
		
		#cur_df = cur_df.iloc[ :, [0,1,2]].copy()
		#cur_df.columns = ['chr', 'start', 'end']
		#cur_df['genomic_class'] = genomic_class
		
		cur_df['aux'] = cur_df.iloc[:, 0] + '-' + cur_df.iloc[:, 1].astype(str) + '-' + cur_df.iloc[:, 1].astype(str)
		print(cur_df.head())
		print(cur_df.shape)


		cur_coords = cur_df['aux'].values
		full_df.loc[ full_df['aux'].isin(cur_coords), 'genomic_class'] = genomic_class
		print(full_df.shape)



	del full_df['aux']
	print('Final full_df:', full_df.shape)	

	full_df.to_csv(full_clinvar_feature_table, sep='\t', index=False)
