import pandas as pd
import sys, os

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from custom_utils import create_out_dir, get_config_params 



def update_clinvar_table(out_dir, patho_benign_sets):

	original_clinvar_table_file = out_dir + '/ml_data/clinvar_feature_tables/full_feature_table.' + patho_benign_sets + '.bed'

	single_nt_bed = original_clinvar_table_file + '.all_chr.single_nt_gwrvis.bed'



	clinvar_df = pd.read_csv(original_clinvar_table_file, sep='\t')
	clinvar_df['aux'] = clinvar_df['chr'].astype(str) + '_' + clinvar_df['start'].astype(str) + '_' + clinvar_df['end'].astype(str)
	print(clinvar_df.head())
	#print(clinvar_df.info())
	print(clinvar_df.shape)


	singlent_gwrvis_df = pd.read_csv(single_nt_bed, sep='\t', header=None)
	singlent_gwrvis_df.columns = ['chr', 'start', 'end', 'singlent_gwrvis']
	singlent_gwrvis_df['aux'] = singlent_gwrvis_df['chr'].astype(str) + '_' + singlent_gwrvis_df['start'].astype(str) + '_' + singlent_gwrvis_df['end'].astype(str)
	singlent_gwrvis_df.drop(['chr', 'start', 'end'], axis=1, inplace=True)
	#print(singlent_gwrvis_df.head())
	print(singlent_gwrvis_df.shape)



	full_df = clinvar_df.merge(singlent_gwrvis_df, left_on='aux', right_on='aux', how='left')
	full_df['gwrvis'] = full_df['singlent_gwrvis'].copy()
	full_df.drop(['aux', 'singlent_gwrvis'], axis=1, inplace=True)

	#print(full_df.info())
	print(full_df.head())
	print(full_df.shape)


	os.system("cp " + original_clinvar_table_file + " " + original_clinvar_table_file+".tiled-window-gwrvis.bed")

	full_df.to_csv(original_clinvar_table_file, sep='\t', index=False)




if __name__ == '__main__':

	config_file = sys.argv[1] 

	config_params = get_config_params(config_file)  

	# patho_benign set
	pathogenic_set = config_params['pathogenic_set']
	benign_set = config_params['benign_set']  
	patho_benign_sets = pathogenic_set + '_' + benign_set

	# output dir
	out_dir = create_out_dir(config_file)

	
	update_clinvar_table(out_dir, patho_benign_sets)
