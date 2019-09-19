from custom_utils import create_out_dir
import pandas as pd
import numpy as np
import subprocess
import sys, os



def merge_gwrvis_and_feature_tables():
	total_features_df = pd.read_csv(out_dir + '/tmp/total_df.Xy.tsv', sep=' ')
	gwrvis_bed_df = pd.read_csv(out_dir + '/gwrvis_scores/full_genome.all_gwrvis.bed', sep='\t')
	

	total_features_df.rename(columns={'y': 'common_variants'}, inplace=True)
	total_features_df.index.name = 'chr_win'
	total_features_df.reset_index(inplace=True)

	gwrvis_bed_df['chr_win'] = gwrvis_bed_df['chr'].astype(str) + '_' + gwrvis_bed_df['win_index'].astype(str)


	gwrvis_features_df = gwrvis_bed_df.merge(total_features_df, how='inner', left_on=['chr_win'], right_on=['chr_win'])
	
	tmp_win_index = gwrvis_features_df['win_index'].copy()
	gwrvis_features_df.drop(['win_index', 'chr_win'], axis=1, inplace=True)
	print(gwrvis_features_df.head())
	print(gwrvis_features_df.shape)
	
	WRITTEN_DF_TO_FILE = False
	while not WRITTEN_DF_TO_FILE:
		try:
			gwrvis_features_bed = out_dir + '/full_genome_out/gwrvis_features_table.bed'
			gwrvis_features_df.to_csv(gwrvis_features_bed, sep='\t', index=None, header=None)
		
			with open(gwrvis_features_bed + '.header', 'w') as fh:
				fh.write('\t'.join(gwrvis_features_df.columns.values))
			WRITTEN_DF_TO_FILE = True
		except Exception as e:
			print(e)
			pass		

	return gwrvis_features_df, gwrvis_features_bed


	
	
def intersect_each_regulatory_feature(gwrvis_features_bed):

	# intersect each feature with all gwRVIS windows
	for regul_elem in regulatory_elements:
		
		print('\n> Extracting features for', regul_elem, '...')
		ensembl_regul_file = '../ensembl/GRCh37-Regulatory_Features/' + cell_line + '.' + regul_elem + '.sorted.bed'

	
		tmp_feature_gwrvis_intersection = scratch_dir + '/gwrvis_' + regul_elem + '.feature_table.bed'
		cmd = 'intersectBed -wao -a ' + gwrvis_features_bed + ' -b ' + ensembl_regul_file + """ | awk '!seen[$1"_"$2] {print} {++seen[$1"_"$2]}' | cut -f1,2,3,23 > """ + tmp_feature_gwrvis_intersection
		#print(cmd)
		
		p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout, stderr = p.communicate()
		res = str(stdout, "utf-8").rstrip()

		
		
def aggreate_all_regulatory_features():

	all_regul_features_df = pd.DataFrame()

	# aggregate all regulatory features across all windows
	print('\n\n Aggregating all regulatory features ...')
	
	cnt = 1
	for regul_elem in regulatory_elements:
		print(regul_elem)
		
		gwrvis_elem_intersection_df = pd.read_csv(scratch_dir + '/gwrvis_' + regul_elem + '.feature_table.bed', sep='\t', header=None)
		gwrvis_elem_intersection_df.fillna(0, inplace=True)
		gwrvis_elem_intersection_df.columns = ['chr', 'start', 'end', 'regul_elem']
		
		inferred_regul_elem = [x for x in gwrvis_elem_intersection_df['regul_elem'].unique() if x != 0][0]
		
		print(inferred_regul_elem)
		gwrvis_elem_intersection_df.replace({inferred_regul_elem: 1}, inplace=True)
		gwrvis_elem_intersection_df.rename(columns={'regul_elem': inferred_regul_elem}, inplace=True)
		
		#gwrvis_elem_intersection_df.index = gwrvis_elem_intersection_df['chr'].astype(str) + '_' + \
		#									gwrvis_elem_intersection_df['start'].astype(str) + '_' + \
		#									gwrvis_elem_intersection_df['end'].astype(str)
		print(gwrvis_elem_intersection_df.head(10))
		print('Non zero elements:', sum(gwrvis_elem_intersection_df[inferred_regul_elem]))
		
		#all_regul_features_df = pd.concat([all_regul_features_df, gwrvis_elem_intersection_df], axis=1)
		
		if all_regul_features_df.shape[0] > 0:
			all_regul_features_df = pd.merge(all_regul_features_df, gwrvis_elem_intersection_df, how='inner', left_on=['chr', 'start', 'end'], right_on=['chr', 'start', 'end'])
		else:
			all_regul_features_df = gwrvis_elem_intersection_df 
		
		print('All regul df:', all_regul_features_df.shape, '\n')
		
		
	print(all_regul_features_df.head())
	print(all_regul_features_df.shape)
		
	return all_regul_features_df
		

		
def merge_all_feature_tables(gwrvis_features_df, all_regul_features_df):

	print(gwrvis_features_df.head())
	print(all_regul_features_df.head())
	
	print('------------\n')
	full_df  = pd.merge(gwrvis_features_df, all_regul_features_df, how='inner', left_on=['chr', 'start', 'end'], right_on=['chr', 'start', 'end'])
	
	full_df.to_csv(out_dir + '/full_genome_out/full_gwrvis_and_regulatory_features.bed', sep='\t', index=None)
	print(full_df.head())
		

if __name__ == '__main__':
	
	config_file = sys.argv[1]
	
	# ------------ Initialisation ------------
	out_dir = create_out_dir(config_file)
	
	scratch_dir = '../scratch'
	if not os.path.exists(scratch_dir):
		os.makedirs(scratch_dir)

	full_genome_dir = out_dir + '/full_genome_out'         
	if not os.path.exists(full_genome_dir):
		os.makedirs(full_genome_dir)

		
	regulatory_elements = ['CTCF_binding_sites', 'Enhancers',
				'Open_chromatin', 'TF_binding_sites', 'H3K27ac',
				'H3K27me3', 'H4K20me1', 'H3K9ac', 'H3K4me1',	
				'H3K4me2', 'H3K4me3', 'H3K36me3']
	cell_line = 'Monocytes_CD14plus'
	# ---------------------------------------------
	
	
	
	gwrvis_features_df, gwrvis_features_bed = merge_gwrvis_and_feature_tables()
	#gwrvis_features_bed = '../out/gnomad-regression_beta-winlen_3000.MAF_0.001.varType_snv.Pop_all-PASS_ONLY-NO_SEGDUP-NO_LCR-high_conf_regions/full_genome_out/gwrvis_features_table.bed'
	
	intersect_each_regulatory_feature(gwrvis_features_bed)
	
	all_regul_features_df = aggreate_all_regulatory_features()

	merge_all_feature_tables(gwrvis_features_df, all_regul_features_df)
