import pandas as pd
import os, sys
pd.set_option('display.max_columns', 40)
pd.set_option('display.width', 1000)



def read_denovodb_dfs(denovo_db_type):

	denovo_columns = ['SampleID', 'StudyName', 'PubmedID', 'NumProbands', 'NumControls', 'SequenceType', 'PrimaryPhenotype', 'Validation', 'Chr', 'Position', 'Variant', 'rsID', 'DbsnpBuild', 'AncestralAllele', '1000GenomeCount', 'ExacFreq', 'EspAaFreq', 'EspEaFreq', 'Transcript', 'codingDnaSize', 'Gene', 'FunctionClass', 'cDnaVariant', 'ProteinVariant', 'Exon/Intron', 'PolyPhen(HDiv)', 'PolyPhen(HVar)', 'SiftScore', 'CaddScore', 'LofScore', 'LrtScore']

	print("\n> Reading original " + denovo_db_type + ' file...')
	df = pd.read_csv('denovo-db.' + denovo_db_type + '-samples.variants.tsv.gz', compression='gzip', sep='\t', comment='#', header=None, low_memory=False)
	df.columns = denovo_columns
	print(df.shape)

	return df



def select_cols_and_convert_to_bed(df, denovo_db_type):

	print("\n> Processing " + denovo_db_type + ' file...')

	df['Start'] = df['Position'] - 1
	df['End'] = df['Position']
	df['Chr'] = 'chr' + df['Chr'].astype(str)

	selected_columns = ['Chr', 'Start', 'End', 'PrimaryPhenotype', 'FunctionClass', 'PubmedID', 'NumProbands', 'NumControls', 'SequenceType', 'Transcript', 'codingDnaSize', 'Gene', 'Variant', 'PolyPhen(HDiv)', 'PolyPhen(HVar)', 'SiftScore', 'CaddScore', 'LofScore', 'LrtScore']

	df = df[selected_columns]

	print(df.shape)

	# Filter out duplicate rows (e.g. with support from multiple samples/studies)
	groupby_cols = ['Chr', 'Start', 'End', 'PubmedID', 'PrimaryPhenotype', 'Gene', 'Variant', 'FunctionClass']
	df = df.groupby(groupby_cols).first().reset_index() 
	print(df.head())

	# Save to BED file
	if not os.path.exists('bed/'):
		os.makedirs('bed/')

	out_file = 'bed/denovo-db.' + denovo_db_type + '-samples.variants.bed'
	print('Saving denovo-db ' + denovo_db_type + ' file to ' + out_file)
	print(df.shape)

	df.to_csv(out_file, sep='\t', header=False, index=False)



if __name__ == '__main__':

	# read original denovo-db tables
	ssc_df = read_denovodb_dfs('ssc')
	non_ssc_df = read_denovodb_dfs('non-ssc')

	# susbset denovo-db files (SSC or non-SSC) and save into BED files
	select_cols_and_convert_to_bed(ssc_df, 'ssc')
	select_cols_and_convert_to_bed(non_ssc_df, 'non-ssc')
