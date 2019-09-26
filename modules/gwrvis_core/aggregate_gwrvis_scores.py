import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import pandas as pd
import numpy as np
import subprocess
import sys, os
import seaborn as sns;
from matplotlib.backends.backend_pdf import PdfPages
import operator

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from custom_utils import create_out_dir


def aggregate_scores_into_csv_per_class():
	
	print(">> Running aggregate_scores_into_csv_per_class()...")

	for cl in input_classes:
		print('>', cl)
		cur_list = []	

		tmp_out_file = full_genome_dir + '/full_genome.' + cl + '.csv' 
		for chr in chroms:
			#print(chr)	
			tmp_file = data_dir + '/chr' + str(chr) + '.' + cl + '.csv'
			if (not os.path.exists(tmp_file)) or (os.stat(tmp_file).st_size == 0):
				print('[!] No data for this class at chr:', chr)
				continue
			tmp_df = pd.read_csv(tmp_file, header=None, index_col=0)
		
			tmp_list = tmp_df[1].tolist()
			#print(len(tmp_list))	

			cur_list.extend(tmp_list)
			#print(len(cur_list))
		full_rvis_scores[cl] = cur_list

		with open(tmp_out_file, 'w') as fh:
			for l in cur_list:
				fh.write(str(l) + '\n')
		
	print(full_rvis_scores.keys())
	with open(full_genome_dir + '/genomic_class_elements.txt', 'w') as fh:
		for cl, vals in full_rvis_scores.items():
			fh.write(cl + '\t' + str(len(vals)) + '\r\n')

# convert P-value to Phred like score
def convert_pval_to_phred(pval):
	if pval == 0:
		pval = 10 ** (-320)

	phred = -10 * np.log10(float(pval))     
	return phred



# Mann-Whitney U Test
def run_mann_whitney_tests():
	
	annot_cols = [genomic_classes_dict[c] for c in input_classes]
	phreds_df = pd.DataFrame(columns = annot_cols)
	medians_dict = dict()

	print(">> In run_mann_whitney_tests()...")
	mann_whitn_u_test_str = ''
	fh = open(full_genome_dir + '/mann_whitney_u_test.txt', 'w')
	for i in range(0, len(input_classes)):
		cl1 = input_classes[i]
		new_row = []
		medians_new_row = []

		#for j in range(i+1, len(input_classes)):
		for j in range(0, len(input_classes)):
			cl2 = input_classes[j]

			pval = None
			phred = None
			if cl1 == cl2:
				pval = 1
				phred = convert_pval_to_phred(pval)
			else:	
				a = full_rvis_scores[cl1]
				b = full_rvis_scores[cl2]
				res = mannwhitneyu(a, b)
			
				if cl1 not in medians_dict:
					medians_dict[cl1] = np.median(a)

				print('\n' + cl1 + ' (median: ' + str(np.median(a)) + ') vs ' + cl2 + ' (median: ' + str(np.median(b)) + '):')
				print(res)

				pval = res.pvalue
				fh.write('> ' + cl1 + ' vs ' + cl2 + ':\n')
				fh.write('> ' + cl1 + ' (median: ' + str(np.median(a)) + ') vs ' + cl2 + ' (median: ' + str(np.median(b)) + '):\n')
				fh.write('[MannwhitneyuResult] statistic: ' + str(res.statistic) + ', P-value: ' + str(pval) + "\n\n")

				phred = convert_pval_to_phred(pval)
			
			new_row.append(phred)

		# append to overall phreds data frame
		phreds_df.loc[i] = new_row
		#print(phreds_df)

	phreds_df.index = annot_cols
	print(phreds_df)
	
	# sort classes by gwRVIS medians
	sorted_medians = sorted(medians_dict.items(), key=operator.itemgetter(1))
	medians = [t[1] for t in sorted_medians]
	classes = [t[0] for t in sorted_medians]
	classes = [genomic_classes_dict[c] for c in classes]
	print(classes)

	# sort indexes and columns in phreds data frame
	print(phreds_df.head())
	phreds_df = phreds_df.loc[:, classes ]
	phreds_df = phreds_df.reindex(index = classes)

	# plot heatmap with Phred like scores
	a4_dims = (14, 10)
	fig1, ax1 = plt.subplots(figsize=a4_dims)
	#ax.set_title("Phred scores of p-values from Mann Whittney U tests\nbetween all pairs of genomic classes")
	fig1.suptitle('Phred scores of p-values from Mann Whittney U tests\nbetween all pairs of genomic classes', fontsize=20)

	mask = np.zeros_like(phreds_df)
	mask[ np.triu_indices_from(mask) ] = True
	sns.heatmap(phreds_df, mask=mask, cmap="GnBu", square=True)

	fig2, ax2 = plt.subplots(figsize=a4_dims)

	y_pos = np.arange(len(medians))

	plt.bar(y_pos, medians)
	plt.xticks(y_pos, classes, rotation=90)
	plt.show()

	pp = PdfPages(full_genome_dir + '/Heatmap.phred_scores_from_mannwhitneyu.pdf') 
	pp.savefig(fig1, bbox_inches='tight')
	pp.savefig(fig2, bbox_inches='tight')
	pp.close()

	# close file with Mann-Whitney U test results
	fh.close()

		

# Concatenate values across all chromosomes for each of the functional classes into separate BED files. [DONE]
# ...to be used for comparison with Orion, CADD and CDTS.
def aggregate_scores_into_bed_files():

	print(">> Running aggregate_scores_into_bed_files()...")
	aggregate_bed_files_dir = full_genome_dir + '/BED'
	if (not os.path.exists(aggregate_bed_files_dir)):
		os.makedirs(aggregate_bed_files_dir)

	all_classes_df = pd.DataFrame()
	all_classes_midpoint_df = pd.DataFrame()

	for cl in input_classes:
		print('> ' + cl)
		cur_df = pd.DataFrame()	

		tmp_out_file = aggregate_bed_files_dir + '/full_genome.' + cl + '.bed' 
		for chr in chroms:
			tmp_file = bed_dir + '/gwrvis_scores_chr' + str(chr) + '.genomic_coords.' + cl + '.bed'
			if (not os.path.exists(tmp_file)) or (os.stat(tmp_file).st_size == 0):
				print('[!] No data for class "' + cl + '" at chr:', chr)
				continue
			tmp_df = pd.read_csv(tmp_file, header=None, sep='\t')

			cur_df = pd.concat([cur_df, tmp_df], axis=0)

		if cur_df.empty:                         
			continue
		
		cur_df.columns = ['chr', 'start', 'end', 'gwrvis']
		cur_df['genomic_class'] = cl
		cur_df.to_csv(tmp_out_file, sep='\t', header=False, index=False)
		all_classes_df = pd.concat([all_classes_df, cur_df], axis=0)
		print(all_classes_df.tail())
		print(all_classes_df.shape)


		middle_point_df = cur_df[ ['chr', 'start'] ].copy()
		middle_point_df.loc[:, 'start'] = (cur_df.loc[:, 'end'] + middle_point_df.loc[:, 'start'])/2
		middle_point_df['start'] = middle_point_df['start'].astype(int)
		middle_point_df['end'] = middle_point_df['start'].copy() + 1
		middle_point_df['gwrvis'] = cur_df['gwrvis'].copy()
		middle_point_df['genomic_class'] = cur_df['genomic_class'].copy()

		# store coordinates from middle points of windows assigned for each class into BED files
		middle_point_out_file = aggregate_bed_files_dir + '/mid_window_point.' + cl + '.bed'
		middle_point_df.to_csv(middle_point_out_file, sep='\t', header=False, index=False)
		all_classes_midpoint_df = pd.concat([all_classes_midpoint_df, middle_point_df], axis=0)

	# remove rows with NaN gwRVIS
	print('Shapes:', all_classes_df.shape, all_classes_midpoint_df.shape)
	all_classes_df.dropna(inplace=True)
	all_classes_midpoint_df.dropna(inplace=True)
	print('Non-NaN gwRVIS shapes:', all_classes_df.shape, all_classes_midpoint_df.shape)

	all_classes_df.to_csv(aggregate_bed_files_dir + '/full_genome.All_genomic_classes.bed', sep='\t', header=False, index=False)
	all_classes_midpoint_df.to_csv(aggregate_bed_files_dir + '/mid_window_point.All_genomic_classes.bed', sep='\t', header=False, index=False)
	



if __name__ == '__main__':
	
	config_file = sys.argv[1]
	input_classes_file = sys.argv[2]


	# ================== Initialisation and Dir definition/creation ==================
	# Read input classes to aggregate 
	input_classes = []
	with open(input_classes_file) as fh:
		for l in fh:
			if not l.startswith('#'):
				input_classes.append( l.rstrip() )
	print(input_classes)

	out_dir = create_out_dir(config_file)
	data_dir = out_dir + '/gwrvis_distribution/data'
	tmp_dir = out_dir + '/gwrvis_distribution/tmp'
	bed_dir = out_dir + '/gwrvis_distribution/BED'


	full_genome_dir = out_dir + '/full_genome_out'
	if not os.path.exists(full_genome_dir):
		os.makedirs(full_genome_dir)
	# ================================================================================



	# ----------- Global variables ------------
	chroms = list(range(1,23))

	genomic_classes_dict = {'ccds': 'rest of CCDS', 'tolerant': '25% most tolerant CCDS', 'intolerant': '25% most intolerant CCDS', 'utr': 'UTRs', 'intron': 'Introns', 'mature_mirna': 'miRNAs', 'ucne': 'UCNEs', 'vista': 'VISTA enhancers', 'omim-HI': 'OMIM-HI', 'pathogenic': 'ClinVar non-coding pathogenic', 'lincrna': 'lincRNAs', 'intergenic': 'Intergenic'}

	full_rvis_scores = {}
	# -----------------------------------------




	# ====================== MAIN ANALYSIS ======================
	sns.set()

	aggregate_scores_into_csv_per_class()

	run_mann_whitney_tests()

	aggregate_scores_into_bed_files()
	# ==========================================================-
