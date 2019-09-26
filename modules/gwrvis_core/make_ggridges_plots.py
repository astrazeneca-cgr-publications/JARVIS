import matplotlib 
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys, os, getopt
from subprocess import call
import colorbrewer
from stats_util import is_outlier
from scipy.stats import mode
from pathlib import Path

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from custom_utils import create_out_dir


def main(argv):

	if len(argv) == 0:
		print("[Error] - Expected call format:\n    python post_process_results.py [ -o <outdir> | -c <config> ]")
		sys.exit(2)

	try:
		opts, args = getopt.getopt(argv,"ho:c:",["outdir=","config="])
	except getopt.GetoptError:
		print("[Error] - Expected call format:\n    python post_process_results.py [ -o <outdir> | -c <config> ]")
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print("python post_process_results.py [ -o <outdir> | -c <config> ]")
			sys.exit() 
		elif opt in ('-o', '--outdir'):
			out_dir = arg
		elif opt in ('-c', '--config'):
			config_file = arg
			out_dir = create_out_dir(config_file)
	return out_dir



def plot_barplot_with_genomic_elements(df):
	"""
	 Make BarPlot with number of elements in each genomic class
	"""

	x_pos = np.arange(df.shape[0])

	fig, ax = plt.subplots() 
	ax.set_title('Windows distribution across studied genomic classes')
	plt.bar(x_pos, df['counts'], color=colorbrewer.brewer['Set1']['9'] )

	for i, v in enumerate(df['counts']):
		ax.text(i - 0.5, v + 10000, str(v), color='black', size=8, family='serif')

	plt.xticks(x_pos, df['genom_class'], rotation='vertical')
	fig.tight_layout()

	plt.savefig(results_dir + '/res_barplot.pdf')



def get_density_plots(gen_df, annot='all'):

	print('Get density plots for: ' + annot)

	df = pd.DataFrame()
	for cl in gen_df['genom_class']:
		tmp_df = pd.read_csv(results_dir + '/full_genome.' + cl + '.csv', header=None, sep='\t', names=[cl])
		print(tmp_df.head())

		# filter  or retain outliers
		vec = tmp_df.ix[:, 0]

		if annot == 'no_outliers':
			vec = vec[ ~is_outlier(vec, 3.5)] 
		elif annot == 'only_outliers':
			vec = vec[ is_outlier(vec, 3.5)] 

		#vec = vec[ vec != intergenic_mode] # [BETA]: remove peak of windows with zero variants 

		tmp_df.ix[:, 0] = vec
		df = pd.concat([df, tmp_df], axis=1)

	print(df.head())
	out_file = 'whole_genome.rvis.data-frame.' + annot + '.csv'
	df.to_csv(results_dir + '/' + out_file, index=False)
	
	
	# Temp fix for native R installation on 'Ubuntu in Windows' that throws 'core dumped' exception
	# during patchwork library installation.
	rscript_path = 'Rscript'
	#rscript_w_patchwork_package = '/home/djifos/R_home/R-3.5.0/bin/Rscript'
	#rscript_w_patchwork_package_path = Path(rscript_w_patchwork_package)

	#if rscript_w_patchwork_package_path.exists():
	#	rscript_path = rscript_w_patchwork_package 

	print(rscript_path, 'gwrvis_core/make_beautiful_density_plots.R', results_dir, out_file, annot)
	call([rscript_path, 'gwrvis_core/make_beautiful_density_plots.R', results_dir, out_file, annot])



if __name__ == '__main__':

	out_dir = main(sys.argv[1:])
	print(out_dir)
	results_dir = out_dir + '/full_genome_out'

	gen_df = pd.read_csv(results_dir + '/genomic_class_elements.txt', header=None, sep='\t')
	gen_df.columns = ['genom_class', 'counts']
	gen_df.sort_values(by=['counts'], inplace=True, ascending=False)
	print(gen_df)

	# Make BarPlot with number of elements in each genomic class
	# DBG - temporarliy commented out:
	#plot_barplot_with_genomic_elements(gen_df)


	# [BETA]
	#cl = 'intergenic'
	#intergenic_df = pd.read_csv(results_dir + '/full_genome.' + cl + '.csv', header=None, sep='\t', names=[cl])
	#intergenic_mode = mode(intergenic_df[cl])[0][0]

	get_density_plots(gen_df, annot='all')
	get_density_plots(gen_df, annot='no_outliers')
	get_density_plots(gen_df, annot='only_outliers')

