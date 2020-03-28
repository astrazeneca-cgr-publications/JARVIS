import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'figure.max_open_warning': 0})
import sys
from subprocess import call
import subprocess
from io import StringIO
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import mode
from stats_util import is_outlier
import os
import seaborn as sns
import matplotlib.patches as mpatches
from copy import deepcopy
from pathlib import Path
import math

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from custom_utils import create_out_dir, get_config_params


def sigmoid(x):
	return 1 / (1 + math.exp(-0.5*x))


def parse_input_classes(input_classes_file):
	"""
	Read input classes to plot in current call
	"""
	input_classes = []
	with open(input_classes_file) as fh:
		for l in fh:
			if l.startswith('#'):
				continue
			else:
				input_classes.append( l.rstrip() )
	return input_classes


def get_sorted_genomic_classes(genomic_classes_log):

	cnt = 0
	genomic_classes = {}
	with open(genomic_classes_log) as fh:
		for l in fh:
			if cnt == 0: # skip header line
				cnt += 1
				continue

			l = l.rstrip()
			name, path, priority = l.split('\t')
			if name not in input_classes:
				continue
			genomic_classes[name] = [path, priority]

	# get sorted list of genomic classes based on their priority score
	sorted_genomic_classes = sorted(genomic_classes.items(), key=lambda x: int(x[1][1]), reverse=True)
	sorted_genomic_classes = [x[0] for x in sorted_genomic_classes]

	return sorted_genomic_classes



def map_genomic_classes_to_colors(sorted_genomic_classes):

	col_palette = sns.color_palette("Set2", 12).as_hex()
	genomic_colors = {'intergenic': '#fdae6b', 'ccds': '#2ca25f', 'cds': '#2ca25f', 'vista': '#807dba', 'omim-HI': '#3182bd', 'intolerant': '#de2d26', 'tolerant': '#737373', 'ucne': '#cab2d6', 'enhancer': '#a6cee3', 'pathogenic': '#fb9a99', 'intron': '#756bb1', 'promoter': '#9a9a9a', 'clingen-hi': '#2171b5', 'clingen-nz': '#df65b0'}

	# assign color to each genomic class
	cnt = 0
	for name in sorted_genomic_classes:
		if name not in genomic_colors:
			genomic_colors[name] = col_palette[cnt]
			cnt += 1

	return genomic_colors


def get_gwrvis_scores(sorted_genomic_classes):
	"""
	     Read gwRVIS scores per genomic class
	     Remove any genomic classes without hits and ammend 'sorted_genomic_classes'
	"""
	gwrvis_lists = {}

	class_hits_dict = {}
	for name in sorted_genomic_classes:
		class_hits_dict[name] = 0
		for chr in chroms:

			cur_class_rvis_file = gwrvis_distr_data_dir + '/chr' + chr + '.' + name + '.csv'	
			if not os.path.exists(cur_class_rvis_file):
				continue

			class_hits_dict[name] = 1
			cur_df = pd.read_csv(cur_class_rvis_file, header=0)
			cur_df_scores = list(cur_df.iloc[:, 1])
			gwrvis_lists[name] = gwrvis_lists.get(name, []) + cur_df_scores

	for name in sorted_genomic_classes:
		gwrvis_lists[name] = np.array(gwrvis_lists[name])

	no_hits_classes = [name for name,val in class_hits_dict.items() if val == 0]
	if len(no_hits_classes) > 0:
		print('[Warning] No hits classes:', no_hits_classes)

	sorted_genomic_classes = [name for name in class_hits_dict if name not in no_hits_classes]

	return gwrvis_lists, sorted_genomic_classes



def get_gwrvis_axis_limits(gwrvis_lists, sorted_genomic_classes):

	max_rvis = -sys.float_info.max
	min_rvis = sys.float_info.max

	x_axis_lim_left = -5
	x_axis_lim_right = 5

	if filter_plot_outliers:
		for name in sorted_genomic_classes:
			cur_rvis = gwrvis_lists[name]
			#print('Genomic class: ' + name)
			#is_outlier(cur_rvis, z_thres)
			
			cur_rvis = cur_rvis[ ~is_outlier(cur_rvis, z_thres) ]
			gwrvis_lists[name] = cur_rvis
		
			cur_max = max(cur_rvis)
			if cur_max > max_rvis:
				max_rvis = cur_max

			cur_min = min(cur_rvis)
			if cur_min < min_rvis:
				min_rvis = cur_min

		x_axis_lim_left = min_rvis
		x_axis_lim_right = max_rvis
	
	return x_axis_lim_left, x_axis_lim_right




def plot_density_plots(gwrvis_lists, sorted_genomic_classes):
	"""
	    Plot Density plots of gwRVIS per genomic class
	"""


	# ****** Plot PMFs (Density plots) ******
	length_per_class = {}
	#plt.style.use('seaborn-whitegrid')
	dens_fig, ax = plt.subplots()
	for name in sorted_genomic_classes[::-1]:
		cur_rvis = gwrvis_lists[name]
		length_per_class[name] = len(cur_rvis)
		sns.kdeplot(cur_rvis, bw=0.5, ax=ax, shade=False, color=genomic_colors[name])

	ax.set_xlim([x_axis_lim_left, x_axis_lim_right])
	#ax.set_ylim([0,0.8])
	dens_fig.suptitle('window: ' + str(win_len) + 'nt, ' + 'MAF: ' + str(MAF_thres), fontsize=12)
	plt.xlabel('rvis score', fontsize=10)
	plt.ylabel('density', fontsize=10)

	# legend
	patches = []
	for name in sorted_genomic_classes[::-1]:
		cur_patch = mpatches.Patch(color=genomic_colors[name], label=name + ' (' + str(length_per_class[name]) + ')')
		patches.append(cur_patch)

	plt.legend(handles = patches, fontsize=8)
	plt.show()
	
	return dens_fig



def plot_cdfs(gwrvis_lists, sorted_genomic_classes):
	"""
	    Plot CDF plots of gwRVIS per genomic class
	"""
	min_val = min([min(v) for k,v in gwrvis_lists.items() ])
	max_val = max([max(v) for k,v in gwrvis_lists.items() ])
	print("min, max:", min_val, max_val)

	binwidth = (max_val - min_val) / 1000
	bins = np.arange(min_val, max_val + binwidth, binwidth)

	cdf_fig, ax = plt.subplots()
	for name in sorted_genomic_classes[::-1]:
		cur_rvis = gwrvis_lists[name]

		#cur_rvis = [sigmoid(x) for x in cur_rvis]
		
		bins = np.arange(min(cur_rvis), max(cur_rvis) + binwidth, binwidth)
		plt.hist(cur_rvis, bins, cumulative=True, density=True, histtype='step', alpha=1.0, linewidth=0.8, label=name, color=genomic_colors[name])

	plt.legend(loc='upper left', frameon=True, facecolor='inherit', fontsize=8)
	ax.set_xlim([x_axis_lim_left, x_axis_lim_right])
	#ax.set_xlim([0, 1])
	plt.show()

	return cdf_fig




def save_full_gwrvis_per_class_df(gwrvis_lists, sorted_genomic_classes):
	"""
	    Save full gwRVIS data frame per genomic class into a csv file.
	    Also return a dictionary with median gwRVIS per genomic class (for use with boxplot plotting)
	"""
	gwrvis_df = pd.DataFrame()
	gwrvis_medians_per_class = {}

	for name in sorted_genomic_classes[::-1]:
		print(name)
		gwrvis_medians_per_class[name] = np.median(gwrvis_lists[name])


	for name, median_gwrvis in sorted(gwrvis_medians_per_class.items(), key=lambda x: x[1]):
		tmp_df = pd.DataFrame({'gwrvis': gwrvis_lists[name], 'genomic_class': name})

		gwrvis_df = pd.concat([gwrvis_df, tmp_df], axis=0)
		print(gwrvis_df.shape)
		print('------\n')

	print(gwrvis_df.head())
	gwrvis_df.to_csv(full_genome_out + '/Whole_genome_gwRVIS_per_class.csv', index=False)

	return gwrvis_df, gwrvis_medians_per_class



def make_boxplots(gwrvis_df, gwrvis_medians_per_class):
	"""
	    Plot Boxplots
	"""
	# Plot in order of increasing median tolerance

	ordered_genomic_classes = gwrvis_df.genomic_class.unique()
	print(ordered_genomic_classes)
	ordered_genomic_colors = [genomic_colors[cl] for cl in ordered_genomic_classes]

	boxplot_fig, ax = plt.subplots()
	sns.boxplot(data=gwrvis_df, 
		    y='gwrvis', x='genomic_class',
		    width=0.5,
		    linewidth=0.4,	
		    fliersize=0,
		    palette=ordered_genomic_colors) #"colorblind")
	ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=6)
	ax.axhline(y=0, linestyle='-', linewidth=0.3, color='red')

	pp = PdfPages(full_genome_out + "/Whole_genome_gwRVIS_Boxplots.pdf")
	pp.savefig(boxplot_fig)
	pp.close()





if __name__ == '__main__':

	config_file = sys.argv[1]	# 'config.log'
	input_classes_file = sys.argv[2]	# 'input_classes.txt'


	# ----- Init directories and global variables -----
	out_dir = create_out_dir(config_file)	# get output dir - already generated
	print('Output dir: ' + out_dir)


	gwrvis_distr_data_dir = out_dir + '/gwrvis_distribution/data'
	full_genome_out = out_dir + '/full_genome_out'
	ucne_dir = '../UCNEbase'
	
	chroms = [str(ch) for ch in list(range(1,23))]


	# ---- Read run parameters from config file and store into a dictionary -----
	run_params = get_config_params(config_file)

	genomic_classes_log = run_params['genomic_classes']
	filter_plot_outliers = run_params['filter_plot_outliers']
	z_thres = run_params['z_thres']
	win_len = run_params['win_len']
	MAF_thres = run_params['MAF_thres']
	filter_zero_variants_wins = run_params['filter_zero_variants_wins']


	# Read input genomic classes and assign colors to them
	input_classes = parse_input_classes(input_classes_file)
	sorted_genomic_classes = get_sorted_genomic_classes(genomic_classes_log)
	genomic_colors = map_genomic_classes_to_colors(sorted_genomic_classes)


	# Get gwRVIS scores
	gwrvis_lists, sorted_genomic_classes = get_gwrvis_scores(sorted_genomic_classes)
	print(sorted_genomic_classes)

	# Plot Density and CDF plots 
	x_axis_lim_left, x_axis_lim_right = get_gwrvis_axis_limits(gwrvis_lists, sorted_genomic_classes)

	dens_fig = plot_density_plots(gwrvis_lists, sorted_genomic_classes)
	cdf_fig = plot_cdfs(gwrvis_lists, sorted_genomic_classes)


	# Get full gwRVIS df per calss & median gwRVIS per class in a dict
	gwrvis_df, gwrvis_medians_per_class = save_full_gwrvis_per_class_df(gwrvis_lists, sorted_genomic_classes)

	# Make boxplots with gwRVIS per class (in descending order of intolerance)
	boxplot_fig = make_boxplots(gwrvis_df, gwrvis_medians_per_class)
	
	pp = PdfPages(full_genome_out + "/Whole_genome_gwRVIS_distribution.pdf")
	pp.savefig(dens_fig)
	pp.savefig(cdf_fig)
	pp.close()

	print('[INFO]: Whole genome RVIS run - complete.')
