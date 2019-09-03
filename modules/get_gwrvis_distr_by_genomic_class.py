import matplotlib 
matplotlib.use('Agg') 
matplotlib.rcParams.update({'figure.max_open_warning': 0})
import sys, os
from subprocess import call
import subprocess
from io import StringIO
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import mode
import seaborn as sns
import matplotlib.patches as mpatches
from custom_utils import create_out_dir, get_config_params, is_outlier



def init_genomic_classes(input_classes_file, genomic_classes_log):
	"""
		Genomic Classes initialisation
	"""
	# Read input classes to plot in current call
	input_classes = []
	with open(input_classes_file) as fh:
		for l in fh:
			if l.startswith('#'):
				continue
			else:
				input_classes.append( l.rstrip() )
	print(input_classes)


	genomic_class_paths = {}
	classes_cnt = -1
	with open(genomic_classes_log) as fh:
		for l in fh:
			if classes_cnt == -1: # skip header line
				classes_cnt += 1
				continue 

			l = l.rstrip()	
			name, path, priority = l.split('\t')
			if name not in input_classes:
				continue
			genomic_class_paths[name] = [path, priority]

			classes_cnt += 1
	print('Num of classes:', classes_cnt)

	# Get sorted list of genomic classes based on their priority score
	sorted_genomic_classes = sorted(genomic_class_paths.items(), key=lambda x: int(x[1][1]), reverse=True) 
	sorted_genomic_classes = [x[0] for x in sorted_genomic_classes]


	# Assign color to each genomic class, if not defined already
	col_palette = sns.color_palette("Set2", 4).as_hex()
	genomic_colors = {'intergenic': '#fdae6b', 'ccds': '#2ca25f', 'cds': '#2ca25f', 'vista': '#807dba', 'omim-HI': '#3182bd', 'intolerant': '#de2d26', 'tolerant': '#737373', 'ucne': '#cab2d6', 'enhancer': '#a6cee3', 'pathogenic': '#fb9a99', 'intron': '#756bb1', 'promoter': '#9a9a9a', 'clingen-hi': '#2171b5', 'clingen-nz': '#df65b0'}
	cnt = 0
	for name in sorted_genomic_classes:
		if name not in genomic_colors:
			genomic_colors[name] = col_palette[cnt]
			cnt += 1 

	return sorted_genomic_classes, genomic_class_paths, genomic_colors





def get_gwrvis_for_subregion(ref_chr_gwrvis_bed, cur_genomic_class_bed, name, next_idx):

	# Perform intersectBed with or without '-wa' option based on input run parameters
	intersect_bed_cmd = './bedtools_wrapper.sh intersect'
	if fixed_win_len_for_all_genomic_elements:
		intersect_bed_cmd = './bedtools_wrapper.sh intersect_wa'

	local_min_overlap_ratio = min_overlap_ratio

	cmd = intersect_bed_cmd + ' ' + ref_chr_gwrvis_bed + ' ' + cur_genomic_class_bed + ' ' + local_min_overlap_ratio
	print(cmd)

	try:
		res = subprocess.check_output(cmd, shell=True)
		res = res.decode('utf-8')
		DATA = StringIO(res)

		# Store windows and respective gwRVIS scores for current genomic class and chromosome into a bed file
		gwrvis_for_cur_class = gwrvis_bed_dir + '/gwrvis_scores_chr' + chr + '.genomic_coords.' + name + '.bed'
		ff = open(gwrvis_for_cur_class, 'w')
		ff.write(DATA.getvalue())
		ff.close()
	except Exception as e:
		# set output string to empty
		DATA = StringIO()

	if DATA.getvalue() == '':
		return ref_chr_gwrvis_bed, -1


	# store all gwrvis values for current class in a list
	df = pd.read_csv(DATA, sep='\t', header=None)
	cur_gwrvis_vec = df[3]


	# Perform subtractBed with or without '-A' option based on input run parameters
	subtract_bed_cmd = './bedtools_wrapper.sh subtract'
	if fixed_win_len_for_all_genomic_elements: # and name == 'ccds':
		subtract_bed_cmd = './bedtools_wrapper.sh subtract_A'
		#local_min_overlap_ratio = '0.1'

	cmd = subtract_bed_cmd + ' ' + ref_chr_gwrvis_bed + ' ' + cur_genomic_class_bed + ' ' + local_min_overlap_ratio
	res = subprocess.check_output(cmd, shell=True)
	res = res.decode('utf-8')
	DATA = StringIO(res)

	new_ref_chr_gwrvis_bed = gwrvis_tmp_dir + '/gwrvis_scores_chr' + chr + '.genomic_coords.' + str(next_idx) + '.bed'
	ff = open(new_ref_chr_gwrvis_bed, 'w')
	ff.write(DATA.getvalue())
	ff.close()

	return new_ref_chr_gwrvis_bed, cur_gwrvis_vec





def extract_gwrvis_progressively_by_genomic_class(sorted_genomic_classes, genomic_class_paths):
	"""
		Annotation of genomic regions of certain classes occurs progressively,
		in a mutually exclusive / priority-based way (e.g. genic classes have higher priority).
		This means that parts of the windows from a chromosome get eliminated as they overlap
		with regions of genomic classes.
		
		However, there might be cases where a window is shared by two genomic classes:
		the one with the higher priority will get assigned the gwRVIS value of the window and
		eliminate the corresponding overlapping region. If the remaining part of the window still
		overlaps with another genomic class, this also gets assigned the gwRVIS value of that window too.
		
		So, calcualtions are mutually exclusive in terms of not counting overlapping regions from two different genomic classes twice,
		but allowing (non-overlapping) regions that share the same gwRVIS assigned to a window.
		As population size increases, window size can decrease more allow for more granularity and separation of such cases between
		different genomic classes.
	"""

	ref_chr_gwrvis_bed = gwrvis_dir + '/gwrvis.chr' + chr + '.genomic_coords.bed'
	print(ref_chr_gwrvis_bed)


	gwrvis_lists = {}
	names_to_del = []
	length_per_class = {}
	cnt = 1
	for name in sorted_genomic_classes:
		print('\n>> Processing class:', name)
		cur_genomic_class_bed = genomic_class_paths[name][0]
		print(cur_genomic_class_bed)

		ref_chr_gwrvis_bed, gwrvis_lists[name] = get_gwrvis_for_subregion(ref_chr_gwrvis_bed, cur_genomic_class_bed, name, cnt)

		if isinstance(gwrvis_lists[name], int) and gwrvis_lists[name]  == -1:
			print('[No elements found] - genomic class: ' + name)
			names_to_del.append(name)
		else:
			print('Length of current class (num. of windows):', len(gwrvis_lists[name]))
			length_per_class[name] = len(gwrvis_lists[name])

		cnt += 1
		print(gwrvis_lists.keys())

	# remove genomic classes with no data
	for name in names_to_del:
		del gwrvis_lists[name]
		sorted_genomic_classes.remove(name)

	return gwrvis_lists, length_per_class



def filter_out_gwrvis_nan(sorted_genomic_classes, gwrvis_lists):
	"""
		Filter out NaNs from RVIS score lists
	"""
	for name in sorted_genomic_classes:
		cur_gwrvis = gwrvis_lists[name]
		cur_gwrvis = cur_gwrvis[ ~np.isnan(cur_gwrvis) ]

		# BETA: Filter out gwRVIS scores from intergenic regions with no variation
		#if name == 'intergenic':
		#	intergenic_mode = mode(cur_gwrvis)[0][0]  # most frequent gwRVIS value in intergenic class
		#	cur_gwrvis = cur_gwrvis[ cur_gwrvis != intergenic_mode ]
		gwrvis_lists[name] = cur_gwrvis

		# BETA: Gives good results, but not scientifically plausible
		#cur_gwrvis = cur_gwrvis[ ~is_outlier(cur_gwrvis, z_thres) ]

		# Save gwRVIS values for current genomic class to a file
		cur_gwrvis.to_csv(gwrvis_data_dir + '/chr' + chr + '.' + name + '.csv', header=False)



def make_plots(sorted_genomic_classes, gwrvis_lists, length_per_class):

	#  Filter out outliers before plotting
	if filter_plot_outliers:
		for name in sorted_genomic_classes:
			cur_gwrvis = gwrvis_lists[name]
			cur_gwrvis = cur_gwrvis[ ~is_outlier(cur_gwrvis, z_thres) ]

			gwrvis_lists[name] = cur_gwrvis

	# PMFs (density plots)
	plt.style.use('seaborn-whitegrid')
	dens_fig, ax1 = plt.subplots()
	for name in sorted_genomic_classes[::-1]:
		cur_gwrvis = gwrvis_lists[name]
		sns.kdeplot(cur_gwrvis, bw=0.5, ax=ax1, shade=False, color=genomic_colors[name])

	#ax1.set_xlim([-2.5, 2.5])
	#ax1.set_ylim([0,0.8])
	ax1.grid(linewidth=0.4)
	dens_fig.suptitle('chr: ' + chr + ', window: ' + str(win_len) + 'nt, ' + 'MAF: ' + str(MAF_thres), fontsize=12) 
	plt.xlabel('RVIS score', fontsize=10) 
	plt.ylabel('Density', fontsize=10)

	# legend
	patches = []
	for name in sorted_genomic_classes[::-1]:
		cur_patch = mpatches.Patch(color=genomic_colors[name], label=name + ' (' + str(length_per_class[name]) + ')')
		patches.append(cur_patch)
	plt.legend(handles = patches)



	# CDFs
	#binwidth = 0.01
	min_val = min([min(v) for k,v in gwrvis_lists.items() ])
	max_val = max([max(v) for k,v in gwrvis_lists.items() ])
	abs_max = max(abs(min_val), abs(max_val))
	print(min_val, max_val)

	binwidth = (max_val - min_val) / 1000

	# min_val = min(min(cds_gwrvis), min(intergenic_gwrvis), min(ucne_gwrvis), min(intolerant_gwrvis), min(tolerant_gwrvis), min(omim_haploinsuf_gwrvis), min(enhancer_gwrvis))
	# max_val = max(max(cds_gwrvis), max(intergenic_gwrvis), max(ucne_gwrvis), max(intolerant_gwrvis), max(tolerant_gwrvis), max(omim_haploinsuf_gwrvis), max(enhancer_gwrvis))

	# DBG:
	#binwidth = max_val / min_val
	bins = np.arange(min_val, max_val + binwidth, binwidth)

	cdf_fig, ax2 = plt.subplots()
	for name in sorted_genomic_classes[::-1]:         
		cur_gwrvis = gwrvis_lists[name]         
		# == BETA == 
		#cur_gwrvis /= abs_max
		#cur_gwrvis.to_csv(gwrvis_data_dir + '/chr' + chr + '.' + name + '.csv')
		# ==========

		plt.hist(cur_gwrvis, bins, cumulative=True, density=True, histtype='step', alpha=1.0, linewidth=0.8, label=name, color=genomic_colors[name]) 

	plt.legend(loc='upper left', frameon=True, facecolor='inherit') 
	cdf_fig.suptitle('chr: ' + chr + ', window: ' + str(win_len) + 'nt, ' + 'MAF: ' + str(MAF_thres), fontsize=12) 
	plt.xlabel('RVIS score', fontsize=10) 
	plt.ylabel('Probability', fontsize=10)
	#ax2.set_xlim([-2.5, 2.5])
	ax2.grid(linewidth=0.4)
	plt.show()


	pp = PdfPages(gwrvis_distr_pdf_dir + "/gwRVIS_distribution.chr" + chr + ".pdf")
	pp.savefig(dens_fig)
	pp.savefig(cdf_fig)
	pp.close()




if __name__ == '__main__':

	args = sys.argv
	if len(args) != 4:
		print("Error: insufficient input arguments. Expected command call:\n> python get_gwrvis_distr_by_genomic_class.py [chr] [config_file] [input_classes]")
		sys.exit()
	 
	chr = args[1]
	config_file = args[2] #'config.log'
	input_classes_file = args[3]


	# ========================= Parameter Initialisation =========================
	# read run parameters from config file and store into a dictionary
	run_params = get_config_params(config_file)
	print(run_params)

	genomic_classes_log = run_params['genomic_classes']
	print(genomic_classes_log)


	filter_plot_outliers = run_params['filter_plot_outliers']
	hg_dir = run_params['hg_dir']
	z_thres = run_params['z_thres']
	win_len = run_params['win_len']
	MAF_thres = run_params['MAF_thres']
	fixed_win_len_for_all_genomic_elements = run_params['fixed_win_len_for_all_genomic_elements']
	min_overlap_ratio = str(run_params['min_overlap_ratio'])


	ucne_dir = '../other_datasets/UCNE_base'

	print('filter_plot_outliers:', filter_plot_outliers)
	print(hg_dir)

	out_dir = create_out_dir(config_file)
	print(out_dir)
	# ===========================================================================


	# ---------------------------- Dir initialisation ---------------------------
	gwrvis_dir = out_dir + '/gwrvis_scores'

	gwrvis_distr_dir = out_dir + '/gwrvis_distribution'
	if not os.path.exists(gwrvis_distr_dir):
		os.makedirs(gwrvis_distr_dir, exist_ok=True)
	gwrvis_tmp_dir = gwrvis_distr_dir + '/tmp'
	if not os.path.exists(gwrvis_tmp_dir):
		os.makedirs(gwrvis_tmp_dir, exist_ok=True)

	# Store BED files with gwRVIS for each genomic class
	gwrvis_bed_dir = gwrvis_distr_dir + '/BED'
	if not os.path.exists(gwrvis_bed_dir):
		os.makedirs(gwrvis_bed_dir, exist_ok=True)

	# save gwRVIS arrays to merge for all chromosomes later on
	gwrvis_data_dir = gwrvis_distr_dir + '/data'
	if not os.path.exists(gwrvis_data_dir):
		os.makedirs(gwrvis_data_dir, exist_ok=True)

	gwrvis_distr_pdf_dir = gwrvis_distr_dir + '/pdf'
	if not os.path.exists(gwrvis_distr_pdf_dir):
		os.makedirs(gwrvis_distr_pdf_dir, exist_ok=True)
	# ---------------------------------------------------------------------------

	


	# =============================== MAIN ANALYSIS ===============================	
	sorted_genomic_classes, genomic_class_paths, genomic_colors = init_genomic_classes(input_classes_file, genomic_classes_log)
	print(sorted_genomic_classes)


	gwrvis_lists, length_per_class = extract_gwrvis_progressively_by_genomic_class(sorted_genomic_classes, genomic_class_paths)
	# Make sure if I keep all windows or exclude those (and their indexes) that have NaN
	

	filter_out_gwrvis_nan(sorted_genomic_classes, gwrvis_lists)

	# TEMP - DEBUG: temporarily commented to speed up total run time of wgs.sh
	#make_plots(sorted_genomic_classes, gwrvis_lists, length_per_class)

	print('[INFO]: Run for chr ' + chr + ' - complete.')
