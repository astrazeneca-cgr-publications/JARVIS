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
from scipy.stats import mannwhitneyu

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from custom_utils import create_out_dir, get_config_params, is_outlier




class GwrvisDistributionPerClass:

	def __init__(self, win_len=3000, class_type='all', input_classes=None):
		"""
			Genomic Classes initialisation
		"""
		# Read input classes to plot in current call

		
		if class_type == 'ccds':
			self.gwrvis_distr_dir = "/projects/cgr/users/kclc950/JARVIS/out/topmed-single-nt-gwrvis-winlen_" + str(win_len) + ".MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/gwrvis_distr_in_ccds_classes" 
		elif class_type == 'all':
			self.gwrvis_distr_dir = "/projects/cgr/users/kclc950/JARVIS/out/topmed-single-nt-gwrvis-winlen_" + str(win_len) + ".MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/gwrvis_distr_in_ALL_classes" 
		
		
		if input_classes is None:
			self.input_classes = [os.path.basename(x[0]) for x in os.walk(self.gwrvis_distr_dir) if 'gwrvis_distr_in' not in os.path.basename(x[0])]
			print(self.input_classes)
		else:
			self.input_classes = input_classes
	
		
		self.sorted_genomic_classes = sorted(self.input_classes)


		# Assign color to each genomic class, if not defined already
		col_palette = sns.color_palette("Set2", 4).as_hex()
		self.genomic_colors = {'intergenic': '#fdae6b', 'ccds': '#2ca25f', 'cds': '#2ca25f', 'vista': '#807dba', 'omim-HI': '#3182bd', 'intolerant': '#de2d26', 'tolerant': '#737373', 'ucne': '#cab2d6', 'enhancer': '#a6cee3', 'pathogenic': '#fb9a99', 'intron': '#756bb1', 'promoter': '#9a9a9a', 'clingen-hi': '#2171b5', 'clingen-nz': '#df65b0'}
		cnt = 0
		for name in self.sorted_genomic_classes:
			if name not in self.genomic_colors:
				self.genomic_colors[name] = col_palette[cnt]
				cnt += 1 
		


	def read_gwrvis_scores_per_class(self, discard_positive_gwrvis=False):

		self.gwrvis_lists = {}

		
		
		for cl in self.sorted_genomic_classes:
			cnt = 0
		
			print('>', cl)
			cur_class_gwrvis_file = self.gwrvis_distr_dir + '/' + cl + '/' + cl + '.gwrvis.bed'	

			cur_scores_list = []
			with open(cur_class_gwrvis_file) as fh:
				for line in fh:
					line = line.rstrip()
					_, _, _, score = line.split('\t')

					cur_scores_list.append(float(score))
					
					# DEBUG
					#cnt += 1
					#if cnt == 1000000:
					#	break

			cur_scores_list = np.array(cur_scores_list)
			cur_scores_list = cur_scores_list[ ~is_outlier(cur_scores_list) ]
			
			if discard_positive_gwrvis:
				cur_scores_list = cur_scores_list[ cur_scores_list <=0 ]
			
			# DEBUG - comment when debugging to reduce mem. requirements
			self.gwrvis_lists[cl] = cur_scores_list

			print('Class:', cl, ' - Median gwRVIS:', np.median(cur_scores_list))
			


			
	def compile_gwrvis_df(self):
		"""
			Compile gwRVIS data frame per genomic class
			Also return a dictionary with median gwRVIS per genomic class (for use with boxplot plotting)
		"""
		self.gwrvis_df = pd.DataFrame()
		self.gwrvis_medians_per_class = {}

		for name in self.sorted_genomic_classes[::-1]:
			print(name)
			self.gwrvis_medians_per_class[name] = np.median(self.gwrvis_lists[name])


		for name, median_gwrvis in sorted(self.gwrvis_medians_per_class.items(), key=lambda x: x[1]):
			tmp_df = pd.DataFrame({'gwrvis': self.gwrvis_lists[name], 'genomic_class': name})

			self.gwrvis_df = pd.concat([self.gwrvis_df, tmp_df], axis=0)
			print(self.gwrvis_df.shape)
			print('------\n')

		print(self.gwrvis_df.head())
		print(self.gwrvis_df.tail())
		print(self.gwrvis_df.shape)

			


		
	def calc_mannwhitney_u_between_classes(self):
	
		for i in range(len(self.sorted_genomic_classes)):
			for j in range(i, len(self.sorted_genomic_classes)):

				left_class = self.sorted_genomic_classes[i]
				right_class = self.sorted_genomic_classes[j]

				if left_class == right_class:
					continue


				_, pval = mannwhitneyu(self.gwrvis_lists[left_class], self.gwrvis_lists[right_class])

				print('Mann-Whitney U - ' + left_class + ' vs ' + right_class + ': ' + str(pval))


		


	def make_plots(self):


		print("PMFs (density plots)")

		plt.style.use('seaborn-whitegrid')
		dens_fig, ax1 = plt.subplots()
		for name in self.sorted_genomic_classes[::-1]:
			cur_gwrvis = self.gwrvis_lists[name]

			#cur_gwrvis = np.array(cur_gwrvis)
			#cur_gwrvis = cur_gwrvis[ ~is_outlier( cur_gwrvis ) ]
			print(cur_gwrvis)

			sns.kdeplot(cur_gwrvis, bw=0.5, ax=ax1, shade=False, color=self.genomic_colors[name])


		ax1.grid(linewidth=0.4)
		dens_fig.suptitle("", fontsize=12) 
		plt.xlabel('gwRVIS score', fontsize=10) 
		plt.ylabel('Density', fontsize=10)

		# legend
		patches = []
		for name in self.sorted_genomic_classes[::-1]:
			cur_patch = mpatches.Patch(color=self.genomic_colors[name], label=name + ' (' + str(len( self.genomic_colors[name] )) + ')')
			patches.append(cur_patch)
		plt.legend(handles = patches)


		dens_fig.savefig(self.gwrvis_distr_dir + '/' + 'gwRVIS_distr.density_plots.pdf', bbox_inches='tight')


	
		print("Plot boxplots in order of increasing median tolerance")

		ordered_genomic_classes = self.gwrvis_df.genomic_class.unique()
		print(ordered_genomic_classes)
		ordered_genomic_colors = [self.genomic_colors[cl] for cl in ordered_genomic_classes]

		boxplot_fig, ax = plt.subplots()
		sns.boxplot(data=self.gwrvis_df, 
				y='gwrvis', x='genomic_class',
				width=0.5,
				linewidth=0.4,	
				fliersize=0,
				palette=ordered_genomic_colors) #"colorblind")
		ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=6)
		ax.axhline(y=0, linestyle='-', linewidth=0.3, color='red')

		boxplot_fig.savefig(self.gwrvis_distr_dir + '/' + 'gwRVIS_distr.Boxplots.pdf', bbox_inches='tight')


				
	
	



if __name__ == '__main__':


	win_len = sys.argv[1]
	class_type = sys.argv[2] # ccds or all

	#input_classes = None
	input_classes = ['lincrna', 'intergenic']

	# =============================== MAIN ANALYSIS ===============================	
	obj = GwrvisDistributionPerClass(win_len=win_len, class_type=class_type, input_classes=input_classes)
	print('Sorted:', obj.sorted_genomic_classes)


	obj.read_gwrvis_scores_per_class()


	obj.compile_gwrvis_df()

	obj.calc_mannwhitney_u_between_classes()
	sys.exit()

	obj.make_plots()
