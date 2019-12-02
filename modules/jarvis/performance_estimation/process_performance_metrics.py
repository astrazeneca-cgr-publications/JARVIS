import matplotlib 
matplotlib.use('agg') 
import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np 
import pandas as pd
import pickle
from glob import glob 

import sys,os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))
import custom_utils




class MetricsBenchmark:

	def __init__(self, config_file, genomic_class):
		config_params = custom_utils.get_config_params(config_file)
		self.genomic_class = genomic_class
		
		self.tables_per_metric = {}

		self.current_palette = sns.color_palette() + sns.color_palette("Paired")
		self.hex_colors = [matplotlib.colors.to_hex(x) for x in self.current_palette]



	def read_metrics_from_saved_files(self):

		# > All but DNN metrics (also includes JARVIS trained with RF; Random Forests)
		with open(clinvar_feature_table_dir + '/all_but_dnn_performance_metrics.' + self.genomic_class + '.pkl', 'rb') as handle:
			self.all_but_dnn_metrics = pickle.load(handle)
		# replace 'jarvis' name with 'JARVIS-rf' (rf: random forest)
		self.all_but_dnn_metrics['JARVIS-rf'] = self.all_but_dnn_metrics.pop('jarvis')


		# > DNN metrics (JARVIS: structured, sequences & both)
		try:
			self.jarvis_metrics = {}
			for input_features in ['structured', 'sequences', 'both']:
				with open(clinvar_feature_table_dir + '/jarvis_performance_metrics.' + input_features + '.' + self.genomic_class + '.pkl', 'rb') as handle:
					self.jarvis_metrics['JARVIS-' + input_features] = pickle.load(handle)

			print(self.all_but_dnn_metrics.keys())
			print(self.jarvis_metrics.keys())
		except:
			print("No DNN-based metrics found.")



	def get_melt_df_from_metric_dict(self, metrics_dict):

		cur_df = pd.DataFrame()

		for score in metrics_dict.keys():
			metric_values = [d[metric] for d in metrics_dict[score]]

			#print('Score:', score)
			#print(metrics_dict[score])
			#print(metric_values)

			tmp_df = pd.DataFrame({'score': score, 'values': metric_values})


			if cur_df.shape[0] > 0:
				cur_df = pd.concat([cur_df, tmp_df])
			else:
				cur_df = tmp_df

		cur_df['genomic_class'] = self.genomic_class
		cur_df.reset_index(inplace=True, drop=True)
		cur_df.fillna(0.5, inplace=True)
		cur_df.replace('NA', 0.5, inplace=True)

		#cur_df['values'] = pd.to_numeric(cur_df['values'])
		
		return cur_df


	def get_order_of_genomic_classes(self, df, order='mean'):

		#print(df)
		if order == 'mean':
			df = df.groupby(['score']).mean().reset_index().sort_values(by='values', ascending=False)
		elif order == 'median':
			df = df.groupby(['score']).median().reset_index().sort_values(by='values', ascending=False)

		return df


	def plot_metrics_per_score(self, metric='auc'):
		"""
		    > Grouped boxplots:
		    https://stackoverflow.com/questions/52392438/plotting-box-plots-of-two-columns-side-by-side-in-seaborn
		"""
		
		all_but_dnn_df = self.get_melt_df_from_metric_dict(self.all_but_dnn_metrics)
		jarvis_df = self.get_melt_df_from_metric_dict(self.jarvis_metrics)


		#print(list(all_but_dnn_df.keys()))
		#print(list(jarvis_df.keys()))

		self.all_scores = sorted(all_but_dnn_df.score.unique().tolist() + jarvis_df.score.unique().tolist())
		self.cur_palette = dict(zip(self.all_scores, self.hex_colors[:len(self.all_scores)]))


		cur_df = pd.concat([all_but_dnn_df, jarvis_df])
		cur_df.reset_index(inplace=True, drop=True)
		#print(cur_df)
		#print(cur_df.shape)

		ordered_classes_df = self.get_order_of_genomic_classes(cur_df, order='mean')
		#print(ordered_classes_df)

		sorted_df = pd.merge(ordered_classes_df.score, cur_df, left_on='score', right_on='score', how='outer')
		#print(sorted_df)


		# Get string for avg. metric values per score
		score_names = ordered_classes_df.score.astype(str).values.tolist()
		avg_values = ordered_classes_df['values'].values.tolist()
		avg_metric_per_score = [score_names[i]+': '+str(avg_values[i]) for i in range(len(score_names))]
		avg_metric_per_score = '> Avg. '+metric+'\n' + '\n'.join(avg_metric_per_score)


		plt.figure() #figsize=(20, 20))
		ax = sns.boxplot(data=sorted_df, x='genomic_class', y='values', hue='score', palette=self.cur_palette, linewidth=0.5)
		plt.legend(title='score', loc='upper left', bbox_to_anchor=(1, 1))
		plt.ylabel('Avg. ' + metric)
		ax.text(0.9, 0, avg_metric_per_score, fontsize=9)

		plot_filepath = clinvar_out_dir + '/Benchmarking-' + metric + '.' + genomic_class + '.pdf'
		ax.get_figure().savefig(plot_filepath, format='pdf', bbox_inches='tight')
		plt.close()


		



def infer_avail_genomic_classes():

	# Infer all run genomic_classes from "jarvis_performance_metrics.[*].[genomic_classes].pkl files
	jarvis_metric_files = glob(clinvar_feature_table_dir + '/jarvis_performance_metrics.both.*.pkl')
	print('All metrics files (with both features run complete):', jarvis_metric_files)

	genomic_classes = []
	for f in jarvis_metric_files:
		vals = f.split('.')
		genomic_classes.append(vals[-2])

	return genomic_classes 



if __name__ == '__main__':

	config_file = sys.argv[1]

	# Read dir structure
	out_dir = custom_utils.create_out_dir(config_file)
	ml_data_dir = out_dir + '/ml_data' 		
	clinvar_feature_table_dir = ml_data_dir + '/clinvar_feature_tables'
	clinvar_out_dir = ml_data_dir + '/clinvar-out'

	# Infer available genomic classes
	genomic_classes = infer_avail_genomic_classes()
		

	for genomic_class in genomic_classes:
		print('\n> Genomic class:', genomic_class)

		bench = MetricsBenchmark(config_file, genomic_class)
		bench.read_metrics_from_saved_files()

		metrics = bench.all_but_dnn_metrics['gwrvis'][0].keys()

		for metric in metrics:
			print('Metric:', metric)
			bench.plot_metrics_per_score(metric)
	
