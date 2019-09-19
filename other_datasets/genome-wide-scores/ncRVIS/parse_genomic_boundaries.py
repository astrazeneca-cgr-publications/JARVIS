import sys, os
import pandas as pd
import numpy as np
import re


def read_ncrvis_scores_per_gene(scores_file):

	ncrvis_per_gene = {}

	cnt = 0
	with open(scores_file) as fh:
		for line in fh:
			if cnt == 0:
				cnt += 1
				continue 

			vals = line.split(',')
	
			# Keeping the CCDS 15 release for gene nomenclature
			gene = vals[1]
			ncrvis = vals[4]
			if ncrvis == '-':
				continue
				
			ncrvis_per_gene[gene] = ncrvis

	return ncrvis_per_gene




def convert_boundaries_to_bed(boundaries_file, ncrvis_per_gene):

	intervals_per_gene = {}

	out_file = 'All_chromosomes.ncRVIS.tmp'
	out_fh = open(out_file, 'w')

	with open(boundaries_file) as fh:
		for line in fh:
			# C9orf163_5prime 9 (139378659..139378659,139378663..139378728,139378730..139378732,139378734..139378735,139378737..139378745,139378790..139378898) 190

			line = line.rstrip()

			gene, chrom, intervals, _ = line.split(' ')
			
			gene = re.sub('_.*', '', gene)
			if gene not in ncrvis_per_gene:
				continue
			else:
				cur_ncrvis = ncrvis_per_gene[gene]

			chrom = 'chr' + chrom


			intervals = re.sub('\(|\)', '', intervals)

			cur_interval_set = intervals.split(',')		
			for interv in cur_interval_set:
				interv_start, interv_end = interv.split('..')

				# convert to 0-based for BED file
				interv_start = str(int(interv_start) - 1)
				new_interv = interv_start + '\t' + interv_end
			
				out_line = chrom + '\t' + new_interv + '\t' + cur_ncrvis
				out_fh.write(out_line + '\n')

	out_fh.close()			
	return out_file


		
def post_process_bed_file(tmp_bed_file):

	merged_tmp_bed_file = tmp_bed_file + '.merged'
	os.system("cat " + tmp_bed_file + " | sortBed | mergeBed -c 4 -o 'collapse' > " + merged_tmp_bed_file)

	df = pd.read_csv(merged_tmp_bed_file, sep='\t', header=None, index_col=None)
	df.columns = ['chr', 'start', 'end', 'ncrvis']
	print(df.head())
	print(df.shape)

	df['ncrvis'] = df['ncrvis'].apply(lambda x: np.mean([float(k) for k in str(x).split(',')]) )
	
	df.to_csv('All_chromosomes.ncRVIS.bed', sep='\t', header=False, index=None)



if __name__ == '__main__':

	scores_file = 'original_ncRVIS_scores_table.csv'
	ncrvis_per_gene = read_ncrvis_scores_per_gene(scores_file)
	

	boundaries_file = 'genomic_boundaries.txt'
	tmp_bed_file = convert_boundaries_to_bed(boundaries_file, ncrvis_per_gene)


	post_process_bed_file(tmp_bed_file)
