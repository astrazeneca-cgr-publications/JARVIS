import pandas as pd
import sys, os
import subprocess
sys.path.insert(1, os.path.join(sys.path[0], '..'))
import custom_utils


def read_all_gwrvis_scores(all_gwrvis_bed_file):

	all_gwrvis_bed_df = pd.read_csv(all_gwrvis_bed_file, header=0, sep='\t')
	print(all_gwrvis_bed_df.head())
	print(all_gwrvis_bed_df.shape)

	# remove NAs
	all_gwrvis_bed_df.dropna(axis=0, inplace=True)
	print(all_gwrvis_bed_df.shape)

	return all_gwrvis_bed_df


def get_most_intol_or_toler_sequences(all_gwrvis_bed_df, tol_type='intolerant', top_ratio=0.1):
	
	
	# sort by gwRVIS value in ascending (for intolerant) or descending (for tolerant) order
	ascending = True
	if tol_type == 'tolerant':
		ascending = False
	all_gwrvis_bed_df.sort_values(by='gwrvis', inplace=True, ascending=ascending)
	all_gwrvis_bed_df.reset_index(drop=True, inplace=True)
	print(all_gwrvis_bed_df.head())
	print(all_gwrvis_bed_df.tail())


	top_seqs_set_size = int(all_gwrvis_bed_df.shape[0] * top_ratio)
	print(top_seqs_set_size)


	windows_file = seq_out_dir + '/most_' + tol_type + '.windows.txt'
	out_fh = open(windows_file, 'w')
	raw_seq_out_file = seq_out_dir + '/most_' + tol_type + '.raw_seqs.out' 
	#out_fh.write('chr\twin_index\tgwrvis\tseq\n'.encode())
	

	seq_list_df = all_gwrvis_bed_df[:top_seqs_set_size]
	seq_list_df = seq_list_df[['chr', 'start', 'end']]

	seq_list_df['chr'] = seq_list_df['chr'].str.replace('chr', '')
	seq_list = seq_list_df['chr'].map(str) + ':' + seq_list_df['start'].astype(str) + '-' + seq_list_df['end'].astype(str)
	out_fh.write("\n".join(seq_list))
	out_fh.close()	

	
	awk_fasta_collapser = "awk '/^>/ {printf(\"\\n%s\\n\",$0);next; } { printf(\"%s\",$0);}  END {printf(\"\\n\");}' | tail -n +2"
	cmd = '../util/twoBitToFa ' + human_ref_genome_2bit + ' -seqList=' + windows_file + " /dev/stdout | " + awk_fasta_collapser + " > " + raw_seq_out_file # | grep -v '>' | tr '\n' ' ' | sed 's/ //g'"
	print('\n' + cmd)

	p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = p.communicate()
	print(str(stderr, "utf-8").rstrip())
	
	



if __name__ == '__main__':

	config_file = sys.argv[1]

	human_ref_genome_2bit = '../../hg19/homo_sapiens_GRCh37_FASTA/hsa37.2bit'

	out_dir = custom_utils.create_out_dir(config_file)
	out_dir = '../' + out_dir

	gwrvis_scores_dir = out_dir + '/gwrvis_scores'
	seq_out_dir = gwrvis_scores_dir + '/raw_seq'
	if not os.path.exists(seq_out_dir):
		os.makedirs(seq_out_dir)

	all_gwrvis_bed_file = gwrvis_scores_dir + '/full_genome.all_gwrvis.bed'

	# Read all gwRVIS scores with BED-style genomic coordinates into a data frame
	all_gwrvis_bed_df = read_all_gwrvis_scores(all_gwrvis_bed_file)	

	
	get_most_intol_or_toler_sequences(all_gwrvis_bed_df)
	get_most_intol_or_toler_sequences(all_gwrvis_bed_df, tol_type='tolerant')
