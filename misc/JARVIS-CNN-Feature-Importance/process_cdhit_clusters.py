import numpy as np
import os
import sys
import pickle


cdhit_clusters = sys.argv[1]   #cdhit_seq_clusters.id80.W7.fa.clstr

cdhit_fasta = 'JARVIS-CNN_pattern_seqs.fasta'  #cdhit_out_identifier + '.fa'


id_to_seq = {}

tmp_id = ''
tmp_seq = ''
with open(cdhit_fasta) as fh:
	for line in fh:
		line = line.rstrip()
	
		if line.startswith('>'):
			tmp_id = line.replace('>', '')
		else:
			tmp_seq = line
			id_to_seq[tmp_id] = tmp_seq




# ===================================
with open('JARVIS-filter_act_per_seq.pkl', 'rb') as handle:
    filter_act_per_seq = pickle.load(handle)


def most_common(lst):
	lst = list(lst)
	return max(set(lst), key=lst.count)


def get_consensus_sequence(seqs_list):

	full_arr = []

	# return sole sequence for singleton clusters
	if len(seqs_list) == 1:
		return seqs_list[0]

	for seq in seqs_list:
		if len(full_arr) > 0:
			full_arr = np.vstack((full_arr, list(seq)))	
		else:
			full_arr = np.array(list(seq))
	#print(full_arr)	

	cons_seq = ''
	for j in range(full_arr.shape[1]):
		cons_seq += most_common(full_arr[:, j])

	#print(cons_seq)
	return cons_seq
	



clust_out_dir = 'clusters-output'
if not os.path.exists(clust_out_dir):
	os.makedirs(clust_out_dir)



seqlists_per_clust = {}
act_sum_per_clust = {}
cons_seq_per_clust = {}

clust_id = None
seqs_in_clust = []
act_per_clust = []


with open(cdhit_clusters) as fh:
	for line in fh:
		line = line.rstrip()

		if line.startswith('>'):

			if len(seqs_in_clust) > 0:
				seqlists_per_clust[clust_id] = seqs_in_clust
				act_sum_per_clust[clust_id] = np.sum(act_per_clust)

				#print(clust_id)
				#print(seqs_in_clust)
				#print(act_per_clust)
				
				cons_seq = get_consensus_sequence(seqs_in_clust)
				cons_seq_per_clust[clust_id] = cons_seq
				#print(cons_seq)
			
				# reset tmp variables
				seqs_in_clust = []
				act_per_clust = []


			clust_id = (line.replace('>', '')).replace(' ', '_')

		else:
			# get sequence id
			_, vals = line.split('\t')
			#print(vals)
			id_vals = vals.split(' ')
			seq_id = id_vals[1]
			seq_id = seq_id.replace('...', '').replace('>', '')
			#print(seq_id)

			cur_seq = id_to_seq[seq_id]
			seqs_in_clust.append(cur_seq)
			act_per_clust.append( filter_act_per_seq[cur_seq] )

			
print(len(seqlists_per_clust))
print(len(act_sum_per_clust))
print(len(cons_seq_per_clust))



for clust_id in seqlists_per_clust.keys():

	cur_seq_out_dir = clust_out_dir + '/' + clust_id
	if not os.path.exists(cur_seq_out_dir):
		os.makedirs(cur_seq_out_dir)

	with open(cur_seq_out_dir + '/sequences.txt', 'w') as fh:
		for seq in seqlists_per_clust[clust_id]:
			fh.write(seq + '\n')

	with open(cur_seq_out_dir + '/summary.txt', 'w') as fh:
		fh.write('consensus_seq\t' + cons_seq_per_clust[clust_id] + '\n')
		fh.write('filter_activation_sum\t' + str(act_sum_per_clust[clust_id]) + '\n')



# Focus on Clusters with highest filter activation sums
print('>> Preparing TomTom input file with consensus sequences...')
seq_out_fh = open('top_jarvis_cons_seqs.fa', 'w')
clust_id_out_fh = open('top_jarvis_cons_seqs.clust_id', 'w')
act_sum_out_fh = open('top_jarvis_cons_seqs.act_sum', 'w')

top_clusters = 100

cnt = 1
for clust_id, act_sum in sorted(act_sum_per_clust.items(), key=lambda x: x[1], reverse=True):
	print(clust_id)
	print(act_sum)
	cons_seq = cons_seq_per_clust[clust_id]
	print(cons_seq)
	
	seq_out_fh.write(cons_seq + '\n')
	clust_id_out_fh.write(clust_id + '\n')
	act_sum_out_fh.write(str(act_sum) + '\n')

	if cnt == 100:
		break
	cnt += 1

# Prepare TomTom input file with cluster __consensus__ sequences
os.system("""awk ' {print;} NR % 1 == 0 { print ""; }' top_jarvis_cons_seqs.fa > top_jarvis_cons_seqs.tomtom_input""")
print('<< Complete - Output file: top_jarvis_cons_seqs.tomtom_input')




# Prepare TomTom input file with all sequences from each cluster
print('>> Preparing TomTom input file with all sequences per cluster...')
tomtom_out_fh = open('top_jarvis_cons_seqs.tomtom_input.all_seqs', 'w')

cnt = 1
for clust_id, act_sum in sorted(act_sum_per_clust.items(), key=lambda x: x[1], reverse=True):
	print(clust_id)
	print(act_sum)

	for seq in seqlists_per_clust[clust_id]:
		tomtom_out_fh.write(seq + '\n')
	tomtom_out_fh.write('\n')
	
	if cnt == 100:
		break
	cnt += 1

print('<< Complete - Output file: top_jarvis_cons_seqs.tomtom_input.all_seqs')

#TODO:
# Store top seqs in fasta file and run with deepbind
# - Select patterns that return the top hits annotations
# - Compile web logos for top hits (based on activation-sums or based on deepbind results)
