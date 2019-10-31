import pandas as pd

mean_depths_file = 'gnomad.r3.0.coverage_means.txt'
#mean_depths_file = 'tmp/means.sample'

df = pd.read_csv(mean_depths_file, header=None, low_memory=False)

out_fh = open('mean_depth_quantiles.tsv', 'w')
out_fh.write('quantile\tvalue\n')


for q in [0.1, 0.2, 0.3, 0.4]:

	depth_q = df.quantile(q).values[0]
	out_fh.write(str(q*100) + '%\t' + str(depth_q) + '\n')

out_fh.close()
	

