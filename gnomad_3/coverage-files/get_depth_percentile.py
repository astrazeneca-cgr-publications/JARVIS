import pandas as pd

mean_depths_file = 'gnomad.r3.0.coverage_means.txt.gz'
#mean_depths_file = 'tmp/means.sample'

print("Reading data frame ...")
df = pd.read_csv(mean_depths_file, header=None, low_memory=False)
print("... Done")


print("Opening file to store quantile results ...")
out_fh = open('mean_depth_quantiles.tsv', 'w')
out_fh.write('quantile\tvalue\n')


print("Calulating quantiles ...")
for q in [0.1, 0.2, 0.3, 0.4]:

	depth_q = df.quantile(q).values[0]
	out_fh.write(str(q*100) + '%\t' + str(depth_q) + '\n')
print("... Done")

out_fh.close()
	

