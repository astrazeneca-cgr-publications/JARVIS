from sys import argv
import pandas as pd

input_file = argv[1] # doi_ejhg.2008.111_haploinsufficient_genes.csv.tmp

df = pd.read_table(input_file, sep=',')
df.columns = ['gene', 'chr', 'start', 'end', 'omim']
df = df[['chr', 'start', 'end', 'gene', 'omim']]
df['chr'] = 'chr' + df['chr']
del df['omim']
print(df.head())

df.to_csv('../hg19/bed/' + input_file.replace('csv.tmp', 'bed'), header=None, index=False, sep='\t')
