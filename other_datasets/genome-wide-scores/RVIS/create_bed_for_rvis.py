import pandas as pd

rvis_df = pd.read_csv("RVIS_scores_per_gene.tsv", sep="\t", header=None)
rvis_df.columns = ['gene', 'rvis']
print(rvis_df.head())


ccds_df = pd.read_csv("CCDS.bed", sep="\t", header=None)
ccds_df.columns = ['chr', 'start', 'end', 'gene']
print(ccds_df.head())


full_df = rvis_df.merge(ccds_df, left_on='gene', right_on='gene', how='inner')
full_df = full_df[ ['chr', 'start', 'end', 'rvis'] ]
full_df['chr'] = 'chr' + full_df['chr']

# remove pathces
full_df = full_df[ ~full_df['chr'].str.contains("chrH") ]
print(full_df.shape)
print(full_df.head())

full_df.dropna(inplace=True)
print(full_df.shape)

full_df.to_csv("All_chromosomes.RVIS.bed", header=False, index=False, sep="\t")
