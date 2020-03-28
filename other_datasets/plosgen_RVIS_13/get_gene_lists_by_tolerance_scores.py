import pandas as pd
from sys import argv

input_file = argv[1]

df = pd.read_csv(input_file)
df.columns = ['gene', 'rvis', 'percentile']
print(df.head())
print(df.shape)


# 25% most intolerant genes
most_intolerant_genes_df = df.loc[ df['percentile'] <= 25, :]
print(most_intolerant_genes_df.head())
print(most_intolerant_genes_df.shape)
most_intolerant_genes_df['gene'].to_csv('most_intolerant_genes_25percentile.csv', index=False)

# 75% most tolerant genes
most_tolerant_genes_df = df.loc[ df['percentile'] > 75, :]
print(most_tolerant_genes_df.head())
print(most_tolerant_genes_df.shape)
most_tolerant_genes_df['gene'].to_csv('most_tolerant_genes_75percentile.csv', index=False)
