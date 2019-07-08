import pandas as pd
import sys

pop_ratios = {'European': 0.4846, 'African': 0.2819, 'Asian': 0.0523, 'Other': 0.1288}

df = pd.read_csv('heptamer_mutability_rates.txt', sep='\t', header=0)
print(df.head())


agg_df = df.groupby(['ref_seq'])['African', 'Asian', 'European'].agg('sum')
agg_df['weighted_sum'] = pop_ratios['European'] * agg_df['European'] + \
			   pop_ratios['African'] * agg_df['African'] + \
			   pop_ratios['Asian'] * agg_df['Asian'] + \
			   pop_ratios['Other'] * agg_df['European']

agg_df['sum'] = (agg_df['European'] + agg_df['African'] + \
		agg_df['Asian'] + agg_df['European']) / 4
print(agg_df.head())


# reverse complement motifs
rev_agg_df =  df.groupby(['rev_complement_ref'])['African', 'Asian', 'European'].agg('sum')
print(rev_agg_df.head())
rev_agg_df['weighted_sum'] = pop_ratios['European'] * rev_agg_df['European'] + \
			   pop_ratios['African'] * rev_agg_df['African'] + \
			   pop_ratios['Asian'] * rev_agg_df['Asian'] + \
			   pop_ratios['Other'] * rev_agg_df['European']

rev_agg_df['sum'] = (rev_agg_df['European'] + rev_agg_df['African'] + \
		rev_agg_df['Asian'] + rev_agg_df['European']) / 4
print(rev_agg_df.head())


# Quick validation:
#print(agg_df.shape)
#print(rev_agg_df.shape)
#print(rev_agg_df[ rev_agg_df.index == 'TTTTTTT' ]) # verified to be the same as agg_df for 'AAAAAAA'


total_agg_df = pd.concat([agg_df, rev_agg_df], axis=0)
print(total_agg_df.head())
print(total_agg_df.shape)




sums_df = total_agg_df[['weighted_sum', 'sum']]
sums_df.index.name = 'ref_seq'
print(sums_df.head())

sums_df.to_csv('heptamer_mutability_rates.processed_sums.txt')

