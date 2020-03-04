import pandas as pd

df = pd.read_csv('jarvis-phastcons-output/mean_scores_by_window.bed', sep='\t', header=None)
df.columns = ['chr', 'star', 'end', 'gwrvis', 'jarvis', 'phastcons']
print(df.head())
