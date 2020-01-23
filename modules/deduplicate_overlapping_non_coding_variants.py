import sys
import pandas as pd
import numpy as np


input_file = sys.argv[1]

df = pd.read_csv(input_file, sep='\t')
print(df.shape)


df['win_index'] = df['start'] / 3000
df['win_index'] = df['win_index'].apply(np.floor)

df.drop_duplicates(subset=['win_index'], inplace=True)
df.drop(columns=['win_index'], axis=1, inplace=True)
print(df.shape)

df.to_csv(input_file, index=False, header=True, sep='\t')
