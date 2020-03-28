import pandas as pd
from sys import argv
import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '../../bin'))
from custom_utils import splitDataFrameList


df = pd.read_table('ucsc_targetscan_mirna_targets_3utrs.bed', header=None)
df.columns = ['chr', 'start', 'end', 'annot']
df[['gene', 'mirnas', 'aux']] = df['annot'].str.split(':', expand=True) 
del df['aux']
del df['annot']

print(df.head())
df['mirnas'] = df['mirnas'].str.replace('/', '/miR-')
print(df.head(50))
sys.exit()
df = splitDataFrameList(df, 4, '/')
df.columns.values[-1] = 'miR'
del df['mirnas']
df = df[['chr', 'start', 'end', 'gene', 'miR']]
print(df.head())
