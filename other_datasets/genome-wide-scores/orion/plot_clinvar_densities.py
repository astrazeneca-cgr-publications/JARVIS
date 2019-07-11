import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys

pathogenic_scores_file = sys.argv[1]
benign_scores_file = sys.argv[2]
include_classes = sys.argv[3]

#try:
pathogenic = pd.read_csv(pathogenic_scores_file, header=None, sep='\t')
pathogenic = pathogenic.dropna().iloc[:, 3].tolist()
benign = pd.read_csv(benign_scores_file, header=None, sep='\t')
benign = benign.dropna().iloc[:, 3].tolist()
#except Exception as e:
#	print("No columns to parse from file for genomic class:", include_classes)

print(type(pathogenic))
print(type(benign))

# Plot
try:
	fig, ax = plt.subplots(figsize=(10, 10))

	sns.distplot(tuple(pathogenic), hist=False, kde=True, label='pathogenic (' + str(len(pathogenic)) + ')')
	sns.distplot(tuple(benign), hist=False, kde=True, label='benign (' + str(len(benign)) + ')')
	plt.title(include_classes)

	fig.savefig('ClinVar_pathogenic_vs_benign.density.' + include_classes + '.pdf')
except:
	print('Insufficient data points for genomic class:', include_classes)
