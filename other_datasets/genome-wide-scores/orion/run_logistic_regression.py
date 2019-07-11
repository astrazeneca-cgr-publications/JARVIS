import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
import sys
import os



def run_logit_regression(df):

	X = df['score'].values.reshape(-1, 1)
	y = df['pathogenic'].values

	# logistic regression
	model = LogisticRegression(C=1e9, solver='lbfgs')
	model.fit(X, y)

	return model, X, y



def plot_roc_curve(model, df, X):

	df['pred_gene_class_prob'] = model.predict_proba(X)[:, 1]

	fpr, tpr, thresholds = roc_curve(df['pathogenic'], df['pred_gene_class_prob'])
	roc_auc = round(auc(fpr, tpr), 2)
	print("Area under the ROC curve : %f" % roc_auc)

	# Plot ROC curve
	fig, ax = plt.subplots(figsize=(10, 10))
	plt.plot(fpr, tpr, color='darkorange',
		 lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
	plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('Receiver operating characteristic')
	plt.legend(loc="lower right", fontsize=14)
	plt.show()

	print('\n-------------------\n', genomic_class, 'AUC:', str(roc_auc))


	pdf_filename = './Logistic_Regression_ROC.AUC_' + str(roc_auc) + '.' + genomic_class + '.pdf'

	fig.savefig(pdf_filename, bbox_inches='tight')




if __name__ == '__main__':

	pathogenic = sys.argv[1]	
	benign = sys.argv[2]	
	genomic_class = sys.argv[3]
	

	try:
		pathogenic_df = pd.DataFrame(pd.read_csv(pathogenic, sep='\t', header=None)[3])
		pathogenic_df.columns = ['score']
		benign_df = pd.DataFrame(pd.read_csv(benign, sep='\t', header=None)[3])
		benign_df.columns = ['score']
	except:
		sys.exit('Insufficient data for class:', genomic_class)
	
	pathogenic_df['pathogenic'] = 1
	benign_df['pathogenic'] = 0
	
	pathogenic_df.dropna(inplace=True)
	benign_df.dropna(inplace=True)


	df = pd.concat([pathogenic_df, benign_df], axis=0)



	# ### Logistic Regression
	#df = df[ ~np.isnan(df.gwrvis) ]
	model, X, y = run_logit_regression(df)

	plot_roc_curve(model, df, X)
