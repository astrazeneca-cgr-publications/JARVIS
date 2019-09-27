import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score


def plot_roc_curve(test_flat, preds_flat, make_plot=False):

	fpr, tpr, _ = roc_curve(test_flat, preds_flat)
	roc_auc = roc_auc_score(test_flat, preds_flat)

	if make_plot:
	    f = plt.figure(figsize=(6, 6))
	    _ = plt.plot(fpr, tpr, label='ROC curve (area = %0.3f)' % roc_auc)
	    _ = plt.plot([0, 1], [0, 1], '--', linewidth=0.5)  # random predictions curve

	    _ = plt.xlim([0.0, 1.0])
	    _ = plt.ylim([0.0, 1.0])
	    _ = plt.title('\nROC (area = %0.3f)' % roc_auc)
	    _ = plt.xlabel('False Positive Rate (1 — Specificity)')
	    _ = plt.ylabel('True Positive Rate (Sensitivity)')
	    plt.grid(True)
	    plt.show()

	    f.savefig("ROC_curve.pdf", bbox_inches='tight')

	return fpr, tpr
