import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.svm import SVR
import sys
import os

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from custom_utils import create_out_dir



def is_outlier(points, thresh=3.5):
    if len(points.shape) == 1:
                points = points[:,None]
       
    median = np.median(points, axis=0)   
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    
    med_abs_deviation = np.median(diff)
    if len(points) == 1 and med_abs_deviation == 0.0:
        return np.array([False])
            
    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh




config_file = "../config.yaml"
out_dir = create_out_dir(config_file)    
base_dir = '../' + out_dir


full_gwrvis_df = pd.read_csv(base_dir + '/full_genome_out/Whole_genome_gwRVIS_per_class.csv')
full_gwrvis_df.head()




full_gwrvis_df.genomic_class.unique()


# In[6]:


#print('Full Df size:', df.shape)
#df = df.loc[ ~is_outlier(df.loc[:,'gwrvis']), :]
#print('Df size (no outlier):', df.shape)


# intolerantintol_class = 'ucne'
intolerant = full_gwrvis_df.loc[full_gwrvis_df.loc[:, 'genomic_class'] == intol_class].copy()

# beta:
intolerant = intolerant.loc[ intolerant['gwrvis'] <= 0, ]

del intolerant['genomic_class']
intolerant['gene_annot'] = 'intolerant'
print(intolerant.head())
intol_elems = intolerant.shape[0]
print('\nIntolerant class:', intol_class, ', size:', intol_elems)



# tolerant
toler_class = 'intergenic'
tolerant = full_gwrvis_df.loc[full_gwrvis_df.loc[:, 'genomic_class'] == toler_class].copy()
del tolerant['genomic_class']

# beta:
tolerant = tolerant.loc[ tolerant['gwrvis'] <= 0, ]


tolerant['gene_annot'] = 'tolerant'
tolerant = tolerant.sample(intol_elems)
print(tolerant.head())
print('\nTolerant class:', toler_class, ', size:', tolerant.shape[0])
tolerant.shape


# In[7]:


df = pd.concat([tolerant, intolerant], axis=0)
df['idx'] = df.index
df.head()
df.tail()


# In[8]:


df.plot(kind='scatter', x='idx', y='gwrvis')

fg = sns.FacetGrid(data=df, hue='gene_annot', aspect=2)
fg.map(plt.scatter, 'idx', 'gwrvis')


# In[9]:


df['gene_class'] = df.gene_annot.map({'tolerant':0, 'intolerant':1})
df.head()


# ### Logistic Regression

# In[10]:


import numpy as np
print(df.shape)
df = df[ ~np.isnan(df.gwrvis) ]
df.shape


# In[11]:


from sklearn.linear_model import LogisticRegression
from sklearn import svm

lr = LogisticRegression(C=1e9)
feature_cols = ['gwrvis']
X = df[feature_cols]
y = df['gene_class']

lr.fit(X, y)

# Support \/ector Regression:
#clf = SVR(C=1, epsilon=0.2)
# clf = svm.SVC(probability=True)
#clf.fit(X,y)

model = lr


# In[12]:


df['pred_gene_class_prob'] = model.predict_proba(X)[:, 1]


# In[ ]:


#df['pred_gene_class_prob']


# In[13]:


# plot the class predictions
plt.scatter(df.gwrvis, df.gene_class)
plt.plot(df.gwrvis, df.pred_gene_class_prob, color='red')
plt.xlabel('gwrvis')
plt.ylabel('gene class')

# df['pred_gene_class_prob'] = np.where(df['pred_gene_class_prob'] > 0.5, 1, 0)
# df['pred_gene_class_prob'].head()


# In[14]:


from sklearn.metrics import roc_curve, auc

fpr, tpr, thresholds = roc_curve(df['gene_class'], df['pred_gene_class_prob'])
roc_auc = auc(fpr, tpr)
print("Area under the ROC curve : %f" % roc_auc)


# ### Plot ROC curve

# In[15]:


plt.figure(figsize=(4, 4), dpi=160)
lw = 2
plt.plot(fpr, tpr, color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.show()

