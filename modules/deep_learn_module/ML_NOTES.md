___ RESULTS with conv[32@50]_mp[5@15]_conv[64@25]_mp[5@15]_fcc[128], batch_size=1024, learning_rate=0.001 (in adam) ___

>>> Non-random:
Epoch 00062: val_loss improved from 0.54272 to 0.53991, saving model to gwrvis_best_model.hdf5
1188/1188 [==============================] - 14s 12ms/sample - loss: 0.4666 - cosine_similarity: 0.8241 - val_loss: 0.5399 - val_cosine_similarity: 0.7895

> Confusion matrix:
 [[127  26]
 [ 54  91]]
TN: 127
FP: 26
FN: 54
TP: 91
> ROC AUC: 0.7288257831868379
> Accuracy: 0.7315436241610739

> Sensitivity: 0.6275862068965518
> Precision: 0.7777777777777778
> Specificity: 0.8300653594771242


>>> Random:
Epoch 00005: val_loss improved from 0.69309 to 0.69300, saving model to gwrvis_best_model.hdf5
1240/1240 [==============================] - 15s 12ms/sample - loss: 0.6933 - cosine_similarity: 0.7070 - val_loss: 0.6930 - val_cosine_similarity: 0.7072


> Confusion matrix:
 [[158   0]
 [152   0]]
TN: 158
FP: 0
FN: 152
TP: 0
train_nn_model.py:140: RuntimeWarning: invalid value encountered in long_scalars
  precision = TP / (TP + FP)
> ROC AUC: 0.5
> Accuracy: 0.5096774193548387

> Sensitivity: 0.0
> Precision: nan
> Specificity: 1.0


<<<<< OR >>>>>
Epoch 00007: val_loss improved from 0.69327 to 0.69315, saving model to gwrvis_best_model.hdf5
1188/1188 [==============================] - 14s 12ms/sample - loss: 0.6931 - cosine_similarity: 0.7071 - val_loss: 0.6932 - val_cosine_similarity: 0.7071

> Confusion matrix:
 [[  0 160]
 [  0 137]]
TN: 0
FP: 160
FN: 0
TP: 137
> ROC AUC: 0.5
> Accuracy: 0.4612794612794613

> Sensitivity: 1.0
> Precision: 0.4612794612794613
> Specificity: 0.0