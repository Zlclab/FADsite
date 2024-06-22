import sys
import pickle
import pandas as pd
import numpy as np
from sklearn import metrics
from sklearn.metrics import matthews_corrcoef
from catboost import CatBoostClassifier

in_file = sys.argv[1]           
label_ = []#########A list for each protein's label
with open('..\\dataset\\'+in_file+'.txt', 'r') as f:
    lines = f.readlines()
for i in range(len(lines)):
    if i % 3 == 2:
        label_.append(lines[i][:-1])
label = [int(i) for i in ''.join(label_)]######### All protein's label
with open('FADsite_seq.pkl', 'rb') as f:
    model = pickle.load(f)
    
ori = pd.read_excel(in_file + '.xlsx')
test_X = ori.iloc[:,:]
test_label = label 
resample_pred = model.predict(np.array(test_X))
y_score = model.predict_proba(np.array(test_X))[:,1]
fpr,tpr,threshold = metrics.roc_curve(test_label, y_score)
roc_auc = metrics.auc(fpr,tpr)
print("The results for FADsite_seq on " +in_file)
print('Pre:',metrics.precision_score(test_label,resample_pred))
print('f1:',metrics.f1_score(test_label,resample_pred))
print('Acc:',metrics.accuracy_score(test_label,resample_pred))
print('AUC:',roc_auc)
print('MCC:',matthews_corrcoef(test_label,resample_pred))    