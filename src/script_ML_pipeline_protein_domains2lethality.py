# -*- coding: utf-8 -*-
# %% [markdown]
# """
# Created on Fri Oct  2 09:38:25 2020
#
# @author: linigodelacruz
# """
#
# # Notebook showing the Machine learning pipeline to predict SL pairs from protein domains "homology"
# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict 
import seaborn as sns
import matplotlib.cm as cm
import scipy as scipy


import os

from python_modules.module_ML_protein_domains_SL import sample_protein_pairs
from python_modules.module_ML_protein_domains_SL import remove_empty_domains
from python_modules.module_ML_protein_domains_SL import feature_building
import time 

t = time.process_time()

# %% [markdown]
# ## Retrieving local data 

# %% Importing datasets

script_dir = os.path.dirname('__file__') #<-- absolute dir the script is in
rel_path_SL = "datasets/data-synthetic-lethals.xlsx"
rel_path_nSL="datasets/data-positive-genetic.xlsx"
rel_path_domains="datasets/proteins-domains-from-Pfam.xlsx"

data_sl=pd.read_excel(os.path.join(script_dir, rel_path_SL),header=0)
data_domains=pd.read_excel(os.path.join(script_dir, rel_path_domains),header=0,index_col='Unnamed: 0')
data_domains=data_domains.dropna()
data_nonsl=pd.read_excel(os.path.join(script_dir, rel_path_nSL),header=0)

domain_id_list=data_domains['domain-name']
#%% selecting size of the protein pairs to analyze

protein_a_list,protein_a_list_non,protein_b_list,protein_b_list_non=sample_protein_pairs(data_domains=data_domains,
                                                                                         data_sl=data_sl,data_nonsl=data_nonsl,sample_size=100)

# %% printing result
print('We are going to analyze',len((protein_a_list)) ,'randomly selected protein pairs, out of',len(data_sl),'SL protein pairs')
print('We are going to analyze',len((protein_a_list_non)) ,'randomly selected protein pairs, out of',len(data_nonsl),'positive protein pairs')

# %% [markdown]
# ### Removing empty protein domains

# %% Postprocesing 1: Removing empty protein domains

protein_a_list_new,protein_b_list_new=remove_empty_domains(protein_a_list,protein_b_list)
protein_a_list_non_new,protein_b_list_non_new=remove_empty_domains(protein_a_list_non,protein_b_list_non)

# %% printing result

print('The empty domain in the SL were:', len(protein_a_list)-len(protein_a_list_new), 'out of', len(protein_a_list),'domains')
print('The empty domain in the nSL were:', len(protein_a_list_non)-len(protein_a_list_non_new), 'out of', len(protein_a_list_non),'domains')
# %% [markdown]
# ### Feature building 

# %% Feature engineering
protein_feat_true_pd=feature_building(protein_b_list_new=protein_b_list_new,protein_a_list_new=protein_a_list_new,domain_id_list=domain_id_list)
protein_feat_non_true_pd=feature_building(protein_b_list_new=protein_b_list_non_new,protein_a_list_new=protein_a_list_non_new,domain_id_list=domain_id_list)


# %% [markdown]
# ### Visualizing the features per type 

# %% printing results

index_2_true=protein_feat_true_pd.where(protein_feat_true_pd==2)
index_2_true_count=index_2_true.count(axis=1).sum()

index_1_true=protein_feat_true_pd.where(protein_feat_true_pd==1)
index_1_true_count=index_1_true.count(axis=1).sum()

index_2_nontrue=protein_feat_non_true_pd.where(protein_feat_non_true_pd==2)
index_2_nontrue_count=index_2_nontrue.count(axis=1).sum()

index_1_nontrue=protein_feat_non_true_pd.where(protein_feat_non_true_pd==1)
index_1_nontrue_count=index_1_nontrue.count(axis=1).sum()


# %% Vizualization

colors = ['#00F28E','#F20064']
plt.bar(['fraction of 2 in the nSL','fraction of 1 in the nSL'],[index_2_nontrue_count/(len(protein_feat_non_true_pd.index)*len(protein_feat_non_true_pd.columns)),index_1_nontrue_count/(len(protein_feat_non_true_pd.index)*len(protein_feat_non_true_pd.columns))],alpha=0.6,color=[colors[0],'black']), 

plt.bar(['fraction of 2 in SL ','fraction of 1 in SL'],[index_2_true_count/(len(protein_feat_true_pd.index)*len(protein_feat_true_pd.columns)),index_1_true_count/(len(protein_feat_true_pd.index)*len(protein_feat_true_pd.columns))],alpha=0.6,color=[colors[1],'black'])

plt.ylabel('Fraction from the population')
plt.yscale('log')
plt.xticks(rotation=40)

# %% Adding the labels (response variables) to each dataset

protein_feat_true_pd['lethality']=np.ones(shape=(len(protein_a_list_new)))
protein_feat_non_true_pd['lethality']=np.zeros(shape=(len(protein_a_list_non_new)))

# %% Joining both datasets

feature_post=pd.concat([protein_feat_true_pd,protein_feat_non_true_pd],axis=0)
feature_post=feature_post.set_index(np.arange(0,len(protein_a_list_new)+len(protein_a_list_non_new)))
print('The number of features are:',feature_post.shape[1])
print('The number of samples are:',feature_post.shape[0])

# %% Postprocessing and exploration of the feature matrix of both datasets

mean=feature_post.T.describe().loc['mean']
std=feature_post.T.describe().loc['std']
lethality=feature_post['lethality']

corr_keys=pd.concat([mean,std,lethality],axis=1)


# %% Viz of the stats
fig, axs = plt.subplots(ncols=2, figsize=(10,5))
a=sns.violinplot(x="lethality", y="mean", data=corr_keys,ax=axs[0],palette='Paired')
a.set_title('How the mean varies with Lethality')
b=sns.violinplot(x="lethality", y="std", data=corr_keys,ax=axs[1],palette='Paired')
b.set_title('How the std varies with Lethality')

pair=sns.pairplot(corr_keys,hue='lethality',diag_kind='kde',kind='reg',palette='Paired')
pair.fig.suptitle('Pairplot to see data dependencies with Lethality',y=1.08)

plt.savefig('Lethality_correlation_with_mean_and_std.png',dpi=60)
# %% P- values and correlation

a=scipy.stats.pearsonr(corr_keys['mean'],corr_keys['lethality'])
p_value_corr=defaultdict(dict)

columns=['mean','std']
for i in columns:
    
    tmp=scipy.stats.pearsonr(corr_keys[i],corr_keys['lethality'])
    p_value_corr[i]['corr with lethality']=tmp[0]
    p_value_corr[i]['p-value']=tmp[1]

p_value_corr_pd=pd.DataFrame(p_value_corr)

# %% Viz of correlations

corr = corr_keys.corr()
fig, axs = plt.subplots(ncols=2, figsize=(10,5))
plt.subplots_adjust(wspace=1)
sns.heatmap(corr, vmax=1,vmin=-1 ,square=True,cmap=cm.PRGn,cbar_kws={'label':'Pearson corr'},ax=axs[0])
axs[1].scatter(x=p_value_corr_pd.loc['corr with lethality','mean'], y=p_value_corr_pd.loc['p-value','mean'],color='b',label='mean',s=100)
axs[1].scatter(x=p_value_corr_pd.loc['corr with lethality','std'], y=p_value_corr_pd.loc['p-value','std'],color='r',label='std',s=100)

axs[1].set_ylabel('p-value')
axs[1].set_xlabel('correlation with lethality')
axs[1].legend()

# %% Separate features from labels to set up the data from the ML workflow

X, y = feature_post.drop(columns=["lethality"]), feature_post["lethality"]

# %% Separating training and testing set
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test =  train_test_split(X,y,test_size = 0.3, random_state= 0)

print ('Train set:', X_train.shape,  y_train.shape)
print ('Test set:', X_test.shape,  y_test.shape)

# %% Choosing the best SVM model

# from sklearn.model_selection import GridSearchCV
# from sklearn.svm import SVC
# parameters = [{'C': [1, 10, 100], 'kernel': ['rbf'], 'gamma': ['auto','scale']}]
# search = GridSearchCV(SVC(), parameters, n_jobs=-1, verbose=1)
# search.fit(X_train, y_train)

# best_parameters = search.best_estimator_
# print(best_parameters)

# %% Training with the best model
print('Training the model...')
from sklearn import svm

clf = svm.SVC(C=10, break_ties=False, cache_size=200, class_weight=None, coef0=0.0,
    decision_function_shape='ovr', degree=3, gamma='scale', kernel='rbf',
    max_iter=-1, probability=False, random_state=None, shrinking=True,
    tol=0.001, verbose=False).fit(X_train, y_train)
print('The  score of the model is =',clf.score(X_test, y_test))

# %% Making predictions
print('Making predictions and evaluating the model...')
from sklearn import metrics
from sklearn.metrics import log_loss
from sklearn.metrics import jaccard_score
from sklearn.metrics import classification_report

y_pred = clf.predict(X_test)
print('Train set Accuracy: ', metrics.accuracy_score(y_train, clf.predict(X_train)))
print('The mean squared error is =',metrics.mean_squared_error(y_test,y_pred))
print('Test set Accuracy: ', metrics.accuracy_score(y_test, y_pred))
print('The Jaccard index is =', jaccard_score(y_test, y_pred))
# Jaccard similarity coefficient, defined as the size of the intersection divided by the size of the union of two label sets. The closer to 1 the better the classifier 
print('The log-loss is =',log_loss(y_test,y_pred))
# how far each prediction is from the actual label, it is like a distance measure from the predicted to the actual , the classifer with lower log loss have better accuracy
print('The f1-score is =',metrics.f1_score(y_test,y_pred))
# The F1 score can be interpreted as a weighted average of the precision and recall, where an F1 score reaches its best value at 1 and worst score at 0. The relative contribution of precision and recall to the F1 score are equal.

# Model Precision: what percentage of positive tuples are labeled as such?
print("Precision:",metrics.precision_score(y_test, y_pred))

# Model Recall: what percentage of positive tuples are labelled as such?
print("Recall:",metrics.recall_score(y_test, y_pred))

print(classification_report(y_test, y_pred, target_names=['NonSl','SL']))
#%% writing a document

# Report training set score
train_score = metrics.accuracy_score(y_train, clf.predict(X_train)) * 100
# Report test set score
test_score = metrics.accuracy_score(y_test, y_pred) * 100
# Report precision
precision=metrics.precision_score(y_test, y_pred) *100

# report Recall
recall=metrics.recall_score(y_test, y_pred) *100

# Write scores to a file
with open("metrics.txt", 'w') as outfile:
        outfile.write("Training accuracy: %2.1f%%\n" % train_score)
        outfile.write("Test accuracy : %2.1f%%\n" % test_score)
        outfile.write("Precision : %2.1f%%\n" % precision)
        outfile.write("Recall : %2.1f%%\n" % recall)
        

# %% ROC curves

fig, axs = plt.subplots(ncols=2, figsize=(10,5))
import sklearn.metrics as metrics
scores=clf.decision_function(X_test)

fpr, tpr, thresholds = metrics.roc_curve(y_test, scores)
area=metrics.auc(fpr,tpr)
axs[0].plot(fpr,tpr,color='darkorange',label='SVM model (area = %0.2f)' % area)
axs[0].plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--',label='Random prediction')
axs[0].set_xlim([0.0, 1.0])
axs[0].set_ylim([0.0, 1.05])
axs[0].set_xlabel('False Positive Rate')
axs[0].set_ylabel('True Positive Rate')
axs[0].set_title('Receiver operating characteristic example')
axs[0].legend(loc="lower right")

# %% Precision vs Recall curves

precision, recall, thresholds = metrics.precision_recall_curve(y_test, scores)
average_precision = metrics.average_precision_score(y_test, scores)
axs[1].plot(precision,recall,color='blue',label='SVM-model')

axs[1].plot([0.5, 1], [1, 0], color='navy', lw=2, linestyle='--',label='Random prediction')

axs[1].set_xlabel('Recall')
axs[1].set_ylabel('Precision')
axs[1].set_title('2-class Precision-Recall curve: '
                   'AP={0:0.2f}'.format(average_precision))
axs[1].legend()


# %% Confusion matrix
class_names=[1,2,3]
fig, ax = plt.subplots()
from sklearn.metrics import confusion_matrix
import sklearn.metrics as metrics

cm = confusion_matrix(y_test, y_pred,normalize="true")

class_names=['SL', 'nSL']

tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)

sns.heatmap(pd.DataFrame(cm), annot=True, cmap="Blues" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('Confusion matrix', y=1.1)
plt.ylabel('Actual label')
plt.xlabel('Predicted label')
plt.savefig("confusion_matrix.png",dpi=60)

# %% Evaluation of the classifier in terms of overfitting : Cross Validation!

print('Cross validation to test overfitting..')

from sklearn.model_selection import StratifiedKFold
import time
import sklearn.metrics as metrics

from sklearn.model_selection import cross_validate


n_samples = X.shape[0]
cv=StratifiedKFold(n_splits=5)


cv_results = cross_validate(clf, X, y, cv=cv)
elapsed_time = time.process_time() - t

# %% Viz of the variation of the test error per fold . If the variation is high , the classifier may be proned to overfitting.

fig, axs = plt.subplots(ncols=1, figsize=(3,3))
sorted(cv_results.keys())

# plt.scatter(['test-1','test-2','test-3','test-4','test-5'],cv_results['test_score'],s=60,alpha=0.7,color='blue')


plt.errorbar(x=['all tests'],y=np.mean(cv_results['test_score']) , yerr=np.std(cv_results['test_score']),capsize=10,fmt='-o',ecolor='black', capthick=2)
plt.ylim(0,1)
plt.title('5-fold crossvalidation result')
plt.ylabel('Accuracy')
plt.savefig('5-fold-crossvalidation.png',dpi=60)

print('The elapsed time was',elapsed_time/60,'minutes')
