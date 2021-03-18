# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 10:25:17 2021
This script will be used to do basic data nalaysis for the satay libraries according 
what is presenting in (Michel et al., 2017)
@author: linigodelacruz
"""
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from collections import defaultdict
import os
import seaborn as sns
import scipy 
#%% Import of dataframes output from the SATAY pipeline

names_libraries={'wt':'WT-merged-all.xlsx','dnrp1':'dnrp1-merged-all.xlsx','wt_agnes':'WT-merged-Agnes.xlsx','wt_strong_alig':'WT-merged-astringent-alignment.xlsx'}
data_library=[]
for i in names_libraries.keys():
 data_library.append(pd.read_excel('datasets/'+names_libraries[i],index_col='Unnamed: 0'))
#%% Creating a big dataframe of the libraries

data_library_pd=pd.concat(data_library,keys=names_libraries.keys(),sort=True)
data_library_pd.fillna(0,inplace=True)
#%% Analysis of transposition sites

def frequency_transposons(data,names_libraries):
    freq=[]
    for i in names_libraries.keys():
        freq.append(data.loc[i]['Nbasepairs'].sum()/data.loc[i]['Ninsertions'].sum())
        
    return freq

def reads_per_transposon(data,names_libraries):
    readspertr=[]
    for i in names_libraries.keys():
        readspertr.append(data.loc[i]['Nreads'].median()/data.loc[i]['Ninsertions'].median())
        
    return readspertr
    
def transposon_density(data,names_libraries):
    density=[]
    for i in names_libraries.keys():
        
        density.append(data.loc[i]['Ninsertions']/data.loc[i]['Nbasepairs'])        
    return density
def median_feature(data,names_libraries,feature):
    median_feature=[]
    for i in names_libraries.keys():
        
        median_feature.append(data.loc[i][feature].median())        
    return median_feature
    
def median_feature_essentials(data,names_libraries,feature):
    median_feature=[]
    for i in names_libraries.keys():
        
        median_feature.append(data.loc[i][data.loc[i]['Essentiality']==1][feature].median())        
    return median_feature
    
def median_feature_nonessentials(data,names_libraries,feature):
    median_feature=[]
    for i in names_libraries.keys():
        
        median_feature.append(data.loc[i][data.loc[i]['Essentiality']==0][feature].median())        
    return median_feature
#%%
freq=frequency_transposons(data_library_pd,names_libraries)
reads_per_transposon=reads_per_transposon(data_library_pd,names_libraries)
tr_density=transposon_density(data_library_pd,names_libraries)
median_insertions=median_feature(data_library_pd,names_libraries,'Ninsertions')
median_insert_essentials=median_feature_essentials(data_library_pd,names_libraries,'Ninsertions')
median_insert_nonessentials=median_feature_nonessentials(data_library_pd,names_libraries,'Ninsertions')
#%%
analysis_libraries=defaultdict(dict)

j=0
for i in names_libraries.keys():
    
    analysis_libraries[i]['one-transposon-per-bp']=freq[j]
    analysis_libraries[i]['median-reads-per-transposons']=reads_per_transposon[j]
    #analysis_libraries[i]['tr-density']=tr_density[j]
    analysis_libraries[i]['median-tr']=median_insertions[j]
    analysis_libraries[i]['median-tr-essentials']=median_insert_essentials[j]
    analysis_libraries[i]['median-tr-non-essentials']=median_insert_nonessentials[j]
    j=j+1
    
analysis_libraries_pd=pd.DataFrame(analysis_libraries)

#%% Defining the dataframes per type 
data_wt=data_library_pd.loc['wt'].copy()
data_wt_agnes=data_library_pd.loc['wt_agnes'].copy()
data_wt_greg2=data_library_pd.loc['wt_strong_alig'].copy()


#%% Transposon density vs genes highlighting the centromere position


data_wt['tr-density']=data_wt['Ninsertions']/data_wt['Nbasepairs']
data_wt_agnes['tr-density']=data_wt_agnes['Ninsertions']/data_wt_agnes['Nbasepairs']
data_wt_greg2['tr-density']=data_wt_greg2['Ninsertions']/data_wt_greg2['Nbasepairs']
#%% Plot transposon density (fig 1B Benoit)

fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)
ax.plot(data_wt['tr-density'],alpha=0.5,color='b')
ax.set_ylabel('transposond density: tn/bp')
ax.set_xlabel('genes')
## annotated centromeres
for i in np.arange(0,len(data_wt)):
    
    if data_wt.loc[i,'Feature_type']=='Centromere': 
   
        ax.vlines(x=i,ymin=0,ymax=0.8,linestyles='--',alpha=0.3)
        ax.text(x=i,y=0.6,s='centromere',rotation=90,fontsize=8)
#%% Transposon density for comparing Greg with Agnes dataset
fig=plt.figure(figsize=(10,9))
grid = plt.GridSpec(2, 1, wspace=0.0, hspace=0.0)
ax = plt.subplot(grid[0,0])
ax2 = plt.subplot(grid[1,0])   

ax.plot(data_wt_greg2['tr-density'],alpha=0.5,color='b')
ax.set_ylabel('transposond density: tn/bp')

## annotated centromeres
for i in np.arange(0,len(data_wt)):
    
    if data_wt.loc[i,'Feature_type']=='Centromere': 
   
        ax.vlines(x=i,ymin=0,ymax=0.8,linestyles='--',alpha=0.3)
        ax.text(x=i,y=0.6,s='centromere',rotation=90,fontsize=8)

ax2.plot(data_wt_agnes['tr-density'],alpha=0.5,color='orange')
ax2.set_ylabel('transposond density: tn/bp Agnes ')
ax2.set_xlabel('genes')
## annotated centromeres
for i in np.arange(0,len(data_wt_agnes)):
    
    if data_wt_agnes.loc[i,'Feature_type']=='Centromere': 
   
        ax2.vlines(x=i,ymin=0,ymax=0.8,linestyles='--',alpha=0.3)
        ax2.text(x=i,y=0.6,s='centromere',rotation=90,fontsize=8)

#%% saving the figure transposon density
fig.savefig('Transposon-density-WT-annotated-centromeres-Greg2_Agnes.png',dpi=300,format='png',transparent=False)
#%% Histograms of number of transposons per gene
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)
ax.set_xlabel('Number of transposons per gene')
ax.set_ylabel('number of genes (CDS)')
ax.hist(data_wt['Ninsertions'],bins=300,color='gray',alpha=0.7)
ax.set_xlim(0,200)
ax.vlines(x=data_wt_greg2['Ninsertions'].median(),ymin=0,ymax=1400,linestyles='--')
ax.text(x=data_wt_greg2['Ninsertions'].median(),y=1200,s='median all genes')
#%% saving the figure
fig.savefig('Number-of-transposons-per-gene-Greg2.png',dpi=300,format='png',transparent=False)


#%% #%% Histograms of number of transposon per gene according their essentiality
### create plot
fig=plt.figure(figsize=(10,9))
grid = plt.GridSpec(2, 1, wspace=0.0, hspace=0.0)

essentials=data_wt_greg2[data_wt_greg2.loc[:,'Essentiality']==1]['Ninsertions']
nonessentials=data_wt_greg2[data_wt_greg2.loc[:,'Essentiality']==0]['Ninsertions']


ax = plt.subplot(grid[0,0])
sns.histplot(nonessentials,binwidth=2,color='gray',alpha=0.7)
max_x = ax.get_xlim()
ax.set_xlim(left=0,right=200)
ax.grid(True)
ax.set_xticklabels([])
ax.vlines(x=nonessentials.median(),ymin=0,ymax=300,linestyles='--')
ax.text(x=nonessentials.median(),y=300,s='median-non-essentials')
ax.set_ylabel('Annotated non essential genes')

ax2 = plt.subplot(grid[1,0])    
sns.histplot(essentials,binwidth=2,color='orange',alpha=0.7)
ax2.invert_yaxis()
ax2.set_xlim(0,200)
ax2.grid(True)
ax2.vlines(x=essentials.median(),ymin=0,ymax=40,linestyles='--')
ax2.text(x=essentials.median(),y=40,s='median-essentials')

ax2.set_xlabel('Number of transposons per gene')
ax2.set_ylabel('Annotated essential genes')

plt.show()
#%% saving the figure
fig.savefig('Number-of-transposons-per-gene-according-essentiality-Greg2.png',dpi=300,format='png',transparent=False)
#%% Looking for regions devoid of transposons 

difficult_genes=data_wt_agnes[data_wt_agnes['Ninsertions']==0]

difficult_genes_essentials=difficult_genes[difficult_genes['Essentiality']==1]

#%% Loooking at correlations with essential genes 
fig = plt.figure(figsize=(11,5))
ax = fig.add_subplot(111)
sns.pairplot(data=data_wt_agnes,vars=['Ninsertions','tr-density'],hue='Essentiality')


#%% Looking at differences between libraries 
fig = plt.figure(figsize=(11,5))

plt.scatter(x=data_wt['tr-density'],y=data_wt_agnes['tr-density'])
plt.plot(np.linspace(0,1),np.linspace(0,1))