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

names_libraries={'wt':'WT-merged-all.xlsx','dnrp1':'dnrp1-merged-all.xlsx'}
data_library=[]
for i in names_libraries.keys():
 data_library.append(pd.read_excel('datasets/'+names_libraries[i],index_col='Unnamed: 0'))
#%% Creating a big dataframe of the libraries

data_library_pd=pd.concat(data_library,keys=names_libraries.keys())
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
#%%
freq=frequency_transposons(data_library_pd,names_libraries)
reads_per_transposon=reads_per_transposon(data_library_pd,names_libraries)
tr_density=transposon_density(data_library_pd,names_libraries)
#%%
analysis_libraries=defaultdict(dict)

j=0
for i in names_libraries.keys():
    
    analysis_libraries[i]['one-transposon-per-bp']=freq[j]
    analysis_libraries[i]['median-reads-per-transposons']=reads_per_transposon[j]
    analysis_libraries[i]['tr-density']=tr_density[j]
    j=j+1
    
analysis_libraries_pd=pd.DataFrame(analysis_libraries)

#%%
data_wt=data_library_pd.loc['wt'].copy()

#data_wt.index=data_wt['chromosome']
tn_median=data_wt['Ninsertions_truncatedgene'].median()
tn_essentials=data_wt[data_wt['Essentiality']==1]['Ninsertions_truncatedgene'].median()
tn_non_essentials=data_wt[data_wt['Essentiality']==0]['Ninsertions_truncatedgene'].median()

data_dnrp1=data_library_pd.loc['dnrp1'].copy()
#data_dnrp1.index=data_dnrp1['chromosome']

#%% Transposon density vs genes highlighting the centromere position


data_wt['tr-density']=data_wt['Ninsertions']/data_wt['Nbasepairs']
data_dnrp1['tr-density']=data_dnrp1['Ninsertions']/data_dnrp1['Nbasepairs']
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
#%% saving the figure transposon density
fig.savefig('Transposon-density-WT-without-annotated-centromeres.png',dpi=300,format='png',transparent=False)
#%% Histograms of number of transposons per gene
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)
ax.set_xlabel('Number of transposons per gene')
ax.set_ylabel('number of genes (CDS)')
ax.hist(data_wt['Ninsertions'],bins=300,color='gray',alpha=0.7)
ax.set_xlim(0,200)
ax.vlines(x=data_wt['Ninsertions'].median(),ymin=0,ymax=1400,linestyles='--')
ax.text(x=data_wt['Ninsertions'].median(),y=1200,s='median all genes')
#%% saving the figure
fig.savefig('Number-of-transposons-per-gene.png',dpi=300,format='png',transparent=False)

#%% Histograms of number of transposon per gene according their essentiality

fig = plt.figure(figsize=(11,5))
ax = fig.add_subplot(121)
ax.set_xlabel('Number of transposons per essential gene')
ax.set_ylabel('number of genes (CDS)')
essentials=data_wt[data_wt.loc[:,'Essentiality']==1]['Ninsertions']
ax.hist(essentials,bins=100,color='orange',alpha=0.7)
ax.set_xlim(0,200)
ax.vlines(x=essentials.median(),ymin=0,ymax=80,linestyles='--')
ax.text(x=essentials.median(),y=50,s='median-essentials')

ax2= fig.add_subplot(122)

ax2.set_xlabel('Number of transposons per non essential gene')
ax2.set_ylabel('number of genes (CDS)')
nonessentials=data_wt[data_wt.loc[:,'Essentiality']==0]['Ninsertions']
ax2.hist(nonessentials,bins=150,color='gray',alpha=0.7)
ax2.set_xlim(0,200)
ax2.vlines(x=nonessentials.median(),ymin=0,ymax=3000,linestyles='--')
ax2.text(x=nonessentials.median(),y=3000,s='median-non-essentials')

#%% saving the figure
fig.savefig('Number-of-transposons-per-gene-according-essentiality.png',dpi=300,format='png',transparent=False)
#%% Looking for regions devoid of transposons 

difficult_genes=data_wt[data_wt['Ninsertions']==0]

difficult_genes_essentials=difficult_genes[difficult_genes['Essentiality']==1]

#%% Loooking at correlations with essential genes 
fig = plt.figure(figsize=(11,5))
ax = fig.add_subplot(111)
sns.pairplot(data=data_wt,vars=['Ninsertions','tr-density'],hue='Essentiality')


#%% Looking at differences between libraries 


plt.scatter(x=data_wt['tr-density'],y=data_dnrp1['tr-density'])
plt.plot(np.linspace(0,1),np.linspace(0,1))