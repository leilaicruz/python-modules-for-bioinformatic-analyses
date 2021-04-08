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

from src.python_modules.module_analysis_transposon_sites import *
#%% Import of dataframes output from the SATAY pipeline

names_libraries={'wt':'WT-merged-all.xlsx','dnrp1':'dnrp1-merged-all.xlsx','wt_agnes':'WT-merged-Agnes.xlsx','wt_strong_alig':'WT-merged-astringent-alignment.xlsx'}
data_library=[]
for i in names_libraries.keys():
 data_library.append(pd.read_excel('datasets/'+names_libraries[i],index_col='Unnamed: 0'))
#%% Creating a big dataframe of the libraries

data_library_pd=pd.concat(data_library,keys=names_libraries.keys(),sort=True)
data_library_pd.fillna(0,inplace=True)


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


#%% Transposon density vs genes 


data_wt['tr-density']=data_wt['Ninsertions']/data_wt['Nbasepairs']
data_wt_agnes['tr-density']=data_wt_agnes['Ninsertions']/data_wt_agnes['Nbasepairs']
data_wt_greg2['tr-density']=data_wt_greg2['Ninsertions']/data_wt_greg2['Nbasepairs']
#%% Reads per transposons
# data_wt['reads-per-tr']=(data_wt['Nreads']/data_wt['Ninsertions'])/data_wt['Nbasepairs']
# data_wt_agnes['reads-per-tr']=(data_wt_agnes['Nreads']/data_wt_agnes['Ninsertions'])/data_wt_agnes['Nbasepairs']

data_wt['reads-per-tr']=(data_wt['Nreads']/data_wt['Ninsertions'])
data_wt_agnes['reads-per-tr']=(data_wt_agnes['Nreads']/data_wt_agnes['Ninsertions'])
#%% Plot transposon density (fig 1B Benoit) highlighting the centromere position

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

ax.plot(data_wt['tr-density'],alpha=0.5,color='b')
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
#%%  Plot reads per transposon  highlighting the centromere position

fig=plt.figure(figsize=(10,9))
grid = plt.GridSpec(2, 1, wspace=0.0, hspace=0.0)
ax = plt.subplot(grid[0,0])
ax2 = plt.subplot(grid[1,0])   

ax.plot(data_wt['reads-per-tr'],alpha=0.7,color='b')
ax.set_ylabel('reads per tr per bp')
ax.set_xlabel('genes')
## annotated centromeres
for i in np.arange(0,len(data_wt)):
    
    if data_wt.loc[i,'Feature_type']=='Centromere': 
   
        ax.vlines(x=i,ymin=0,ymax=80,linestyles='--',alpha=0.3)
        #ax.text(x=i,y=4000,s='centromere',rotation=90,fontsize=8)
    elif data_wt.loc[i,'reads-per-tr']>15:
        ax.vlines(x=i,ymin=0,ymax=80,linestyles='-',alpha=0.5)
        ax.text(x=i,y=70,s=data_wt.loc[i,'Standard_name'],rotation=90,fontsize=8)
        

ax2.plot(data_wt['tr-density'],alpha=0.7,color='b')
ax2.set_ylabel('transposond density: tn/bp')

## annotated centromeres
for i in np.arange(0,len(data_wt)):
    
    if data_wt.loc[i,'Feature_type']=='Centromere': 
   
        ax2.vlines(x=i,ymin=0,ymax=0.8,linestyles='--',alpha=0.3)
        ax2.text(x=i,y=0.6,s='centromere',rotation=90,fontsize=8)
#%% saving the figure reads per transposon density
fig.savefig('Reads-per-tr-merged-WT-Greg-along-genome.png',dpi=300,format='png',transparent=False)

#%% determine the local variation of transposons along te genome 
## Data per chromosome
mean_wt_chrom=data_wt.groupby(by='chromosome')['Ninsertions'].mean()
std_wt_chrom=data_wt.groupby(by='chromosome')['Ninsertions'].std()

mean_wt_chrom_trdensity=data_wt.groupby(by='chromosome')['tr-density'].mean()
std_wt_chrom_trdensity=data_wt.groupby(by='chromosome')['tr-density'].std()

mean_wt_chrom_readspertr=data_wt.groupby(by='chromosome')['reads-per-tr'].mean()
std_wt_chrom_readspertr=data_wt.groupby(by='chromosome')['reads-per-tr'].std()

fig=plt.figure(figsize=(10,9))
grid = plt.GridSpec(3, 1, wspace=0.0, hspace=0.0)
ax = plt.subplot(grid[0,0])
ax.errorbar(data_wt.loc[:,'chromosome'].unique(), mean_wt_chrom, std_wt_chrom, marker='s', mfc='red',
         mec='green', ms=10, mew=1,capsize=4)
ax.set_ylabel('Ninsertions')

ax2 = plt.subplot(grid[1,0])
ax2.errorbar(data_wt.loc[:,'chromosome'].unique(), mean_wt_chrom_trdensity, std_wt_chrom_trdensity, marker='s', mfc='red',
         mec='green', ms=10, mew=1,capsize=4)
ax2.set_ylabel('Ninsertions per bp')

ax3 = plt.subplot(grid[2,0])
ax3.errorbar(data_wt.loc[:,'chromosome'].unique(), mean_wt_chrom_readspertr,std_wt_chrom_readspertr, marker='s', mfc='red',
         mec='green', ms=10, mew=1,capsize=4)
ax3.set_ylabel('Nreads per Ninsertions')
#%% assesing local variation per chromosome. Viz per chromosome

magnitudes=['Ninsertions','tr-density','reads-per-tr']
chromosomes=data_wt.loc[:,'chromosome'].unique()

windows=10
chrom=chromosomes[4]


fig=plt.figure(figsize=(10,9))
grid = plt.GridSpec(3, 1, wspace=0.0, hspace=0.0)
ax = plt.subplot(grid[0,0])
ax.set_title('Errorbars over chrom='+str(chrom)+', every '+str(windows)+'genes')
mean_over_chromI,std_over_chromI=local_variation(chrom=chrom, windows=windows, data=data_wt,column=magnitudes[0])
chrom_data=pd.DataFrame([mean_over_chromI,std_over_chromI],index=['mean','std'])


ax.errorbar(np.arange(0,len(data_wt[data_wt.loc[:,'chromosome']==chrom]),windows), chrom_data.loc['mean',:], chrom_data.loc['std',:], marker='s', 
            mfc='red', mec='green', ms=10, mew=1,capsize=4)
ax.set_ylabel(magnitudes[0])

mean_over_chromI,std_over_chromI=local_variation(chrom=chrom, windows=windows, data=data_wt,column=magnitudes[1])
chrom_data=pd.DataFrame([mean_over_chromI,std_over_chromI],index=['mean','std'])

ax2= plt.subplot(grid[1,0])
ax2.errorbar(np.arange(0,len(data_wt[data_wt.loc[:,'chromosome']==chrom]),windows), chrom_data.loc['mean',:], chrom_data.loc['std',:], marker='s', mfc='red',
         mec='green', ms=10, mew=1,capsize=4)
ax2.set_ylabel(magnitudes[1])

mean_over_chromI,std_over_chromI=local_variation(chrom=chrom, windows=windows, data=data_wt,column=magnitudes[2])
chrom_data=pd.DataFrame([mean_over_chromI,std_over_chromI],index=['mean','std'])

ax3= plt.subplot(grid[2,0])
ax3.errorbar(np.arange(0,len(data_wt[data_wt.loc[:,'chromosome']==chrom]),windows), chrom_data.loc['mean',:], chrom_data.loc['std',:], marker='s', mfc='red',
         mec='green', ms=10, mew=1,capsize=4)
ax3.set_ylabel(magnitudes[2])
#ax.set_xlabel('genes along the windows')

#%% assesing local variation per chromosome. Visualizing more than one chromosome

magnitudes=['Ninsertions','tr-density','reads-per-tr']
chromosomes=data_wt.loc[:,'chromosome'].unique()

windows=10

chroms=[chromosomes[0],chromosomes[1],chromosomes[2],chromosomes[3]]

fig=plt.figure(figsize=(20,20))
grid = plt.GridSpec(3, len(chroms), wspace=0.0, hspace=0.0)

for i in np.arange(0,len(chroms)):
    for j in np.arange(0,len(magnitudes)): 
        ax = plt.subplot(grid[j,i])
       
        mean_over_chromI,std_over_chromI=local_variation(chrom=chroms[i], windows=windows, data=data_wt,column=magnitudes[j])
        chrom_data=pd.DataFrame([mean_over_chromI,std_over_chromI],index=['mean','std'])
        
        
        ax.errorbar(np.arange(0,len(data_wt[data_wt.loc[:,'chromosome']==chroms[i]]),windows), chrom_data.loc['mean',:], chrom_data.loc['std',:], marker='s', 
                    mfc='red', mec='green', ms=10, mew=1,capsize=4)
        
        
    ax.set_title('Errorbars over chrom='+str(chroms[i])+', every '+str(windows)+'genes')

fig.text(0.08, 0.5, magnitudes[1], va='center', rotation='vertical',fontsize=15)  
fig.text(0.08, 0.75, magnitudes[0], va='center', rotation='vertical',fontsize=15)   
fig.text(0.08, 0.25, magnitudes[2], va='center', rotation='vertical',fontsize=15)  
# axlabel.set_ylabel(magnitudes[0],fontsize=15)     
#ax.set_xlabel('genes along the windows')

#%% saving the figure
fig.savefig('local-variation-merged_wt-windows-'+str(windows)+'genes-'+ str(chroms)+'.png',dpi=300,format='png',transparent=False)
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