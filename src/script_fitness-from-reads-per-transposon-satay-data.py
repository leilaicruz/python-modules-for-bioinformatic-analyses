# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 16:27:36 2021

@author: linigodelacruz
"""

import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import os
import seaborn as sns
import scipy 

from src.python_modules.module_intergenic_model import getting_r,pd_to_compare_libraries,filtering_significant_genes,viz_data
#%% # Reads and transposon dataset
path='datasets/Leila_10Feb21_PerGeneFiles/'
files=os.listdir(path) # name of all the files 
datasets=[]
for i in np.arange(0,len(files)):
    datasets.append(pd.read_csv(path+files[i],delimiter=" "))

## Removing the ADE2 and URA3 gene reads and insertions from the population
for i in np.arange(0,len(datasets)):
    datasets[i]=datasets[i][datasets[i].gene_name != 'ADE2']
    datasets[i]=datasets[i][datasets[i].gene_name != 'URA3']

#%%  From reads to fitness from a intergenic competition model
## asumming we start with zero reads so one copy per cell
# r=K/T(K-N)
# K=total number of reads/total of transposons (carrying capacity) per library- excluding the ade2 and URA3 gene
# T=90
# N=number of reads per gene/transposons (excluding the ade2 gene)

## Getting the rates from the function : getting_r

rates,reads_per_transposons=getting_r(datasets,time=False)
#%% Preparing datasets for analyses

rates_pd=pd_to_compare_libraries(rates,datasets)
readpertransposon_pd=pd_to_compare_libraries(reads_per_transposons,datasets)

#%% Plotting the differences from dnrp1 fitness vs WT fitness levels 

fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)
#ax.scatter(readpertransposon_pd['average_wt'],readpertransposon_pd['average_dnrp1'],alpha=0.5)
ax.set_xlabel('average reads per transposon WT')
ax.set_ylabel('average reads per transposon dnrp1')

# ax2 = fig.add_subplot(122)
# ax2.hist(readpertransposon_pd['average_dnrp1'])
# ax2.set_xlabel('average reads per transposon dnrp1')
# ax2.set_ylabel('counts')
for i in np.arange(0,len(rates_pd)):
    x=readpertransposon_pd['average_wt'][i]
    y=readpertransposon_pd['average_dnrp1'][i]
    ax.scatter(x,y,alpha=0.5,color='b')
    if np.abs(x-y)>35:
        ax.text(x*(1-0.02) ,y*(1+0.02) , readpertransposon_pd.index[i] , fontsize=3)

ax.set_xlim(0,150)
ax.set_ylim(0,150)
#%% saving the plot
fig.savefig('output_images/dnrp1-vs-wt-comparison-reads-per-transposons.png',format='png',dpi=300,transparent=True)
#%% significant genes
test_rates=[]
test_readpertransposon=[]

for i in np.arange(0,len(rates_pd)):
    test_rates.append(scipy.stats.ttest_ind(rates_pd.loc[rates_pd.index[i],['nrp1_1_a','nrp1_1_b','nrp1_2_a','nrp2_2_b']], rates_pd.loc[rates_pd.index[i],['wt_a','wt_b']]))
    test_readpertransposon.append(scipy.stats.ttest_ind(readpertransposon_pd.loc[readpertransposon_pd.index[i],['nrp1_1_a','nrp1_1_b','nrp1_2_a','nrp2_2_b']], readpertransposon_pd.loc[readpertransposon_pd.index[i],['wt_a','wt_b']]))
### FC and p-value or volcano on fitness 

for i in np.arange(0,len(rates_pd)):
    if rates_pd['average_wt'][i]!=0 and rates_pd['average_dnrp1'][i]!=0:
        rates_pd.loc[rates_pd.index[i],'FC']=np.divide(rates_pd.loc[rates_pd.index[i],'average_wt'],rates_pd.loc[rates_pd.index[i],'average_dnrp1'])
        value=rates_pd['FC'][i]
        rates_pd.loc[rates_pd.index[i],'fc_log2']=np.log2(value)
        rates_pd.loc[rates_pd.index[i],'-log_p']=-1*np.log10(test_rates[i][1])
    else:
        rates_pd.loc[rates_pd.index[i],'FC']=0
        rates_pd.loc[rates_pd.index[i],'fc_log2']=0
        rates_pd.loc[rates_pd.index[i],'-log_p']=0
    if readpertransposon_pd['average_wt'][i]!=0 and readpertransposon_pd['average_dnrp1'][i]!=0:
        readpertransposon_pd.loc[readpertransposon_pd.index[i],'FC']=np.divide(readpertransposon_pd.loc[rates_pd.index[i],'average_wt'],readpertransposon_pd.loc[readpertransposon_pd.index[i],'average_dnrp1'])
        value=readpertransposon_pd['FC'][i]
        readpertransposon_pd.loc[readpertransposon_pd.index[i],'fc_log2']=np.log2(value)
        readpertransposon_pd.loc[readpertransposon_pd.index[i],'-log_p']=-1*np.log10(test_readpertransposon[i][1])
    else:
        readpertransposon_pd.loc[readpertransposon_pd.index[i],'FC']=0
        readpertransposon_pd.loc[readpertransposon_pd.index[i],'fc_log2']=0
        readpertransposon_pd.loc[readpertransposon_pd.index[i],'-log_p']=0
        
        

rates_pd['fc_log2'].fillna(0,inplace=True)
rates_pd['-log_p'].fillna(0,inplace=True)

readpertransposon_pd['fc_log2'].fillna(0,inplace=True)
readpertransposon_pd['-log_p'].fillna(0,inplace=True)
#%% volcano

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(121)
ax.plot(rates_pd['fc_log2'],rates_pd['-log_p'],'bo')
ax.set_xlim(-10,10)
ax.set_title('fitness-from-intergenic-model')

ax2 = fig.add_subplot(122)
ax2.plot(readpertransposon_pd['fc_log2'],readpertransposon_pd['-log_p'],'bo')
ax2.set_xlim(-10,10)
ax2.set_title('reads-per-transposons')

for i in np.arange(0,len(rates_pd)):
    x=rates_pd['fc_log2'][i]
    y=rates_pd['-log_p'][i]
    if np.abs(x)>5 and y>2:
        ax.plot(x, y, 'ro')
    ax.text(x*(1-0.02) ,y*(1+0.02) , rates_pd.index[i] , fontsize=5)
    
for i in np.arange(0,len(readpertransposon_pd)):
    x=readpertransposon_pd['fc_log2'][i]
    y=readpertransposon_pd['-log_p'][i]
    if np.abs(x)>4 and y>1: 
        ax2.plot(x, y, 'ro')
    ax2.text(x*(1-0.02) ,y*(1+0.02) , readpertransposon_pd.index[i] , fontsize=5)

#%% saving figure

fig.savefig('volcano-fitness-alldnrp1v-s-WT.png', format='png',dpi=300)
#%% exporting datasets

readpertransposon_pd.to_excel('datasets/reads-per-transposons-agnes-sequencing.xlsx')
rates_pd.to_excel('datasets/fitness-from-intergenic-model-agnes-sequencing.xlsx')

#%% Filtering significant genes

#%% Vixualizating the most significant genes for all the data inside Volcano folders
log_p_value_th=2
log2_fc_th=4
viz_data(readpertransposon_pd,log_p_value_th,log2_fc_th,filenames='reads-per-transposon')

#%%
viz_data(rates_pd,log_p_value_th,log2_fc_th,filenames='fitness-from-intergenic-model')
# log_p_value_th_reads=1
# log2_fc_th_reads=3.3
# viz_data(datasets_reads,log_p_value_th_reads,log2_fc_th_reads,files_reads)

#%% Fitness plots

plt.plot(rates_pd.loc[:,'norm2max_wt'],rates_pd.loc[:,'norm2max_dnrp1'],'bo',alpha=0.2)

plt.plot(np.linspace(0,1),np.linspace(0,1)*rates_pd.loc['NRP1','norm2max_wt'],linewidth=4,alpha=0.4,color='k')
plt.plot(np.linspace(0,rates_pd.loc['NRP1','norm2max_wt']),np.linspace(0,rates_pd.loc['NRP1','norm2max_wt']),linewidth=4,alpha=0.4,color='k')
plt.hlines(rates_pd.loc['NRP1','norm2max_wt'],rates_pd.loc['NRP1','norm2max_wt'],1,linewidth=4,alpha=0.4,color='k')

plt.xlim(0,1)
plt.ylim(0,1)
plt.xlabel('single knockout fitness')
plt.ylabel('double knockout fitness: dnrp1dgenex')
#plt.savefig('dnrp1-fitness-map-from-satay.png',dpi=300,transparent=True,format='png')

#%% fitness annotating the significant genes
log_p_value_th=2
log2_fc_th=3
ref,exp=filtering_significant_genes(rates_pd,log_p_value_th,log2_fc_th)
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)
for i in range(len(rates_pd)):
    x = rates_pd.loc[rates_pd.index[i],'norm2max_wt']
    y = rates_pd.loc[rates_pd.index[i],'norm2max_dnrp1']
    ax.plot(x, y, 'bo',alpha=0.2)
    if rates_pd.index[i] in ref.index :
        ax.text(x*(1-0.02) , y*(1+0.02)  ,rates_pd.index[i] , fontsize=6,rotation=45)
    elif rates_pd.index[i] in exp.index:
        ax.text(x*(1-0.02) , y*(1+0.02)  ,rates_pd.index[i] , fontsize=6)
        
#%% saving fitness plots
ax.set_xlim(0,1)
ax.set_ylim(0,1)
ax.set_xlabel('single knockout fitness')
ax.set_ylabel('double knockout fitness: dnrp1dgenex')
fig.savefig('significant-dnrp1-fitness-map-from-satay.png',dpi=300,transparent=False,format='png')

#%% Intergenic model per gene in the population

K_wt=np.sum(readpertransposon_pd.loc[:,'average_wt'])# carrying capacity
gene='BEM3'
#gene=MED11 essential for WT and not for nrp1
# gene = MPM1  # max for wt
t=np.linspace(0,90)
rm_wt=np.abs(rates_pd.loc[gene,'average_wt'])
N_wt=np.exp(rm_wt*t)/(1+np.exp(rm_wt*t)/K_wt)

K_dnrp1=np.sum(readpertransposon_pd.loc[:,'average_dnrp1'])# carrying capacity
rm_dnrp1=np.abs(rates_pd.loc[gene,'average_dnrp1'])
N_dnrp1=np.exp(rm_dnrp1*t)/(1+np.exp(rm_dnrp1*t)/K_dnrp1)
# gene = PFA5 # max for dnrp1 
fig = plt.figure(figsize=(13,5))
ax = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax.scatter(x=t,y=N_wt,label='WT_backg',color='k')
ax.set_ylim(0,K_wt/100+100)
ax.hlines(K_wt/100,0,t[-1],label='K/100')
ax2.scatter(x=t,y=N_dnrp1,label='dnrp1_back',color='blue')
ax2.set_ylim(0,K_dnrp1/100+100)
ax2.hlines(K_dnrp1/100,0,t[-1],label='K/100')


for axes in [ax,ax2]:
    axes.legend()
    axes.set_ylabel('Reads per transposons')
    axes.set_xlabel('Time-hours')
    axes.set_title(gene)

#%% saving the figure

fig.savefig('intergenic_model_growth_gene_'+gene+'.png',dpi=300,format='png')