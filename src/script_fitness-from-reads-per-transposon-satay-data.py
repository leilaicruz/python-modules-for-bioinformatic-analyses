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

def getting_r(datasets): 
    total_reads=[]
    T=90
    r=[]
    reads_per_transposons=[]
    for i in np.arange(0,len(datasets)):
        total_reads=np.sum(datasets[i]['number_of_read_per_gene'])
        total_tn=np.sum(datasets[i]['number_of_transposon_per_gene'])
        K=total_reads/total_tn
        N=datasets[i]['number_of_read_per_gene']/datasets[i]['number_of_transposon_per_gene']
        r.append(np.log(N/(1-N/K))/T)
        reads_per_transposons.append(N)
        
       
    
    return r,reads_per_transposons

#%% Getting the rates from the function

rates,reads_per_transposons=getting_r(datasets)
#%% Reads per transposon dataframe 
readpertransposon_pd=pd.DataFrame(reads_per_transposons,index=np.arange(0,len(files)))
readpertransposon_pd=readpertransposon_pd.T
readpertransposon_pd.columns=['nrp1_1_a','nrp1_1_b','nrp1_2_a','nrp2_2_b','wt_a','wt_b']
readpertransposon_pd.index=datasets[0]['gene_name']
readpertransposon_pd.fillna(0,inplace=True)
readpertransposon_pd.replace(float('-inf'),0,inplace=True)

readpertransposon_pd['average_dnrp1']=np.mean(readpertransposon_pd.loc[:,['nrp1_1_a','nrp1_1_b','nrp1_2_a','nrp2_2_b']],axis=1)
readpertransposon_pd['std_dnrp1']=np.std(readpertransposon_pd.loc[:,['nrp1_1_a','nrp1_1_b','nrp1_2_a','nrp2_2_b']],axis=1)

readpertransposon_pd['average_wt']=np.mean(readpertransposon_pd.loc[:,['wt_a','wt_b']],axis=1)
readpertransposon_pd['std_wt']=np.std(readpertransposon_pd.loc[:,['wt_a','wt_b']],axis=1)

readpertransposon_pd['average_differences_wt-mutant']=readpertransposon_pd['average_wt']-readpertransposon_pd['average_dnrp1']
#%% rates dataframe

rates_pd=pd.DataFrame(rates,index=np.arange(0,len(files)))
rates_pd=rates_pd.T
rates_pd.columns=['nrp1_1_a','nrp1_1_b','nrp1_2_a','nrp2_2_b','wt_a','wt_b']
rates_pd.index=datasets[0]['gene_name']
rates_pd.fillna(0,inplace=True)
rates_pd.replace(float('-inf'),0,inplace=True)
rates_normalized2HO=np.divide(rates_pd,rates_pd.loc['HO',:])


rates_pd['average_dnrp1']=np.mean(rates_pd.loc[:,['nrp1_1_a','nrp1_1_b','nrp1_2_a','nrp2_2_b']],axis=1)
rates_pd['std_dnrp1']=np.std(rates_pd.loc[:,['nrp1_1_a','nrp1_1_b','nrp1_2_a','nrp2_2_b']],axis=1)

rates_pd['average_wt']=np.mean(rates_pd.loc[:,['wt_a','wt_b']],axis=1)
rates_pd['std_wt']=np.std(rates_pd.loc[:,['wt_a','wt_b']],axis=1)

rates_pd['average_differences_wt-mutant']=rates_pd['average_wt']-rates_pd['average_dnrp1']

rates_pd['norm2max_wt']=rates_pd.loc[:,'average_wt']/rates_pd.loc[:,'average_wt'].max()
rates_pd['norm2max_dnrp1']=rates_pd.loc[:,'average_dnrp1']/rates_pd.loc[:,'average_dnrp1'].max()

#%% Plotting the differences from dnrp1 fitness vs WT fitness levels 

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)

ax.hist(rates_pd['average_differences_wt-mutant'])

ax.set_xlabel('WT_fitness-dnrp1_fitness')
#%% significant genes
test_rates=[]
test_readpertransposon=[]

for i in np.arange(0,len(rates_pd)):
    test_rates.append(scipy.stats.ttest_ind(rates_pd.loc[rates_pd.index[i],['nrp1_1_a','nrp1_1_b','nrp1_2_a','nrp2_2_b']], rates_pd.loc[rates_pd.index[i],['wt_a','wt_b']]))
    test_readpertransposon.append(scipy.stats.ttest_ind(readpertransposon_pd.loc[readpertransposon_pd.index[i],['nrp1_1_a','nrp1_1_b','nrp1_2_a','nrp2_2_b']], readpertransposon_pd.loc[readpertransposon_pd.index[i],['wt_a','wt_b']]))
#%% FC and p-value or volcano on fitness 

rates_pd['FC']=np.divide(rates_pd['average_wt'],rates_pd['average_dnrp1'])

readpertransposon_pd['FC']=np.divide(readpertransposon_pd['average_wt'],readpertransposon_pd['average_dnrp1'])
## excluding zeros from FC
for i in np.arange(0,len(rates_pd)):
    if rates_pd['FC'][i]!=0:
        value=rates_pd['FC'][i]
        rates_pd.loc[rates_pd.index[i],'fc_log2']=np.log2(value)
        rates_pd.loc[rates_pd.index[i],'-log_p']=-1*np.log10(test_rates[i][1])
    if readpertransposon_pd['FC'][i]!=0:
        value=readpertransposon_pd['FC'][i]
        readpertransposon_pd.loc[readpertransposon_pd.index[i],'fc_log2']=np.log2(value)
        readpertransposon_pd.loc[readpertransposon_pd.index[i],'-log_p']=-1*np.log10(test_readpertransposon[i][1])
        
        

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

#%% function to filter the data 
def filtering_significant_genes(dataset,p_value_th,fc_th):
    dataset_highfc=dataset[dataset.loc[:,'fc_log2']>fc_th] 
    dataset_highp=dataset[dataset.loc[:,'-log_p']>p_value_th]
    dataset_significant_genes_wt= pd.merge(dataset_highfc, dataset_highp, how='inner', on=['-log_p'],left_index=True, right_index=True)
    dataset_significant_genes_wt.drop(columns=['fc_log2_y'],inplace=True)
    dataset_significant_genes_wt.rename(columns={'fc_log2_x': 'fc_log2'},inplace=True)
    
    
    dataset_highfc=dataset[dataset.loc[:,'fc_log2']<-fc_th] 
    dataset_highp=dataset[dataset.loc[:,'-log_p']>p_value_th] 
    dataset_significant_genes_mutant= pd.merge(dataset_highfc, dataset_highp, how='inner', on=['-log_p'],left_index=True, right_index=True)
    dataset_significant_genes_mutant.drop(columns=['fc_log2_y'],inplace=True)
    dataset_significant_genes_mutant.rename(columns={'fc_log2_x': 'fc_log2'},inplace=True)
    
    return dataset_significant_genes_wt,dataset_significant_genes_mutant


def viz_data(dataset,p_value_th,fc_th,filenames,saving=True,plotting=True):
    
    results=[]
    if type(dataset)==list:
        
        for i in np.arange(0,len(dataset)):
            
            ref,exp=filtering_significant_genes(dataset[i],p_value_th,fc_th)
            results.append([ref,exp])
    else : 
        
        ref,exp=filtering_significant_genes(dataset,p_value_th,fc_th)
        results.append([ref,exp])        
### Plotting all datasets 
    for data in np.arange(0,len(results)):
        fig = plt.figure(figsize=(10,5))
        ax = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        
        data_ref=results[data][0]
        data_exp=results[data][1]
        
        for i in range(len(data_ref)):
            x = data_ref.loc[data_ref.index[i],'fc_log2']
            y = data_ref.loc[data_ref.index[i],'-log_p']
            ax.plot(x, y, 'bo')
        
            ax.text(x*(1-0.02) , y*(1+0.02)  , data_ref.index[i], fontsize=4)
         
        for i in range(len(data_exp)):
            x = data_exp.loc[data_exp.index[i],'fc_log2']
            y = data_exp.loc[data_exp.index[i],'-log_p']
            ax2.plot(x, y, 'bo')
           
            ax2.text(x*(1-0.02) , y*(1+0.02)  , data_exp.index[i], fontsize=4)
        
        for axes in [ax,ax2]:
            axes.set_ylim((0, 10))
            axes.set_xlabel('fc_log2')
            axes.set_ylabel('-log_p')
            if axes==ax:
                axes.set_xlim((1, 5*fc_th))
                axes.set_title('Signif for WT-'+ filenames,fontsize=8)
            else:
                axes.set_xlim((-5*fc_th, -1))
                axes.set_title('Signif for mutant-'+ filenames,fontsize=8)
                
        if plotting==True:
            plt.show()
        if saving==True:
            fig.savefig(filenames + '.png',dpi=300,format='png',transparent=False)

#%%
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

plt.plot(rates_pd.loc[:,'norm2max_wt'],rates_pd.loc[:,'norm2max_dnrp1'],'bo',alpha=0.3)

plt.plot(np.linspace(0,1),np.linspace(0,1)*rates_pd.loc['NRP1','norm2max_wt'],linewidth=4,alpha=0.4,color='k')
plt.plot(np.linspace(0,1),np.linspace(0,1),linewidth=4,alpha=0.4,color='k')

plt.xlim(0,1)
plt.ylim(0,1)
plt.xlabel('single knockout fitness')
plt.ylabel('double knockout fitness: dnrp1dgenex')
plt.savefig('dnrp1-fitness-map-from-satay.png',dpi=300,transparent=True,format='png')

#%% fitness annotating the significant genes
log_p_value_th=2
log2_fc_th=7
ref,exp=filtering_significant_genes(rates_pd,log_p_value_th,log2_fc_th)
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)
for i in range(len(rates_pd)):
    x = rates_pd.loc[rates_pd.index[i],'norm2max_wt']
    y = rates_pd.loc[rates_pd.index[i],'norm2max_dnrp1']
    ax.plot(x, y, 'bo')
    if rates_pd.index[i] in ref.index :
        ax.text(x*(1-0.02) , y*(1+0.02)  ,rates_pd.index[i] , fontsize=4,rotation=45)
    elif rates_pd.index[i] in exp.index:
        ax.text(x*(1-0.02) , y*(1+0.02)  ,rates_pd.index[i] , fontsize=4)
        
#%% saving fitness plots
ax.set_xlim(0,1)
ax.set_ylim(-0.25,1)
ax.set_xlabel('single knockout fitness')
ax.set_ylabel('double knockout fitness: dnrp1dgenex')
fig.savefig('significant-dnrp1-fitness-map-from-satay.png',dpi=300,transparent=False,format='png')

#%% Intergenic model per gene in the population

K=np.sum(readpertransposon_pd.loc[:,'average_wt'])
gene='MRPL44'

t=np.linspace(0,100)
rm_wt=np.abs(rates_pd.loc[gene,'norm2max_wt'])
N_wt=np.exp(rm_wt*t)/(1+np.exp(rm_wt*t)/K)


rm_dnrp1=np.abs(rates_pd.loc[gene,'norm2max_dnrp1'])
N_dnrp1=np.exp(rm_dnrp1*t)/(1+np.exp(rm_dnrp1*t)/K)

fig = plt.figure(figsize=(12,5))
ax = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax.scatter(x=t,y=N_wt,label='WT_backg')
ax2.scatter(x=t,y=N_dnrp1,label='dnrp1_back',color='k')

for axes in [ax,ax2]:
    axes.legend()
    axes.set_ylabel('Reads per transposons')
    axes.set_xlabel('Time-hours')

