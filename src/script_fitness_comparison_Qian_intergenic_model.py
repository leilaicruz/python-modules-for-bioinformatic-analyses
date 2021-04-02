# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 13:21:50 2021

@author: linigodelacruz
"""
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import os
import seaborn as sns
import scipy 

#%% Importing datasets

conversion=pd.read_csv('datasets/from_systematic_genenames2standard_genenames.csv',delimiter=',')
conversion.columns=['systematic_name','standard_name']

qian_dataset=pd.read_excel('datasets/Gowthrate_Qian.xls',header=1)
qian_dataset_SC=qian_dataset.loc[:,['ORF','SC fitness']]

fitness_intergenic=pd.read_excel('datasets/fitness-from-intergenic-model-agnes-sequencing.xlsx')
fitness_intergenic.index=fitness_intergenic['gene_name']
#%% Replacing the systematic names from QIan to standard names
j=0
for i in qian_dataset_SC['ORF']:
    value=conversion[conversion.loc[:,'systematic_name']==i]['standard_name']
    if len(value)!=0 :
        qian_dataset_SC.loc[j,'ORF_standard']=value.tolist()[0]
    j=j+1
#%%
qian_dataset_SC.dropna(inplace=True)
qian_dataset_SC.index=qian_dataset_SC['ORF_standard']

#%% Comparing relative fitness by numbers
column_wt='wt_rates_intergenic_merged_hand'
column_dnrp1='dnrp1_rates_intergenic_merged_hand'
norm_wt=fitness_intergenic.loc['HO',column_wt]

norm_wt_fitness=fitness_intergenic[column_wt]/norm_wt


qian_dataset_SC_fitness=qian_dataset_SC['SC fitness']

similar_fitness=[]


for genes in qian_dataset_SC.index:
    if genes  in norm_wt_fitness.index: 
        if np.isclose(qian_dataset_SC_fitness[genes],norm_wt_fitness[genes],0.05) :
            similar_fitness.append(genes)

print('The percentage of genes that have similar fitnesses are:', 100*len(similar_fitness)/len(qian_dataset_SC_fitness))
#%% Comparing with the avergae of all datasets 
average_all_datasets_intergenic=np.mean(fitness_intergenic.loc[:,['average_wt','average_dnrp1']],axis=1)
average_all_datasets_intergenic_normed=average_all_datasets_intergenic/average_all_datasets_intergenic['HO']
average_merged=np.mean(fitness_intergenic.loc[:,[column_dnrp1,column_wt]],axis=1)
average_merged_normed=average_merged/average_merged['HO']
#%% Comparing relative fitness scatter

fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)

for genes in qian_dataset_SC.index:
    if genes in norm_wt_fitness.index:
        #ax.scatter(x=norm_wt_fitness[genes],y=qian_dataset_SC_fitness[genes],color='blue',alpha=0.6)
        ax.scatter(x=average_all_datasets_intergenic_normed[genes],y=qian_dataset_SC_fitness[genes],color='blue',alpha=0.6)

ax.set_ylim(0,1.5)
ax.set_xlim(0,1.5)
ax.set_ylabel('Qian  SC fitness values ')
ax.set_xlabel('Intergenic model fitness relative values to HO ')

ax.plot(np.linspace(0,1.5),np.linspace(0,1.5),color='black',linewidth=5,alpha=0.6)
#%% histograms
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(121)
#a=fitness_intergenic[fitness_intergenic[column_dnrp1]!=-np.inf][column_dnrp1]
a=average_all_datasets_intergenic[average_all_datasets_intergenic!=-np.inf]
a_new=a[a!=np.inf]
a_new_norm=a_new/norm_wt
ax.hist(a_new_norm)
ax.set_xlabel('Satay to HO')

ax1 = fig.add_subplot(122)
a=qian_dataset_SC['SC fitness']
# a_new=a[a!=np.inf]
# a_new_norm=a_new/norm_wt
ax1.hist(a)
ax1.set_xlabel('Qian')

#%% Fitness differences among replicates
std_ho=fitness_intergenic.loc['HO','std_wt']
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)


x=fitness_intergenic
y=fitness_intergenic
#ax.scatter(x=x['nrp1_1_a']/norm_wt,y=y['nrp1_1_b']/norm_wt,color='blue',alpha=0.6)
#ax.scatter(x=x['average_wt']/norm_wt,y=y['std_wt']/std_ho,color='blue',alpha=0.6)
#ax.set_xlabel('nrp1_1_a-nrp1_1_b')
# ax.set_ylabel('nrp1_1_b')
#ax.plot(np.linspace(0,2.5),np.linspace(0,2.5),color='black',linewidth=5,alpha=0.6)

#ax.hist(x['nrp1_1_a']/norm_wt-y['nrp1_1_b']/norm_wt,bins=50,label='wt_a',color='blue',alpha=0.6)
#ax.hist(average_merged_normed,color='orange',alpha=0.7,label='wt_b')
# ax.legend()
#%%
replicate_differences=np.abs(x['wt_a']/norm_wt-x['wt_b']/norm_wt)

#%% save this figure

fig.savefig('average hist_all_datasets_fitness_Qian_vs_intergenic.png',format='png',dpi=300,transparent=False)
