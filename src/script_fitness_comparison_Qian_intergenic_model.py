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

fitness_intergenic=pd.read_excel('datasets/data_wt_rates.xlsx')
fitness_intergenic.index=fitness_intergenic['Standard_name']
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
#column_wt='wt_rates_intergenic_merged_hand'
#column_dnrp1='dnrp1_rates_intergenic_merged_hand'


#norm_wt=fitness_intergenic.loc['HO',column_wt]

#norm_wt_fitness=fitness_intergenic[column_wt]/norm_wt
norm_wt_fitness=fitness_intergenic['rates-intergenic']/fitness_intergenic.loc['HO','rates-intergenic']
#norm_dnrp1_fitness=fitness_intergenic[column_dnrp1]
qian_dataset_SC_fitness=qian_dataset_SC['SC fitness']
#%% Comparing relative fitness scatter


qian_fitness=[]
satay_fitness=[]
for genes in qian_dataset_SC.index:
    if genes in norm_wt_fitness.index:
        
       
        satay_fitness.append(np.mean(norm_wt_fitness[genes]))
        qian_fitness.append(qian_dataset_SC.loc[genes,'SC fitness'])
       
#%% dataframe of the fitness
data=pd.DataFrame([satay_fitness,qian_fitness])
data=data.T
data.columns=['satay-fitness', 'Qian-fitness']
#%%  Plot comparison
p=sns.jointplot('satay-fitness','Qian-fitness',data=data,kind="reg",height=5, ratio=2, marginal_ticks=True,color='black')
p.fig.suptitle("Fitness from Qian et al vs SATAY")

p.fig.tight_layout()
p.fig.subplots_adjust(top=0.95) # Reduce plot to make room        


#%% saving plot

p.savefig('fitness-from-satay-vs-qian.png', format='png',dpi=300)

#%% Comparing relative fitness scatter

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)


for genes in qian_dataset_SC.index:
    if genes in norm_wt_fitness.index:
        
       
        x=norm_wt_fitness[genes]
       
        ax.scatter(x.mean(),y=qian_dataset_SC_fitness[genes],color='blue',alpha=0.4,label='WT')

# ax.set_ylim(0.2,1.2)
# ax.set_xlim(0.2,1.2)
ax.set_ylabel('Qian  SC fitness values of mutants relative to WT ')
ax.set_xlabel('Intergenic model fitness from WT relative values to HO ')


ax.plot(np.linspace(0,1.5),np.linspace(0,1.5),color='black',linewidth=10,alpha=0.2)
#%% saving the figure

fig.savefig('scatter-qian-fitness-satay.png',format='png',dpi=300)
#%% histograms
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(121)

ax.hist(norm_wt_fitness,color='blue')
ax.set_xlabel('fitness-relative-to-HO-satay')
ax2=fig.add_subplot(122)
ax2.hist(qian_dataset_SC_fitness,color='gray')
ax2.set_xlabel('qian-fitness-relative-to-wt')
#%%

fig.savefig('histograms-qian-fitness-satay.png',format='png',dpi=300)

