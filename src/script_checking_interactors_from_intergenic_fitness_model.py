# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 10:23:15 2021

@author: linigodelacruz
"""


import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import os
import seaborn as sns
import scipy 

#%%
interactions_biogrid=pd.read_excel('datasets/data-BioGrid-Yeast.xlsx')
#%%
data=pd.read_csv('datasets/NRP1_genetic_interactions_filtered_by_Costanzo.txt', delimiter = "\t",header=7)
#%%
fitness_values_wt=pd.read_excel('datasets/data_wt_rates.xlsx')
fitness_values_wt.index=fitness_values_wt['Standard_name']
fitness_values_dnrp1=pd.read_excel('datasets/data_dnrp1_rates.xlsx')
fitness_values_dnrp1.index=fitness_values_dnrp1['Standard_name']
#%% Looking for nrp1 existing interactors
interactions_biogrid=data
#nrp1_interactors=interactions_biogrid[interactions_biogrid['gene-query-name']=='NRP1']
nrp1_interactors=interactions_biogrid[interactions_biogrid['Interactor']=='NRP1']
# nrp1_negative=nrp1_interactors[nrp1_interactors['interaction-type']=='Negative Genetic']
# nrp1_positive=nrp1_interactors[nrp1_interactors['interaction-type']=='Positive Genetic']
nrp1_negative=nrp1_interactors[nrp1_interactors['Assay']=='Negative Genetic']
nrp1_positive=nrp1_interactors[nrp1_interactors['Assay']=='Positive Genetic']
nrp1_negative=pd.unique(nrp1_negative['Interactor.1'])
nrp1_positive=pd.unique(nrp1_positive['Interactor.1'])



#%% Plotting the constanzo SGA scores vs my scores out of fitness values 

column='rates-intergenic'

#score_fitness=fitness(dgenednrp1)-fitness(nrp1)fitness(gene)
cte=fitness_values_wt.loc['NRP1',column]/fitness_values_wt.loc['HO',column]
norm_wt=fitness_values_wt[column]/fitness_values_wt.loc['HO',column]
norm_dnrp1=fitness_values_dnrp1[column]/fitness_values_wt.loc['HO',column]

#%%
fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(111)

# for i in data['Interactor.1']:
#     scores.append(fitness_values.loc[i,'average_dnrp1']/norm_dnrp1-cte*fitness_values.loc[i,'average_wt']/norm_wt**2)
#     scores_sga.append(nrp1_interactors[nrp1_interactors['Interactor.1']==i]['SGA score'].tolist()[0])

fitness_values_wt.fillna(0,inplace=True)
fitness_values_dnrp1.fillna(0,inplace=True)
       
        
true_scores=0
for i in nrp1_positive:
    scores=(fitness_values_dnrp1.loc[i,column]-cte*fitness_values_wt.loc[i,column])#normalized to HO locus from WT 
    scores_sga=(nrp1_interactors[nrp1_interactors['Interactor.1']==i]['SGA score'].tolist()[0])
    ax.scatter(x=scores_sga,y=scores,label='positive by Constanzo',color='green',alpha=0.4)
    if scores>0:
        ax.scatter(x=scores_sga,y=scores,label='positive by Constanzo',color='green')
        ax.text(x=scores_sga,y=scores,s=i,fontsize=7.5,rotation=50)
        true_scores=true_scores+1




for i in nrp1_negative:
    scores=(fitness_values_dnrp1.loc[i,column].mean()-cte*fitness_values_wt.loc[i,column].mean())
    scores_sga=(nrp1_interactors[nrp1_interactors['Interactor.1']==i]['SGA score'].tolist()[0])
    ax.scatter(x=scores_sga,y=scores,label='negative by Constanzo',color='purple',alpha=0.4)
    if scores<0:
        ax.scatter(x=scores_sga,y=scores,label='negative by Constanzo',color='purple')
        ax.text(x=scores_sga,y=scores,s=i,fontsize=7.5,rotation=50)
        true_scores=true_scores+1


ax.set_title('Constanzo interactors')
ax.set_xlabel('Scores_SGA')
ax.grid()
ax.set_ylabel('Scores_intergenic model')
ax.set_ylim(-0.05,0.05)
ax.set_xlim(-0.5,0.2)



fig.savefig('merged_rel_to_HO_scores_intergenic_vs_constanzo_scores_nrp1-Greg.png',format='png',dpi=300,transparent=False)


