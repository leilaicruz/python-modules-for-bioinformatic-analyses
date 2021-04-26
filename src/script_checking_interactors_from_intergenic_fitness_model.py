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
fitness_values=pd.read_excel('datasets/fitness-from-intergenic-model-greg-sequencing.xlsx')
fitness_values.index=fitness_values['gene_name']
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
#%% Looking into the fitness values by the intergenic model

nrp1_intergenic_fitness_positive=[]
nrp1_intergenic_pvalue_positive=[]
for i in nrp1_positive:
    nrp1_intergenic_fitness_positive.append(fitness_values.loc[i,'fc_log2'])
    nrp1_intergenic_pvalue_positive.append(fitness_values.loc[i,'-log_p'])

nrp1_intergenic_fitness_negative=[]
nrp1_intergenic_pvalue_negative=[]    
for j in nrp1_negative:
    nrp1_intergenic_fitness_negative.append(fitness_values.loc[j,'fc_log2'])
    nrp1_intergenic_pvalue_negative.append(fitness_values.loc[j,'-log_p'])
    

fig = plt.figure(figsize=(12,5))
ax = fig.add_subplot(121)
ax2=fig.add_subplot(122)


ax.scatter(x=nrp1_positive,y=nrp1_intergenic_fitness_positive)
ax.set_title('Positive Genetic')
ax2.scatter(x=nrp1_negative,y=nrp1_intergenic_fitness_negative)
ax2.set_title('Negative Genetic')


#%%
fig = plt.figure(figsize=(20,5))
ax = fig.add_subplot(121)
ax2=fig.add_subplot(122)
# Generate a custom diverging colormap
from matplotlib import cm

#cmap = sns.diverging_palette(133, 10, as_cmap=True)
cmap=cm.PRGn
sns.heatmap([nrp1_intergenic_fitness_positive,nrp1_intergenic_pvalue_positive],vmin=-2,vmax=2,cmap=cmap,
            annot=True,xticklabels=nrp1_positive,yticklabels=['log2FC','-log10p'],ax=ax)
sns.heatmap([nrp1_intergenic_fitness_negative,nrp1_intergenic_pvalue_negative],vmin=-2,vmax=2,cmap=cmap,
            annot=True,xticklabels=nrp1_negative,yticklabels=['log2FC','-log10p'],ax=ax2)


    
ax.set_xlabel('Existing positive interactors (expecting negative log2FC)')
ax2.set_xlabel('Existing negative interactors (expecting positive log2FC)')

#fig.savefig('heatmap-existing-interactors-vs-log2fc.png',format='png',dpi=300,transparent=False)
#%% saving the figure

fig.savefig('heatmap-existing-interactors-vs-log2fc.png',format='png',dpi=300,transparent=False)

#%% Plotting the constanzo SGA scores vs my scores out of fitness values 

column_wt='wt_rates_intergenic_merged_hand'
column_dnrp1='dnrp1_rates_intergenic_merged_hand'
#score_fitness=fitness(dgenednrp1)-fitness(nrp1)fitness(gene)
cte=fitness_values.loc['NRP1',column_wt]
norm_wt=fitness_values.loc['HO',column_wt]
norm_dnrp1=fitness_values.loc['HO',column_dnrp1]


fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(111)

# for i in data['Interactor.1']:
#     scores.append(fitness_values.loc[i,'average_dnrp1']/norm_dnrp1-cte*fitness_values.loc[i,'average_wt']/norm_wt**2)
#     scores_sga.append(nrp1_interactors[nrp1_interactors['Interactor.1']==i]['SGA score'].tolist()[0])

       
        
true_scores=0
for i in nrp1_positive:
    scores=(fitness_values.loc[i,column_dnrp1]/norm_wt-cte*fitness_values.loc[i,column_wt]/norm_wt**2)#normalized to HO locus from WT 
    scores_sga=(nrp1_interactors[nrp1_interactors['Interactor.1']==i]['SGA score'].tolist()[0])
    ax.scatter(x=scores_sga,y=scores,label='positive by Constanzo',color='green',alpha=0.4)
    if scores>0:
        ax.scatter(x=scores_sga,y=scores,label='positive by Constanzo',color='green')
        ax.text(x=scores_sga,y=scores,s=i,fontsize=7.5,rotation=50)
        true_scores=true_scores+1




for i in nrp1_negative:
    scores=(fitness_values.loc[i,'average_dnrp1']/norm_dnrp1-cte*fitness_values.loc[i,'average_wt']/norm_wt**2)
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
ax.set_ylim(-0.6,0.35)
ax.set_xlim(-0.6,0.35)



fig.savefig('merged_rel_to_HO_scores_intergenic_vs_constanzo_scores_nrp1-Greg.png',format='png',dpi=300,transparent=False)


