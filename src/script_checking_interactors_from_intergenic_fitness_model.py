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
fitness_values=pd.read_excel('datasets/fitness-from-intergenic-model-agnes-sequencing.xlsx')
fitness_values.index=fitness_values['gene_name']
#%% Looking for nrp1 existing interactors

nrp1_interactors=interactions_biogrid[interactions_biogrid['gene-query-name']=='NRP1']
nrp1_negative=nrp1_interactors[nrp1_interactors['interaction-type']=='Negative Genetic']
nrp1_positive=nrp1_interactors[nrp1_interactors['interaction-type']=='Positive Genetic']

#%% Looking into the fitness values by the intergenic model

nrp1_intergenic_fitness_positive=[]
nrp1_intergenic_pvalue_positive=[]
for i in nrp1_positive['gene-target-name']:
    nrp1_intergenic_fitness_positive.append(fitness_values.loc[i,'fc_log2'])
    nrp1_intergenic_pvalue_positive.append(fitness_values.loc[i,'-log_p'])

nrp1_intergenic_fitness_negative=[]
nrp1_intergenic_pvalue_negative=[]    
for j in nrp1_negative['gene-target-name']:
    nrp1_intergenic_fitness_negative.append(fitness_values.loc[j,'fc_log2'])
    nrp1_intergenic_pvalue_negative.append(fitness_values.loc[j,'-log_p'])
    

fig = plt.figure(figsize=(12,5))
ax = fig.add_subplot(121)
ax2=fig.add_subplot(122)


ax.scatter(x=nrp1_positive['gene-target-name'],y=nrp1_intergenic_fitness_positive)
ax.set_title('Positive Genetic')
ax2.scatter(x=nrp1_negative['gene-target-name'],y=nrp1_intergenic_fitness_negative)
ax2.set_title('Negative Genetic')


#%%
fig = plt.figure(figsize=(12,5))
ax = fig.add_subplot(121)
ax2=fig.add_subplot(122)
# Generate a custom diverging colormap
from matplotlib import cm

#cmap = sns.diverging_palette(133, 10, as_cmap=True)
cmap=cm.PRGn
sns.heatmap([nrp1_intergenic_fitness_positive,nrp1_intergenic_pvalue_positive],vmin=-2,vmax=2,cmap=cmap,
            annot=True,xticklabels=nrp1_positive['gene-target-name'],yticklabels=['log2FC','-log10p'],ax=ax)
sns.heatmap([nrp1_intergenic_fitness_negative,nrp1_intergenic_pvalue_negative],vmin=-2,vmax=2,cmap=cmap,
            annot=True,xticklabels=nrp1_negative['gene-target-name'],yticklabels=['log2FC','-log10p'],ax=ax2)


    
ax.set_xlabel('Existing positive interactors')
ax2.set_xlabel('Existing negative interactors')

fig.savefig('heatmap-existing-interactors-vs-log2fc.png',format='png',dpi=300,transparent=False)