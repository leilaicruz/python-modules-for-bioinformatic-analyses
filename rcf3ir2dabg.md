```
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 10:33:19 2020
@author: linigodelacruz
"""
#%% Libraries 

import numpy as np

import seaborn as sns

import pandas as pd


import matplotlib.pyplot as plt
from python_modules.module_essential_and_non_essential_genes_interactors import how_many_interactors

#%%
## Importing datasets
essential_list_2= np.loadtxt(r'../datasets/Cervisiae_EssentialGenes_List_2.txt',skiprows=3,dtype='U4')
data_interactions=pd.read_excel(r'../datasets/data-BioGrid-Yeast.xlsx')

# %% Preparing the datasets
essential_list_pd=pd.DataFrame(essential_list_2)
essential_list_pd.columns=['essential-genes']

data_interactions['essential-genes']=essential_list_pd
essential_list_array=np.array(essential_list_pd)

#%% Executing the function

output=how_many_interactors(data_interactions,essential_list_array)
#%% saving to excel
output.to_excel('../datasets/interactors-of-esential-and-non-essential-genes.xlsx')

#%% Vizualization
interactors=output
sns.distplot(interactors.loc['non-essentials'],label='non-essential')
sns.distplot(interactors.loc['essentials'],label='essential')
plt.legend()
plt.xlabel('Number of total interactors')
plt.ylabel('normalized density')
plt.savefig('output_images/essential-and-not-essential-genes-number-of-interactors.png',format='png',dpi=300,transparent=True)

```