#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 24 10:24:08 2021

@author: linigodelacruz
"""

import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import os
import seaborn as sns
import scipy 

from numpy import inf
#%% Reading from a post processed dataset exported to an excel sheet

names_libraries={'wt':'data_wt_filtered_reads_per_tr.xlsx','dnrp1':'data_nrp1_filtered_reads_per_tr.xlsx'}
data_library=[]
for i in names_libraries.keys():
 data_library.append(pd.read_excel('datasets/'+names_libraries[i],index_col='Unnamed: 0'))

datasets=data_library

del data_library,i

#%% Removing the ADE2 and URA3 gene reads and insertions from the population

for i in np.arange(0,len(datasets)):
    datasets[i]=datasets[i][datasets[i].Standard_name != 'ADE2']
    datasets[i]=datasets[i][datasets[i].Standard_name != 'URA3']
#%% Preparing datasets for analyses

data_library_pd=pd.concat(datasets,keys=names_libraries.keys(),sort=True)
data_library_pd.fillna(0,inplace=True)

data_wt=data_library_pd.loc['wt']
data_nrp1=data_library_pd.loc['dnrp1']

data_wt.index=data_wt['Standard_name']
data_nrp1.index=data_nrp1['Standard_name']

#%% Opjinen fitness model , for single knockouts fitness

n_0=1/len(data_wt) # normalized number of reads at the start of the reseed, same to all

normalization_n1=np.sum(data_wt['Nreads'])
d=33 # expansion factor of the population from t=0 to t=90 , OD_0=0.25 OD_90=8.3

factor=d/n_0
factor_1=d/(1-n_0)

w_opjinen=[]
a=np.log(factor*np.mean(data_wt.loc['HO','reads-per-tr'])/normalization_n1)
b=np.log(factor_1*(1-np.mean(data_wt.loc['HO','reads-per-tr'])/normalization_n1))
w_ref=a/b


for i in data_wt.index:
    n1=np.mean(data_wt.loc[i,'reads-per-tr'])/normalization_n1
    value_up=np.log(n1*factor)
    value_down=np.log((1-n1)*factor_1)
    w_opjinen.append(np.abs(value_up/value_down))
    
w_opjinen_norm=w_opjinen/np.abs(w_ref)

#%% Readjusting values
w_opjinen_norm[w_opjinen_norm == inf] = 1

w_opjinen_norm[w_opjinen_norm > 2] = 1

#%%comparing with intergenic model 

fitness_intergenic=pd.read_excel('datasets/data_wt_rates.xlsx')

fitness_intergenic.index=fitness_intergenic['Standard_name']

#%%

fitness_satay=fitness_intergenic['rates-intergenic']/fitness_intergenic.loc['HO','rates-intergenic'].tolist()

#%%
data=pd.DataFrame([fitness_satay.tolist(),w_opjinen_norm])
data=data.T
data.columns=['intergenic-fitness', 'Opjinen-fitness']
#%%  Plot comparison
p=sns.jointplot('intergenic-fitness','Opjinen-fitness',data=data,kind="reg",height=5, ratio=2, marginal_ticks=True,color='black')
p.fig.suptitle("Fitness from Opjinen et al vs intergenic model")

p.fig.tight_layout()
p.fig.subplots_adjust(top=0.95) # Reduce plot to make room       

#%%
p.savefig('fitness-from-intergenic-vs-opjinen.png', format='png',dpi=300)