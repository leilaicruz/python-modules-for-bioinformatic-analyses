#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 12:01:37 2021

@author: linigodelacruz
"""


import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from collections import defaultdict
import os
import seaborn as sns
import scipy 
from functools import reduce

from src.python_modules.module_analysis_transposon_sites import *
#%% Import of dataframes output from the SATAY pipeline

names_libraries={'wt_a':'data_wt_a.xlsx','wt_b':'data_wt_b.xlsx',
                 'dnrp1_1_a':'dnrp1_1_1_a_unmerged.xlsx','dnrp1_1_b':'dnrp1_1_1_b_unmerged.xlsx',
                 'dnrp1_2_a':'dnrp1_1_2_a_unmerged.xlsx','dnrp1_2_b':'dnrp1_1_2_b_unmerged.xlsx'}
data_library=[]
for i in names_libraries.keys():
 data_library.append(pd.read_excel('datasets/'+names_libraries[i],index_col='Unnamed: 0'))
#%%
#### Creating a big dataframe of the libraries

datasets=data_library

for i in np.arange(0,len(datasets)):
    datasets[i]=datasets[i][datasets[i].Standard_name != 'ADE2']
    datasets[i]=datasets[i][datasets[i].Standard_name != 'URA3']

data_library_pd=pd.concat(datasets,keys=names_libraries.keys(),sort=True)
data_library_pd.fillna(0,inplace=True)

#%%

################# Computing measures ######################

freq=frequency_transposons(data_library_pd,names_libraries)
reads_per_tr=reads_per_transposon(data_library_pd,names_libraries)
tr_density=transposon_density(data_library_pd,names_libraries)
median_insertions=median_feature(data_library_pd,names_libraries,'Ninsertions')
mean_insertions=mean_feature(data_library_pd,names_libraries,'Ninsertions')
mean_reads_per_insrt=mean_feature(data_library_pd,names_libraries,'Nreadsperinsrt')
std_reads_per_insrt=std_feature(data_library_pd,names_libraries,'Nreadsperinsrt')
std_insertions=std_feature(data_library_pd,names_libraries,'Ninsertions')
median_insert_essentials=median_feature_essentials(data_library_pd,names_libraries,'Ninsertions')
median_insert_nonessentials=median_feature_nonessentials(data_library_pd,names_libraries,'Ninsertions')


#%% Assembling the masure into a dataframe 
analysis_libraries=defaultdict(dict)

j=0
for i in names_libraries.keys():
    
    analysis_libraries[i]['one-transposon-per-bp']=freq[j]
    analysis_libraries[i]['mean-insertions']=mean_insertions[j]
    analysis_libraries[i]['std-insertions']=std_insertions[j]
    analysis_libraries[i]['mean-reads-per-tr']=mean_reads_per_insrt[j]
    analysis_libraries[i]['std-reads-per-tr']=std_reads_per_insrt[j]
    analysis_libraries[i]['Ninsertions']=data_library_pd.loc[i]['Ninsertions']
    analysis_libraries[i]['Nreads']=data_library_pd.loc[i]['Nreads']
    analysis_libraries[i]['median-reads-per-tr']=reads_per_tr[j]
    analysis_libraries[i]['median-tr']=median_insertions[j]
    analysis_libraries[i]['median-tr-essentials']=median_insert_essentials[j]
    analysis_libraries[i]['median-tr-non-essentials']=median_insert_nonessentials[j]
    analysis_libraries[i]['tr-density']=data_library_pd.loc[i]['Ninsertions']/data_library_pd.loc[i]['Nbasepairs']
    analysis_libraries[i]['reads-density']=data_library_pd.loc[i]['Nreads']/data_library_pd.loc[i]['Nbasepairs']
    
    j=j+1
    
del i,j
analysis_libraries_pd=pd.DataFrame(analysis_libraries)

#%%  Plot comparison with vectors 
strain_1='dnrp1_1_a'
strain_2='dnrp1_2_a'
variable='Ninsertions'
p=sns.jointplot(x=analysis_libraries_pd.loc[variable,strain_1],y=analysis_libraries_pd.loc[variable,strain_2],data=analysis_libraries_pd,kind="reg",height=5, ratio=3, marginal_ticks=True,color='black')
p.fig.suptitle(strain_1 + " vs "+ strain_2)

p.fig.tight_layout()
p.fig.subplots_adjust(top=0.95) # Reduce plot to make room  
#%%
p.savefig(strain_1 + " vs "+ strain_2 + "_"+variable+'.png',format='png',dpi=300)

#%%
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)
variable='reads-per-tr'
yerr=analysis_libraries_pd.loc['std-reads-per-tr',:]
ax.errorbar(x=analysis_libraries_pd.columns,y=analysis_libraries_pd.loc['median-reads-per-tr',:],yerr=yerr,capsize=10, color='black', 
            elinewidth=2)
ax.set_ylabel(variable)
#%%

 
fig.savefig('mean-and-std-reads-per-tr-all-libraries.png',format='png',dpi=300) \