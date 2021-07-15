# -*- coding: utf-8 -*-

# import libraries
#from combat.pycombat import pycombat
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#%% EXAMPLe with benoit data

dataset_1 = pd.read_csv("datasets/Processed_dataset - WildType1_pergene.txt",sep='\t') # datasets can also be stored in csv, tsv, etc files
dataset_2 = pd.read_csv("datasets/Processed_dataset - WildType2_pergene.txt",sep='\t')

for datasets in [dataset_1,dataset_2]:
    datasets.drop(columns='number_of_read_pergene',inplace=True)
    datasets.columns=['tr-per-gene', 'reads-per-gene']
    
#%%
# we merge all the datasets into one, by keeping the common genes only
df_expression = pd.concat([dataset_1,dataset_2],axis=1)

# plot raw data
#plt.boxplot(df_expression.transpose())

#%%
# we generate the list of batches
batch = []

datasets = [dataset_1,dataset_2]

for j in range(len(datasets)):
    batch.extend([j for _ in range(len(datasets[j].columns))])

df_expression.fillna(1,inplace=True)
df_expression.replace(0,1,inplace=True)

# run pyComBat
df_corrected = pycombat(df_expression,batch)