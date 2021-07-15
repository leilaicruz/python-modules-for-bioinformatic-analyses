# -*- coding: utf-8 -*-
#pyComBat (Behdenna et al, 2020) 
# import libraries
from combat.pycombat import pycombat
import pandas as pd
import matplotlib.pyplot as plt
#%%
# prepare data
# the datasets are dataframes where:
    # the indexes correspond to the gene names
    # the column names correspond to the sample names
# Any number (>=2) of datasets can be treated
dataset_1 = pd.read_pickle("datasets/batch_effects_data_example/GSE18520.pickle") # datasets can also be stored in csv, tsv, etc files
dataset_2 = pd.read_pickle("datasets/batch_effects_data_example/GSE66957.pickle")
dataset_3 = pd.read_pickle("datasets/batch_effects_data_example/GSE69428.pickle")

# we merge all the datasets into one, by keeping the common genes only
df_expression = pd.concat([dataset_1,dataset_2,dataset_3],join="inner",axis=1)

# plot raw data
plt.boxplot(df_expression.transpose())
#%%
# we generate the list of batches
batch = []
datasets = [dataset_1,dataset_2,dataset_3]
for j in range(len(datasets)):
    batch.extend([j for _ in range(len(datasets[j].columns))])

# run pyComBat
df_corrected = pycombat(df_expression,batch)

# visualise results
plt.boxplot(df_corrected.transpose())


