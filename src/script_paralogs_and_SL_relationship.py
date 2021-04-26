# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 15:10:06 2020

@author: linigodelacruz
"""

import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


#%%

data_sl=pd.read_excel('../datasets/data-synthetic-lethals.xlsx',header=0)
all_paralogs_from_sl=pd.read_excel('../datasets/paralogs-all-unique-SL-pairs.xlsx')
query_paralogs_pd=all_paralogs_from_sl.drop(columns='Unnamed: 0')
query_paralogs_pd.columns=['name-gene','name-paralogue']
query_paralogs_pd_withoutnan=query_paralogs_pd.dropna()
query_paralogs_pd_withoutnan.index=np.arange(0,len(query_paralogs_pd_withoutnan))
query_paralogs_pd=query_paralogs_pd_withoutnan

del query_paralogs_pd_withoutnan,all_paralogs_from_sl

#%% Search inside the all paralogs list which pair are also synthetic lethals
indexes_sl_query=[]

for i in np.arange(0,len(query_paralogs_pd)):
    paralog_target=query_paralogs_pd[query_paralogs_pd['name-gene']==query_paralogs_pd['name-gene'][i]]['name-paralogue'].tolist()[0]
    list_targets_sl=data_sl[data_sl['gene-query-name']==query_paralogs_pd['name-gene'][i]]['gene-target-name'].tolist()
    
    if paralog_target in list_targets_sl:
        indexes_sl_query.append(query_paralogs_pd[query_paralogs_pd['name-paralogue']==paralog_target].index[0])

## Putting 1's if the paralog is also SL
sL_values=np.zeros_like(query_paralogs_pd['name-gene'])
for i in np.arange(0,len(query_paralogs_pd)):
    if i in indexes_sl_query:
        sL_values[i]=1
query_paralogs_pd['sL']=sL_values

#%% Assigning a new dataframe 
paralogs_sl_pd=query_paralogs_pd
del query_paralogs_pd

sl_that_are_paralogs=paralogs_sl_pd[paralogs_sl_pd['sL']==1]
sl_that_are_paralogs.set_index(np.arange(0,len(sl_that_are_paralogs)))


#%% What is the contribution of paralogs to SL that share protein domains?

pairs_sL=np.load('../datasets/pairs-sL-that-share-domains.npy')

shared_sL_paralogs=[]
for i in np.arange(0,len(sl_that_are_paralogs)):
    for j in np.arange(0,len(pairs_sL)):
        if set(sl_that_are_paralogs.iloc[i,0:2].tolist())==set(pairs_sL[j]):
            shared_sL_paralogs.append(pairs_sL[j])

print('The contribution of paralogs to the SL pairs that shared domains is =', 100*len(shared_sL_paralogs)/len(pairs_sL),'%')
print('Only',len(shared_sL_paralogs),'out of',len(sl_that_are_paralogs),'paralogs that are SL , share annotated protein domains.')

print('The number of paralogs that are SL out of the total number of paralogs is',len(sl_that_are_paralogs),'out of',len(paralogs_sl_pd),'=',100*len(sl_that_are_paralogs)/len(paralogs_sl_pd),'%')

print('The contribution of paralogs to the total number of SL pairs is  =', 100*len(sl_that_are_paralogs)/17871,'%')
print('The number of SL that share domains out of the total number of SL pairs is =',100*len(pairs_sL)/17871,'%')

#%%  Vizualization proportion of paralogs that are also SL 

labels = ['Paralogs SL','Paralogs non SL']
sizes = [len(sl_that_are_paralogs)/len(paralogs_sl_pd),1-len(sl_that_are_paralogs)/len(paralogs_sl_pd)]

colors = ['#00F28E','#F20064']
 
fig1, ax1 = plt.subplots()
patches, texts, autotexts = ax1.pie(sizes, colors = colors, labels=labels, autopct='%1.1f%%', startangle=90)
for text in texts:
    text.set_color('black')
for autotext in autotexts:
    autotext.set_color('black')
# Equal aspect ratio ensures that pie is drawn as a circle
ax1.axis('equal')  
plt.tight_layout()

#%% save figure

fig1.savefig('../output_images/one-quarter-of-paralogs-are-SL.png',format='png',dpi=300,transparent=True)