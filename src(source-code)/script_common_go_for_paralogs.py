# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 14:36:59 2020

@author: linigodelacruz
"""
import numpy as np
from collections import defaultdict 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from python_modules.module_common_measures import common_go_paralogs

#%% getting the data

data_go=pd.read_excel('../datasets/slim-goterms-filtered-data.xlsx')
data_paralogs=pd.read_excel('../datasets/paralogs-SL-individual-genes.xlsx')
data_paralogs=data_paralogs.drop(columns=['Unnamed: 0'])
data_paralogs.columns=['query','paralogue-name']
data_go.columns=['Gene','gene-id','go-aspect','go-term','go-id','feature-type' ]

#%% calling the function
paralogs_go=common_go_paralogs(data_go=data_go,data_paralogs=data_paralogs)

#%% viz
fig, axes=plt.subplots(1,1)
sns.set(style="ticks", color_codes=True)
colors = ['#00F28E','#F20064']
sns.distplot(paralogs_go['fraction-of-common-go'],color=colors[1])
plt.ylabel('normalized counts')

plt.title('Evidence that paralogs has functionally diverged')

#%% saving the figure 
#fig.savefig('../output_images/functional-diversification-of-paralogs.png',format='png',dpi=300,transparent=True)