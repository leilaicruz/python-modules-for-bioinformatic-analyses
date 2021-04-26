# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 13:01:45 2020

@author: linigodelacruz
"""

#%% Libraries

import numpy as np
from collections import defaultdict 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from python_modules.module_common_measures import common_partners

#%% getting the data

data=pd.read_excel('../datasets/data-BioGrid-Yeast.xlsx')
#%%
query=['CLA4']

#%% Calling the function

common_partners_data=common_partners(query,data)

#%% Postprocessing the data

sl=common_partners_data[common_partners_data['Type']=='Synthetic Lethality']
pg=common_partners_data[common_partners_data['Type']=='Positive Genetic']
ng=common_partners_data[common_partners_data['Type']=='Negative Genetic']

common_partners_data.loc[sl.index,'score']='SL'
common_partners_data.loc[pg.index,'score']='PG'
common_partners_data.loc[ng.index,'score']='NG'

#common_partners_data.fillna(0,inplace=True)

#%% Vizualising the pattern of common_interactors vs Type of interaction

sns.set(style="ticks", color_codes=True)
plot=sns.pairplot(common_partners_data,hue='score',vars=['fraction-of-common-partners','number of partners of pairB'],palette='dark')
plot.fig.suptitle(query[0])
#%% Saving the figure

plot.savefig('../output_images/common-interactors-of-'+ query[0]+'-based-on-their-type.png',dpi=300,format='png',transparent=True)