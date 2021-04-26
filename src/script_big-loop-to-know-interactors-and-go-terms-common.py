# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 16:53:29 2020

@author: linigodelacruz
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


#%% getting the data

data=pd.read_excel('../datasets/common-go-terms-and-interactors-of-many-genes.xlsx')

#%% postprocessing of the data
sl=data[data['Type']=='Synthetic Lethality']
pg=data[data['Type']=='Positive Genetic']
ng=data[data['Type']=='Negative Genetic']

data.loc[sl.index,'score']='SL'
data.loc[pg.index,'score']='PG'
data.loc[ng.index,'score']='NG'

#%%viz

sns.set(style="ticks", color_codes=True)
plot=sns.pairplot(data,hue='score',vars=['fraction-of-common-go','fraction-of-common-partners'],palette='dark')

plot.fig.suptitle("No correlation with the common go terms and interactors and the type of GI", y=1.08) # y= some height>1

#%% Saving the figure

plot.savefig('../output_images/no-correlation-based-on-the-common-interact-and-go-terms-type.png',dpi=300,format='png',transparent=True)

