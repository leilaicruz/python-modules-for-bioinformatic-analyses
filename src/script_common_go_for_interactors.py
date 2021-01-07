# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 15:51:35 2020

@author: linigodelacruz
"""
import numpy as np
from collections import defaultdict 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import norm
import os
os.chdir('D:/Users/Thomas/Studie/MEP/python-modules-for-bioinformatic-analyses/src')
from python_modules.module_common_measures import common_go

from python_modules.module_common_measures import common_partners

#%% getting the data

data_go=pd.read_excel('../datasets/slim-goterms-filtered-data.xlsx')
# data_go=pd.read_excel('../datasets/testDataGo.xlsx')

data_go.columns=['Gene','gene-id','go-aspect','go-term','go-id','feature-type']

data=pd.read_excel('../datasets/data-BioGrid-Yeast-doubled.xlsx')
# data=pd.read_excel('../datasets/testDataInteractions.xlsx')
#%% query
query = pd.read_excel('../datasets/genesCellPolarity_SGD_amigo.xlsx')
query = np.unique(query).tolist()
# query=['BEM2','BEM3']
# query=['geneA','GeneB']

#%% Calling the function common_partners

common_partners_data=common_partners(query,data)

common_go=common_go(data_go=data_go,data_common_partners=common_partners_data)


#%% Postprocessing the data

sl=common_partners_data[common_partners_data['Type']=='Synthetic Lethality']
pg=common_partners_data[common_partners_data['Type']=='Positive Genetic']
ng=common_partners_data[common_partners_data['Type']=='Negative Genetic']

common_go.loc[sl.index,'score']='SL'
common_go.loc[pg.index,'score']='PG'
common_go.loc[ng.index,'score']='NG'

common_go.loc[sl.index,'common_interactors']=common_partners_data.loc[sl.index,'fraction-of-common-partners']
common_go.loc[pg.index,'common_interactors']=common_partners_data.loc[pg.index,'fraction-of-common-partners']
common_go.loc[ng.index,'common_interactors']=common_partners_data.loc[ng.index,'fraction-of-common-partners']

#%% Messing around: Correlation coefficient etc.
# covar = np.cov(ng['fraction-of-common-partners'].astype(float),ng['fraction-of-common-partners'].astype(float))
# corcoeff = np.corrcoef(ng['fraction-of-common-partners'].astype(float),ng['fraction-of-common-partners'].astype(float))
# testdata = np.random.normal(loc=5.0, scale=2.0, size=1000)
# mean,std=norm.fit(data)

#%% viz
sns.set(style="ticks", color_codes=True)
plot=sns.pairplot(common_go,hue='score',vars=['fraction-of-common-go','common_interactors'],palette='dark')
# plt.title(query[0])
plt.title('CellPolarity Genes')


#%% Saving the figure

# plot.savefig('../output_images/common-go-terms-of-'+ query[0]+'-based-on-their-type.png',dpi=300,format='png',transparent=True)