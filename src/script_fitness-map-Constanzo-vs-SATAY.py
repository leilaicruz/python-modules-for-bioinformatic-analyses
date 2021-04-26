# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 10:30:46 2020

@author: linigodelacruz
"""
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from tkinter import Tk
from tkinter import filedialog

#%% Calling the module to compute the request

gene_name='DPL1'


#%%
# Data from fitness Constanzo 

root = Tk()
root.filename =  filedialog.askopenfilename(title = "choose constanzo fitness dataset ",filetypes = (("excel files","*.xlsx"),("all files","*.*")))
filename_fitness_constanzo=root.filename
root.withdraw()

# Data from fitness SATAY
root = Tk()
root.filename =  filedialog.askopenfilename(title = "choose SATAY fitness dataset ",filetypes = (("excel files","*.xlsx"),("all files","*.*")))
filename_fitness_satay=root.filename
root.withdraw()

datawithfitness_constanzo=pd.read_excel(filename_fitness_constanzo)
datawithfitness_constanzo.columns=['query-allele-name','array-allele-name','score','p-value','query-fitness','array-fitness','double-fitness','double-fitness-std']

#%%
datawithfitness_satay=pd.read_excel(filename_fitness_satay,index_col='Unnamed: 0')
datawithfitness_satay.columns=['query-allele-name','array-allele-name','query-fitness','array-fitness','double-fitness','score']
#%%
#%%
dpl1_fitness=datawithfitness_satay[datawithfitness_satay['array-allele-name']=='DPL1']['array-fitness'].tolist()[0]
# %% Fitness plot SATAY
interactions_pd=datawithfitness_satay
fig, axes=plt.subplots(1,2)
plt.subplots_adjust(right=1,wspace=0.6)
max=10
axes[0].scatter(x=interactions_pd['array-fitness'],y=interactions_pd['double-fitness'],alpha=0.4)
axes[0].set_xlim([0,max])
axes[0].set_xlabel('Single mutant fitness-b')
axes[0].set_ylabel('double mutant fitness-ab')
axes[0].set_ylim([0,max])
x = np.linspace(0, max)
axes[0].plot(x, dpl1_fitness*x,linestyle='solid',color='black');
x_masking = np.linspace(0, dpl1_fitness)
axes[0].plot(x, x,linestyle='solid',color='black');

trianglex=[0,1,dpl1_fitness,0]
triangley=[0,dpl1_fitness,dpl1_fitness,0]

axes[0].fill(trianglex, triangley,alpha=0.5,color='gray')

axes[0].set_title('satay-dpl1/HO')
#plt.savefig('dpl1_data_using_reads_normalized_with_HO_as_fitness.png',format='png',dpi=300,transparent=True)

# Fitness plot Constanzo
interactions_pd=datawithfitness_constanzo[datawithfitness_constanzo['array-allele-name']=='dpl1']
#interactions_pd=datawithfitness_constanzo[datawithfitness_constanzo['query-allele-name']=='bem1']
interactions_pd.index=np.arange(0,len(interactions_pd))

## Redefinig huge values of fitness (>10000)
threshold=1
for i in np.arange(0,len(interactions_pd)):
    if interactions_pd.loc[i,'query-fitness'] > threshold:
        interactions_pd.loc[i,'query-fitness']=threshold
    if interactions_pd.loc[i,'array-fitness'] > threshold:
        interactions_pd.loc[i,'array-fitness']=threshold
    if interactions_pd.loc[i,'double-fitness'] > threshold:
        interactions_pd.loc[i,'double-fitness']=threshold



axes[1].scatter(x=interactions_pd['query-fitness'],y=interactions_pd['double-fitness'],alpha=0.4)
axes[1].set_xlim([0,threshold])
axes[1].set_xlabel('Single mutant fitness-b')
axes[1].set_ylabel('double mutant fitness-ab')
axes[1].set_ylim([0,threshold])
x = np.linspace(0, threshold)
axes[1].plot(x, interactions_pd.loc[0,'array-fitness']*x,linestyle='solid',color='black');

axes[1].plot(x, x,linestyle='solid',color='black');

trianglex=[0,threshold,interactions_pd.loc[0,'array-fitness'],0]
triangley=[0,(threshold)*interactions_pd.loc[0,'array-fitness'],interactions_pd.loc[0,'array-fitness'],0]

axes[1].fill(trianglex, triangley,alpha=0.2,color='gray')

axes[1].set_title('Constanzo' )

#%% save figure

fig.savefig('output_images/constanzo-vs-satay-dpl1-fitness-map.png',format='png',dpi=300,transparent=True)