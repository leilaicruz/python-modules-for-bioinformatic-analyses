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
import seaborn as sns

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
root.filename =  filedialog.askopenfilename(title = "choose SATAY WT fitness dataset ",filetypes = (("excel files","*.xlsx"),("all files","*.*")))
filename_fitness_satay=root.filename
root.withdraw()

root = Tk()
root.filename =  filedialog.askopenfilename(title = "choose SATAY mutant fitness dataset ",filetypes = (("excel files","*.xlsx"),("all files","*.*")))
filename_fitness_satay_mutant=root.filename
root.withdraw()

datawithfitness_constanzo=pd.read_excel(filename_fitness_constanzo)
datawithfitness_constanzo.columns=['query-allele-name','array-allele-name','score','p-value','query-fitness','array-fitness','double-fitness','double-fitness-std']

#%%
datawithfitness_satay=pd.read_excel(filename_fitness_satay,index_col='Unnamed: 0')
datawithfitness_satay.index=datawithfitness_satay['Standard_name']
datawithfitness_satay_mutant=pd.read_excel(filename_fitness_satay_mutant,index_col='Unnamed: 0')
datawithfitness_satay_mutant.index=datawithfitness_satay_mutant['Standard_name']

#datawithfitness_satay.columns=['query-allele-name','array-allele-name','query-fitness','array-fitness','double-fitness','score']
#%% Compare fitness SGA vs SATAY
names_query=datawithfitness_constanzo.loc[:,'query-allele-name'].unique()
names_array=datawithfitness_constanzo.loc[:,'array-allele-name'].unique()
ref=datawithfitness_satay.loc['HO','rates-intergenic']
satay_fitness=[]
sga_fitness=[]
for genes in names_query:
    if genes.upper() not in datawithfitness_satay.index:
        satay_fitness.append(0)
        sga_fitness.append(0)
    else: 
        value=datawithfitness_satay.loc[genes.upper(),'rates-intergenic']/ref
        value_sga=datawithfitness_constanzo[datawithfitness_constanzo['query-allele-name']==genes]['query-fitness'].unique()[0]
        if value_sga>2:
            sga_fitness.append(2)
        else:
            sga_fitness.append(value_sga)

        # if value.dtype=='float64':
        #     satay_fitness.append(value)
        # else:
        satay_fitness.append(np.mean(value))
        
for genes in names_array:
    if genes.upper() not in datawithfitness_satay.index:
        satay_fitness.append(0)
        sga_fitness.append(0)
    else: 
        value=datawithfitness_satay.loc[genes.upper(),'rates-intergenic']/ref
        value_sga=datawithfitness_constanzo[datawithfitness_constanzo['array-allele-name']==genes]['array-fitness'].unique()[0]
        if value_sga>2:
            sga_fitness.append(2)
        else:
            sga_fitness.append(value_sga)

        # if value.dtype=='float64':
        #     satay_fitness.append(value)
        # else:
        satay_fitness.append(np.mean(value))
        
#%% dataframe of the fitness
data=pd.DataFrame([satay_fitness,sga_fitness])
data=data.T
data.columns=['satay-fitness', 'sga-fitness']
#%%  Plot comparison
p=sns.jointplot('satay-fitness','sga-fitness',data=data,kind="reg",height=5, ratio=2, marginal_ticks=True,color='black')
p.fig.suptitle("Fitness from SGA vs SATAY")

p.fig.tight_layout()
p.fig.subplots_adjust(top=0.95) # Reduce plot to make room 

#%% saving plot

p.savefig('fitness-from-satay-vs-sga.png', format='png',dpi=300)
#%%
#dpl1_fitness=datawithfitness_satay[datawithfitness_satay['array-allele-name']=='DPL1']['array-fitness'].tolist()[0]
nrp1_fitness=datawithfitness_satay.loc['NRP1','rates-intergenic']/datawithfitness_satay.loc['HO','rates-intergenic']
#%%
interactions_pd=pd.DataFrame()
interactions_pd['array-fitness']=nrp1_fitness*datawithfitness_satay['rates-intergenic']
interactions_pd['double-fitness']=datawithfitness_satay_mutant['rates-intergenic']
# %% Fitness plot SATAY

fig, axes=plt.subplots(1,2)
plt.subplots_adjust(right=1,wspace=0.6)
max=1
axes[0].scatter(x=interactions_pd['array-fitness'],y=interactions_pd['double-fitness'],alpha=0.4)
axes[0].set_xlim([0,max])
axes[0].set_xlabel('Single mutant fitness-b')
axes[0].set_ylabel('double mutant fitness-ab')
axes[0].set_ylim([0,max])
x = np.linspace(0, max)
axes[0].plot(x, nrp1_fitness*x,linestyle='solid',color='black');
x_masking = np.linspace(0, nrp1_fitness)
axes[0].plot(x, x,linestyle='solid',color='black');

trianglex=[0,1,nrp1_fitness,0]
triangley=[0,nrp1_fitness,nrp1_fitness,0]

axes[0].fill(trianglex, triangley,alpha=0.5,color='gray')

axes[0].set_title('satay-dpl1/HO')
#plt.savefig('dpl1_data_using_reads_normalized_with_HO_as_fitness.png',format='png',dpi=300,transparent=True)

## FIX THIS!!
# Fitness plot Constanzo
interactions_pd_constanzo=datawithfitness_constanzo[datawithfitness_constanzo['array-allele-name']=='dpl1']
#interactions_pd=datawithfitness_constanzo[datawithfitness_constanzo['query-allele-name']=='bem1']
interactions_pd_constanzo.index=np.arange(0,len(interactions_pd_constanzo))

## Redefinig huge values of fitness (>10000)
threshold=1
for i in np.arange(0,len(interactions_pd_constanzo)):
    if interactions_pd_constanzo.loc[i,'query-fitness'] > threshold:
        interactions_pd_constanzo.loc[i,'query-fitness']=threshold
    if interactions_pd_constanzo.loc[i,'array-fitness'] > threshold:
        interactions_pd_constanzo.loc[i,'array-fitness']=threshold
    if interactions_pd_constanzo.loc[i,'double-fitness'] > threshold:
        interactions_pd_constanzo.loc[i,'double-fitness']=threshold



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