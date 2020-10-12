# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 11:57:58 2020

@author: linigodelacruz
"""

import numpy as np
import matplotlib.pyplot as plt

from python_modules.module_merging_fitness_with_interactions_and_go_terms import merging_datasets
#%% Calling the module to compute the request

gene_name='BEM1'
data=merging_datasets(gene_name)
score_synt_pd=data
#%% Adjust the type on which interaction type was used 
type_int='All interactions'

#%% Plots for the multiplicative model

fig, axes = plt.subplots(1, 2, figsize=(7,3), dpi=100, sharex=True, sharey=True)
plt.subplots_adjust(right=2,wspace=2)
max=1
plt.subplot(121)
plt.scatter(score_synt_pd['query-fitness']*score_synt_pd['array-fitness'],score_synt_pd['double-fitness'],alpha=0.8, c=score_synt_pd['score'])
plt.plot(np.linspace(0,max),np.linspace(0,max))
plt.xlabel('Multiplicative model , query*array fitness')
plt.ylabel('Double-fitness (relative to WT)')
plt.title( type_int +'_'+ gene_name)
plt.colorbar(label='interaction score SGA')
plt.tight_layout()
plt.ylim([0,max])
plt.xlim([0,max])


plt.subplot(122)
plt.scatter(score_synt_pd['query-fitness']*score_synt_pd['array-fitness'],score_synt_pd['double-fitness'],alpha=0.8, c=score_synt_pd['p-value'],cmap='inferno')
plt.plot(np.linspace(0,max),np.linspace(0,max))
plt.xlabel('Multiplicative model , query*array fitness')
plt.ylabel('Double-fitness (relative to WT)')
plt.title( type_int +'_'+ gene_name)
plt.colorbar(label='p-value SGA')
plt.tight_layout()
plt.ylim([0,max])
plt.xlim([0,max])

#%% saving the figure

plt.savefig('output_images/'+ gene_name + '_data_from_constanzo-check-of-the-scores.png',dpi=300,format='png',transparent=True)

#%% PLot Ocurrence of GO Terms for the specific searched interactors

from collections import Counter

go_term_all=[]


merge_data_numeric=score_synt_pd

for i in np.arange(0,len(merge_data_numeric)):
        go_term_all.append(merge_data_numeric.loc[:,'go-term-filtered-gene'][i])
    
go_new=[]     
for i in go_term_all:
    if i!='gene not found':
        go_new.append(i)
        

out_go = np.concatenate(go_new).ravel()


word_list_go = out_go
counts_go = Counter(word_list_go)
labels_go, values_go = zip(*counts_go.items())

# sort your values in descending order
indSort_go = np.argsort(values_go)[::-1]


# rearrange your data
labels_go = np.array(labels_go)[indSort_go]
values_go = np.array(values_go)[indSort_go]
indexes_go = np.arange(len(labels_go))
bar_width = 0.1

fig, axes = plt.subplots(1, 1, figsize=(5,5), dpi=100, sharex=True, sharey=True)

bins=plt.hist(values_go,alpha=0.6,color='blue')

ticks=[]
a=bins[1]
b=bins[0]
for i in np.arange(0,len(b)):
    
    if b[i]!=0:
        
        ticks.append(i)

plt.xticks(a[ticks],labels_go,rotation=70)
plt.title('Ocurrence of GO-terms for '+ type_int + ' of-' + gene_name)
plt.tight_layout()

#%% saving the figure

fig.savefig('output_images/'+ gene_name + '_go_terms_for' + type_int+'.png',dpi=300,format='png',transparent=True)