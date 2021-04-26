# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 08:55:02 2020

@author: linigodelacruz
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict 
import seaborn as sns
import matplotlib.cm as cm
import scipy as scipy
import random
import os


def sample_protein_pairs(data_domains,data_sl,data_nonsl,sample_size):
    
    """
    This function randomly samples over the population of protein pairs,
    the number of pairs is wished to analyze in the ML pipeline.
    input: 
        data_domains: dataframe with the list of domains per protein , their description , position etc from P-famA 
        data_sl: dataframe with the SL gene pairs
        data_nsl: dataframe with the nSL gene pairs , in this case we can use positive genetic genetic pairs 
    output:
        list of the SL pairs with their respective domains(a nd b)
        list of the nSL pairs with their respective domains (non a and non b)
    
    """
    # Selecting the meaningful columns in the respective dataset
    query_gene=data_sl['gene-query-name']
    target_gene=data_sl['gene-target-name']
    query_gene_nonlethal=data_nonsl['gene-query-name']
    target_gene_nonlethal=data_nonsl['gene-target-name']
    
    
    
    # Initialising the arrays
    protein_a_list=[]
    protein_b_list=[]
    protein_a_list_non=[]
    protein_b_list_non=[]
    
    population = np.arange(0,len(data_sl))
    
    # For loop for 10000 pairs sampled randomly from the SL/nSl pair list , and creating a big array of protein domains id per protein pair
    random.seed() # initialize the random number generation
    for m in random.sample(list(population), sample_size):
        protein_a=data_domains[data_domains['name']==query_gene[m]]
        protein_b=data_domains[data_domains['name']==target_gene[m]]
        protein_a_list.append(protein_a['domain-name'].tolist())
        protein_b_list.append(protein_b['domain-name'].tolist())
    
        protein_a_non=data_domains[data_domains['name']==query_gene_nonlethal[m]]
        protein_b_non=data_domains[data_domains['name']==target_gene_nonlethal[m]]
        protein_a_list_non.append(protein_a_non['domain-name'].tolist())
        protein_b_list_non.append(protein_b_non['domain-name'].tolist())

    return protein_a_list,protein_a_list_non,protein_b_list,protein_b_list_non

def remove_empty_domains(protein_list_search,protein_list_pair):
  
    """
    Function that remove empty domains from the protein pair . 
    Parameters
    ----------
        protein_list_search- list
        one of both protein pairs from the same type of interaction
        protein_list_pair- list
        the other pair list
        
    Returns
    -------
        protein_list_search_new: updated protein pairs without empty domains
        protein_list_pair_new: updated protein pairs without empty domains
    """

    index=[]
    for i in np.arange(0,len(protein_list_search)):
        if protein_list_search[i]==[] or protein_list_pair[i]==[]:
            index.append(i) ## index of empty values for the protein_a_list meaning they dont have any annotated domain

    y=[x for x in np.arange(0,len(protein_list_search)) if x not in index] # a list with non empty values from protein_a list

    protein_list_search_new=[]
    protein_list_pair_new=[]
    for i in y:
        protein_list_search_new.append(protein_list_search[i])
        protein_list_pair_new.append(protein_list_pair[i])
    return protein_list_search_new,protein_list_pair_new

def feature_building(protein_a_list_new,protein_b_list_new,domain_id_list):
    """
    Function that builds the features from which the ML will learn from.
    In this case is based on whether a protein domain is shared or not between the proteins pairs . 
    Example: if protein domain A is shared in the pair then the position of that domain in the list of domains will have a 2.
        if domain A is only found in one of the pairs then , it will have a 1 in the list of domains.
        if domain A is not found in any of the pairs then it will have a 0 in th list of domains. 

    Parameters
    ----------
    protein_a_list_new : list
        the list of one of the pairs protein domains WITHOUT EMPTY ONES.
    protein_b_list_new : list
        the list of the other pair of protein domains WITHOUT EMPTY ONES.
    domain_id_list : list,array , series
        List of all domains 

    Returns
    -------
    protein_feat_pd : dataframe
        dataframe with size protein-pairs X domain list size , full of 0 , 1 or 2 accordingly. 

    """
    x = np.unique(domain_id_list)
    ## To avoid taking repeated domains from one protein of the pairs , lets reduced the domains of each protein from the pairs to their unique members
    protein_a_list_unique=[]
    protein_b_list_unique=[]
    for i in np.arange(0,len(protein_a_list_new)):
        protein_a_list_unique.append(np.unique(protein_a_list_new[i]))
        protein_b_list_unique.append(np.unique(protein_b_list_new[i]))
        
    protein_feat=np.zeros(shape=(len(x),len(protein_a_list_unique)))
    pair_a_b_array=[]
    get_indexes = lambda x, xs: [i for (y, i) in zip(xs, range(len(xs))) if x == y] # a function that give the index of whether a value appear in array or not
    for i in np.arange(0,len(protein_a_list_unique)):
                
        pair=[protein_a_list_unique[i],protein_b_list_unique[i]]
        pair_a_b=np.concatenate(pair).ravel()
        pair_a_b_array.append(pair_a_b)

    j=0
    for i in pair_a_b_array:  
        array,index,counts=np.unique(i,return_index=True,return_counts=True)
        
        for k,m in zip(counts,array):
            if k ==2:
                protein_feat[get_indexes(m,x),j]=2
                
            if k==1:
                protein_feat[get_indexes(m,x),j]=1
        j=j+1
    protein_feat_pd=pd.DataFrame(protein_feat.T)
    return protein_feat_pd