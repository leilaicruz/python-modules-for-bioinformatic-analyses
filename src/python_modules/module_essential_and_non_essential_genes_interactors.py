# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 11:28:00 2020

@author: linigodelacruz

"""

import numpy as np
from collections import defaultdict 
import pandas as pd

def how_many_interactors(data_of_interactors,who,excluding_type):
    
    """
    function to compute how many interactors have the essential genes
    compared with non essential genes
    
    input: data_of_interactors=dataframe with all the genes and their interactors from budding yeast 
           who=  array containing the essential genes in WT
           excluding_type= Type of interactions you will like to exclude from the analysis, for example "physical"
    output : dataframe with the total number of interactors for the essential and non essential genes 
    
    """
    
    total_interactors=defaultdict(dict)
    
    ## excluding physical interactions
    
    if excluding_type==None:

        for i in np.arange(0,len(who)):
          
            tmp_target=data_of_interactors[data_of_interactors['gene-target-name']==who[i][0]]['gene-query-name'].unique()
            tmp_query=data_of_interactors[data_of_interactors['gene-query-name']==who[i][0]]['gene-target-name'].unique()
            tmp=np.unique(np.concatenate([tmp_target,tmp_query]))
            total_interactors['essentials',who[i][0]]['total of interactors']=len(tmp)
        
        unique_list_non_essential_genes=np.unique(np.concatenate([data_of_interactors['gene-query-name'].unique(),data_of_interactors['gene-target-name'].unique()]))
        for k in np.arange(0,len(unique_list_non_essential_genes)):
            
            if unique_list_non_essential_genes[k] not in who:
                tmp_query=data_of_interactors[data_of_interactors['gene-query-name']==unique_list_non_essential_genes[k]]['gene-target-name'].unique()
                tmp_target=data_of_interactors[data_of_interactors['gene-target-name']==unique_list_non_essential_genes[k]]['gene-query-name'].unique()
                tmp=np.unique(np.concatenate([tmp_target,tmp_query]))
                total_interactors['non-essentials',unique_list_non_essential_genes[k]]['total of interactors']=len(tmp)

    else:
        filtered_data=data_of_interactors[data_of_interactors['interaction-type']!=excluding_type]
        for i in np.arange(0,len(who)):
            
            tmp_target=filtered_data[filtered_data['gene-target-name']==who[i][0]]['gene-query-name'].unique()
            tmp_query=filtered_data[filtered_data['gene-query-name']==who[i][0]]['gene-target-name'].unique()
            tmp=np.unique(np.concatenate([tmp_target,tmp_query]))
            total_interactors['essentials',who[i][0]]['total of interactors']=len(tmp)
    
        unique_list_non_essential_genes=np.unique(np.concatenate([filtered_data['gene-query-name'].unique(),filtered_data['gene-target-name'].unique()]))
        for k in np.arange(0,len(unique_list_non_essential_genes)):
            
            if unique_list_non_essential_genes[k] not in who:
                tmp_query=filtered_data[filtered_data['gene-query-name']==unique_list_non_essential_genes[k]]['gene-target-name'].unique()
                tmp_target=filtered_data[filtered_data['gene-target-name']==unique_list_non_essential_genes[k]]['gene-query-name'].unique()
                tmp=np.unique(np.concatenate([tmp_target,tmp_query]))
                total_interactors['non-essentials',unique_list_non_essential_genes[k]]['total of interactors']=len(tmp)



    total_interactors_pd=pd.DataFrame(total_interactors)

    return total_interactors_pd.T