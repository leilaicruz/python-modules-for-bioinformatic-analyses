# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 10:18:38 2020

@author: linigodelacruz
"""
import numpy as np
from collections import defaultdict 
import pandas as pd
from tkinter import Tk
from tkinter import filedialog
from sys import exit


def merging_datasets(gene_name=None,datawithfitness=None,datawithgoterms=None,datawithinteractions=None):
    """
    Function that allows filtering the fitness dataset from Constanzo 2016, with data from 
    interactions and go terms , genewise. 
    Inputs:
gene_name= gene on which you want to know the data available , it should be introduced in capitals e.g. 'BEM1'
data_interactions= here is a dataframe on the data on interactors for example we would like to merge with the fitness one to filtered  it 
type_of_interaction_selected=string containing the info of which dataset are we using to filter  the fitness data e.g 'synthetic-lethals'
    Output: a dataframe with the merge information 
    """
   
    if gene_name==None:
        gene_name='BEM1' # do this gene as an example
        
    if datawithfitness == None:
    
        root = Tk()
        root.filename =  filedialog.askopenfilename(title = "choose your fitness dataset",filetypes = (("excel files","*.xlsx"),("all files","*.*")))
        filename_fitness=root.filename
        root.withdraw()
        
    if datawithgoterms == None:
    
        root = Tk()
        root.filename =  filedialog.askopenfilename(title = "choose your go terms dataset",filetypes = (("excel files","*.xlsx"),("all files","*.*")))
        filename_go=root.filename
        root.withdraw()
        
    if datawithinteractions== None:
    
        root = Tk()
        root.filename =  filedialog.askopenfilename(title = "choose your interaction type dataset",filetypes = (("excel files","*.xlsx"),("all files","*.*")))
        filename_interactors=root.filename
        root.withdraw()
    
    
    datawithfitness=pd.read_excel(filename_fitness)
    datawithfitness.columns=['query-allele-name','array-allele-name','score','p-value','query-fitness','array-fitness','double-fitness','double-fitness-std']
    
    datawithgoterms=pd.read_excel(filename_go)
    datawithgoterms.columns=['Gene','gene-id','go-aspect','go-term','go-id','feature-type' ]
    
    datawithinteractions=pd.read_excel(filename_interactors)
    datawithinteractions.columns=['Gene', 'Interactor', 'Assay', 'citation']
    
   
    data_fitness_sga=datawithfitness
    data_raw_slim_go=datawithgoterms
    
    
      
    
    merge_data=defaultdict(dict)
    
    interactors_gene_name=datawithinteractions[datawithinteractions['Gene']==gene_name]['Interactor'].unique()
    constanzo_data=datawithfitness[datawithfitness['query-allele-name']==gene_name.casefold()]['array-allele-name'].tolist()
    
    if len(constanzo_data)==0:
        
        constanzo_data=datawithfitness[datawithfitness['array-allele-name']==gene_name.casefold()]['query-allele-name'].tolist()
    
    if len(constanzo_data)==0:
        print('The specified interactors for', gene_name,'are not in the fitness dataset of Constanzo')
        exit()
        
    genes_to_analyze=[]
    for i in interactors_gene_name:
        if i.casefold() in constanzo_data:
            genes_to_analyze.append(i)
            
   
    filtered_genes=genes_to_analyze
    data_fitness_sga_subset=data_fitness_sga[data_fitness_sga['query-allele-name']==gene_name.casefold()]
    
    merge_data=defaultdict(dict)
    for i in np.arange(0,len(filtered_genes)):
        
       
        tmp=data_fitness_sga_subset[data_fitness_sga_subset['array-allele-name']==filtered_genes[i].casefold()]
        tmp.index=np.arange(0,len(tmp))
        if len(tmp)!=0:
               
            if tmp['query-fitness'].all() > 2 or tmp['array-fitness'].all() > 2 or tmp['double-fitness'].all() >2 : ## these are relative growth rates , so this is to avoid large numbers in the datasets, that are non -interpretable data
                tmp['query-fitness']=2
                tmp['array-fitness']=2
                tmp['double-fitness']=2
            
            merge_data['gene_name'][filtered_genes[i]]=gene_name
            merge_data['array_name'][filtered_genes[i]]=gene_name
            merge_data['query-gene'][filtered_genes[i]]=filtered_genes[i]
            merge_data['score'][filtered_genes[i]]=tmp['score'].tolist()[0]
            merge_data['p-value'][filtered_genes[i]]=tmp['p-value'].tolist()[0]
            merge_data['double-fitness'][filtered_genes[i]]=tmp['double-fitness'].tolist()[0]
            merge_data['array-fitness'][filtered_genes[i]]=tmp['array-fitness'].tolist()[0]
            merge_data['query-fitness'][filtered_genes[i]]=tmp['query-fitness'].tolist()[0]
            
    
            data_go=data_raw_slim_go[data_raw_slim_go['Gene']==filtered_genes[i]]
    
        
        else:
            
            tmp=data_fitness_sga_subset[data_fitness_sga_subset['query-allele-name']==filtered_genes[i].casefold()]
            tmp.index=np.arange(0,len(tmp))
            
            if tmp['query-fitness'].all() > 2 or tmp['array-fitness'].all() > 2 or tmp['double-fitness'].all() >2 : ## these are relative growth rates , so this is to avoid large numbers in the datasets, that are non -interpretable data
                tmp['query-fitness']=2
                tmp['array-fitness']=2
                tmp['double-fitness']=2
            
            merge_data['gene_name'][filtered_genes[i]]=gene_name
            merge_data['array_name'][filtered_genes[i]]=gene_name
            merge_data['query-gene'][filtered_genes[i]]=filtered_genes[i]
            merge_data['score'][filtered_genes[i]]=tmp['score'].tolist()[0]
            merge_data['p-value'][filtered_genes[i]]=tmp['p-value'].tolist()[0]
            merge_data['double-fitness'][filtered_genes[i]]=tmp['double-fitness'].tolist()[0]
            merge_data['array-fitness'][filtered_genes[i]]=tmp['array-fitness'].tolist()[0]
            merge_data['query-fitness'][filtered_genes[i]]=tmp['query-fitness'].tolist()[0]
    
        if len(data_go)==0:
            merge_data['go-term-filtered-gene'][filtered_genes[i]]='gene not found'
        else:
            merge_data['go-term-filtered-gene'][filtered_genes[i]]=data_go.iloc[:,3].tolist()
    
    
       
     
    merge_data_pd=pd.DataFrame(merge_data)
    merge_data_numeric=merge_data_pd[pd.to_numeric(merge_data_pd['score'] ,errors='coerce').notnull()]
    
            

    return merge_data_numeric


