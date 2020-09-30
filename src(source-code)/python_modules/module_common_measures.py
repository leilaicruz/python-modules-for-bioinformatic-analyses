# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 12:30:12 2020

@author: linigodelacruz
"""
import numpy as np
from collections import defaultdict 
import pandas as pd



## Function to find common (how many) interaction partners of one gene with its interactors


def common_partners(query,data):
    
    """
    A function that finds the common interactors partners between the query genes provided from the data provided.
    query: dataframe with the target genes to look at
    data:dataframe with the info to take from
    """
        
    d2 = defaultdict(dict)
   
    # giant for loop
    
    for genes_names in query:
        #filtering the table just for the value of the query
        q1 = data[data['gene-query-name']==genes_names]
        q1_interact=q1['gene-target-name'].unique()
        
            
        for query2 in q1_interact:

        
            q2=data[data['gene-query-name']==query2] #these are get_query(q1[i])

            q2_interact=q2['gene-target-name'].unique()
        

            d = defaultdict(int)
            common = []
            for genes_names2  in q2_interact:
                if genes_names2 in q1_interact: # if a gene interactor of the query1 is in interactors of query 2 
                    common.append(genes_names2)
                    # tmp=np.array(q1[q1['gene-target-name']==genes_names2]['interaction-type'].tolist())
                    # type_interaction.append(tmp.ravel().ravel().ravel())
                    d[genes_names2] += 1
                    
            d2[query2]['query']=genes_names

            tmp=q1[q1['gene-target-name']==query2]['interaction-type'].tolist()
            #tmp1.append(tmp)
            tmp=np.unique(np.array(tmp))
           
            if "Synthetic Lethality" in tmp:
                d2[query2]["Type"] = "Synthetic Lethality"
                        
            else:
                d2[query2]["Type"] = tmp[0]

            d2[query2]["common"] = common
            d2[query2]["names of genes"]=query2
            d2[query2]["n_common"] = len(common)
            d2[query2]["number of partners of pairA"] = len(q1_interact)
            d2[query2]["number of partners of pairB"] = len(q2_interact)

            
            
            if len(q1)==0 :
                d2[query2]["fraction-of-common-partners"] = 0
            else:
                d2[query2]["fraction-of-common-partners"] = len(d)/len(q1_interact) *100     
       

        df=pd.DataFrame(d2).T
        #df.set_index=query2
        if len(df)==0:
            df_sorted=[]
        else:
            df_sorted=df.sort_values(by=["fraction-of-common-partners"])
            df_sorted=df_sorted[::-1]

        

    return df_sorted


def common_go_paralogs(data_go,data_paralogs):

 """
 a function that finds the common go terms for paralogs genes
 paralogs: paralogs are homologous genes that have evolved by duplication and code for protein with similar, but not identical functions.
 input: data_go= dataframe where all to go terms are 
        data_paralogs= dataframe where all the paralogs are
 output: dataframe with the fraction of common go for each paralog pair
 """
 
 query=np.unique(np.array(data_paralogs['query']))
 # big for loop for each gene analyzed in common partners
 d2=defaultdict(dict)

    
 for genes,i in zip(data_paralogs['paralogue-name'],query):
    d2[genes]['query']=i
    d2[genes]['names of paralogue']=genes

    tmp=data_go[data_go['Gene']==i]['go-term'].tolist()
    tmp=np.unique(tmp).tolist()

    

    tmp2=data_go[data_go['Gene']==genes]['go-term'].tolist()
    tmp2=np.unique(tmp2).tolist()

                


    d2[genes]['common-go-terms']=np.intersect1d(tmp,tmp2)
    if len(tmp)==0:
        d2[genes]['fraction-of-common-go']=0
    else:
        d2[genes]['fraction-of-common-go']=len(np.intersect1d(tmp,tmp2))/len(tmp) *100

    
    
 df=pd.DataFrame(d2).T
 df_sorted=df.sort_values(by='fraction-of-common-go',ascending=False)
      
    
    

 return df_sorted


def common_go(data_go,data_common_partners):
    """"
    function that computes the common go terms or interactors genes
    input: data_go= dataframe with all go terms per gene
    data_common_partners= dataframe output of the function common_partners
    
    """
 
    query=np.unique(np.array(data_common_partners['query']))
    # big for loop for each gene analyzed in common partners
    for i in np.arange(0,len(query)):
        partners=data_common_partners[data_common_partners['query']==query[i]]['names of genes']

        d2=defaultdict(dict)

        for genes in partners:
            d2[genes]['query']=query[i]
            d2[genes]['names of genes']=genes

            tmp=data_go[data_go['Gene']==query[i]]['go-term'].tolist()
            tmp=np.unique(tmp).tolist()

            

            tmp2=data_go[data_go['Gene']==genes]['go-term'].tolist()
            tmp2=np.unique(tmp2).tolist()

                        

       
            d2[genes]['common-go-terms']=np.intersect1d(tmp,tmp2)
            if len(tmp)==0:
                d2[genes]['fraction-of-common-go']=0
            else:
                d2[genes]['fraction-of-common-go']=len(np.intersect1d(tmp,tmp2))/len(tmp) *100

        

    df=pd.DataFrame(d2).T
    df_sorted=df.sort_values(by='fraction-of-common-go',ascending=False)
  



    return df_sorted
