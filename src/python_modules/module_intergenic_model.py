# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 11:31:37 2021

@author: linigodelacruz
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 

def getting_r(datasets,time=False,count_for_zero_tr=False): 
    """
    This function computes the maximum rates per gene,effectively the fitness for 
    a strain with a knockout in gene X, given an intergenic model . The intergenic model assumes a 
    fitine population growth until the carrying capacity (K). This population growth describes a 
    logistic growth. 
    
    Parameters
    ----------
    datasets : TYPE - list of dataframes
        These are dataframes containing the information on insertions per gene and reads per gene 
        given by the transposon data sequencing.
        The columns should have the following names:
            number_of_read_per_gene
            number_of_transposon_per_gene
            

    Returns
    -------
    r : array of length of the list given as input 
        array containing the maximum rates per gene per dataset in the population according an intergenic model
    reads_per_transposons : 
        array of length of the list given as input 
        array containing the reads per transposons per dataset. 

    """
    ## add a for loop over times and compute the rates over those times and averaged them out and std 
    total_reads=[]
    T=90
    interval=np.linspace(1,90,10)
    r=[]
    reads_per_transposons=[]
    r_t=np.zeros(shape=(len(datasets),len(datasets[0]),len(interval)))
    
    
    
    for i in np.arange(0,len(datasets)):
        tr_exceptions=np.where((datasets[i]['number_of_transposon_per_gene']==1) | (datasets[i]['number_of_transposon_per_gene']==0) )[0]
        reads_exceptions=np.where((datasets[i]['number_of_read_per_gene']==0) )[0]
        exceptions=np.intersect1d(tr_exceptions,reads_exceptions)
        total_reads=np.sum(datasets[i]['number_of_read_per_gene'])
        total_tn=np.sum(datasets[i]['number_of_transposon_per_gene'])
        K=total_reads/total_tn # it will compute a carrying capacity per dataset 
        N=datasets[i]['number_of_read_per_gene']/(datasets[i]['number_of_transposon_per_gene']-1)#there is one more transposon from the maximum reads
        
        reads_per_transposons.append(N)
        
        r.append(np.log(N/(1-N/K))/T)
        
        if count_for_zero_tr==True: ## not working yet think ina diiferent way to do it 
            for j in datasets[i].index:
                                            
                if j not in exceptions: 
                    reads_per_transposons.append(datasets[i]['number_of_read_per_gene'][j]/(datasets[i]['number_of_transposon_per_gene'][j]-1))#there is one more transposon from the maximum reads
                
                
                else:
                    reads_per_transposons.append(25/(5-1))
    
            
       
        
            
        if time==True: # NOT WORKING YET 
            for t in np.arange(0,len(interval)):
                 
                 r_t.append(np.log(reads_per_transposons/(1-reads_per_transposons/K))/interval[t])
                 
                 
                
        
       
    
    return r,reads_per_transposons

### Configuring the dataframes for analyses
# def pd_to_compare_libraries(data_list,datasets_for_index,norm=True):
#     pd_list=pd.DataFrame(data_list,index=np.arange(0,len(datasets_for_index)))
#     pd_list=pd_list.T
       
    
#     pd_list.columns=['nrp1_1_a','nrp1_1_b','nrp1_2_a','nrp2_2_b','wt_a','wt_b']
    
#     pd_list.index=datasets_for_index[0]['gene_name']
#     pd_list.fillna(0,inplace=True)
#     pd_list.replace(float('-inf'),0,inplace=True)
    
#     pd_list['average_dnrp1']=np.mean(pd_list.loc[:,['nrp1_1_a','nrp1_1_b','nrp1_2_a','nrp2_2_b']],axis=1)
#     pd_list['std_dnrp1']=np.std(pd_list.loc[:,['nrp1_1_a','nrp1_1_b','nrp1_2_a','nrp2_2_b']],axis=1)
    
#     pd_list['merged_dnrp1']=np.sum(pd_list.loc[:,['nrp1_1_a','nrp1_1_b','nrp1_2_a','nrp2_2_b']],axis=1)
    
    
#     pd_list['average_wt']=np.mean(pd_list.loc[:,['wt_a','wt_b']],axis=1)
#     pd_list['std_wt']=np.std(pd_list.loc[:,['wt_a','wt_b']],axis=1)
    
#     pd_list['merged_wt']=np.sum(pd_list.loc[:,['wt_a','wt_b']],axis=1)
    
    
#     pd_list['average_ratio_mutant_vs_wt']=pd_list['average_dnrp1']/pd_list['average_wt']
#     if norm==True:
        
#         pd_list['norm2max_wt']=pd_list.loc[:,'average_wt']/pd_list.loc[:,'average_wt'].max()
#         pd_list['norm2max_dnrp1']=pd_list.loc[:,'average_dnrp1']/pd_list.loc[:,'average_dnrp1'].max()
        

        
#     return pd_list


def pd_to_compare_libraries(data_list,datasets_for_index,filesnames,norm=True):
    pd_list=pd.DataFrame(data_list,index=np.arange(0,len(datasets_for_index)))
    pd_list=pd_list.T
       
    columns=[]
    
    for i in np.arange(0,len(filesnames)):
        if filesnames[i].startswith('d'): #dnrp1
            columns.append(filesnames[i][0:7]+filesnames[i][27:29])
        elif filesnames[i].startswith('W'):
            columns.append(filesnames[i][0:2]+filesnames[i][22:24])
            
        
    pd_list.columns=columns
    
    pd_list.index=datasets_for_index[0]['gene_name']
    pd_list.fillna(0,inplace=True)
    pd_list.replace(float('-inf'),0,inplace=True)
    
    pd_dnrp1=pd_list.filter(regex='nrp')
    
    pd_list['average_dnrp1']=np.mean(pd_dnrp1,axis=1)
    pd_list['std_dnrp1']=np.std(pd_dnrp1,axis=1)
    
    pd_list['merged_dnrp1']=np.sum(pd_dnrp1,axis=1)
    
    pd_wt=pd_list.filter(regex='WT')
    pd_list['average_wt']=np.mean(pd_wt,axis=1)
    pd_list['std_wt']=np.std(pd_wt,axis=1)
    
    pd_list['merged_wt']=np.sum(pd_wt,axis=1)
    
    
    pd_list['average_ratio_mutant_vs_wt']=pd_list['average_dnrp1']/pd_list['average_wt']
    if norm==True:
        
        pd_list['norm2max_wt']=pd_list.loc[:,'average_wt']/pd_list.loc[:,'average_wt'].max()
        pd_list['norm2max_dnrp1']=pd_list.loc[:,'average_dnrp1']/pd_list.loc[:,'average_dnrp1'].max()
        

        
    return pd_list
### function to filter the data 
def filtering_significant_genes(dataset,p_value_th,fc_th):
    dataset_highfc=dataset[dataset.loc[:,'fc_log2']>fc_th] 
    dataset_highp=dataset[dataset.loc[:,'-log_p']>p_value_th]
    dataset_significant_genes_wt= pd.merge(dataset_highfc, dataset_highp, how='inner', on=['-log_p'],left_index=True, right_index=True)
    dataset_significant_genes_wt.drop(columns=['fc_log2_y'],inplace=True)
    dataset_significant_genes_wt.rename(columns={'fc_log2_x': 'fc_log2'},inplace=True)
    
    
    dataset_highfc=dataset[dataset.loc[:,'fc_log2']<-fc_th] 
    dataset_highp=dataset[dataset.loc[:,'-log_p']>p_value_th] 
    dataset_significant_genes_mutant= pd.merge(dataset_highfc, dataset_highp, how='inner', on=['-log_p'],left_index=True, right_index=True)
    dataset_significant_genes_mutant.drop(columns=['fc_log2_y'],inplace=True)
    dataset_significant_genes_mutant.rename(columns={'fc_log2_x': 'fc_log2'},inplace=True)
    
    return dataset_significant_genes_wt,dataset_significant_genes_mutant


def viz_data(dataset,p_value_th,fc_th,filenames,saving=True,plotting=True):
    
    results=[]
    if type(dataset)==list:
        
        for i in np.arange(0,len(dataset)):
            
            ref,exp=filtering_significant_genes(dataset[i],p_value_th,fc_th)
            results.append([ref,exp])
    else : 
        
        ref,exp=filtering_significant_genes(dataset,p_value_th,fc_th)
        results.append([ref,exp])        
### Plotting all datasets 
    for data in np.arange(0,len(results)):
        fig = plt.figure(figsize=(10,5))
        ax = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        
        data_ref=results[data][0]
        data_exp=results[data][1]
        
        for i in range(len(data_ref)):
            x = data_ref.loc[data_ref.index[i],'fc_log2']
            y = data_ref.loc[data_ref.index[i],'-log_p']
            ax.plot(x, y, 'bo')
        
            ax.text(x*(1-0.02) , y*(1+0.02)  , data_ref.index[i], fontsize=4)
         
        for i in range(len(data_exp)):
            x = data_exp.loc[data_exp.index[i],'fc_log2']
            y = data_exp.loc[data_exp.index[i],'-log_p']
            ax2.plot(x, y, 'bo')
           
            ax2.text(x*(1-0.02) , y*(1+0.02)  , data_exp.index[i], fontsize=4)
        
        for axes in [ax,ax2]:
            axes.set_ylim((0, 10))
            axes.set_xlabel('fc_log2')
            axes.set_ylabel('-log_p')
            if axes==ax:
                axes.set_xlim((1, 5*fc_th))
                axes.set_title('Signif for WT-'+ filenames,fontsize=8)
            else:
                axes.set_xlim((-5*fc_th, -1))
                axes.set_title('Signif for mutant-'+ filenames,fontsize=8)
                
        if plotting==True:
            plt.show()
        if saving==True:
            fig.savefig(filenames + '.png',dpi=300,format='png',transparent=False)
