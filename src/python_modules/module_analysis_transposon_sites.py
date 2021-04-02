# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 15:17:00 2021

@author: linigodelacruz
"""

import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from collections import defaultdict
import os
import seaborn as sns
import scipy 
#%% Analysis of transposition sites

def frequency_transposons(data,names_libraries):
    freq=[]
    for i in names_libraries.keys():
        freq.append(data.loc[i]['Nbasepairs'].sum()/data.loc[i]['Ninsertions'].sum())
        
    return freq

def reads_per_transposon(data,names_libraries):
    readspertr=[]
    for i in names_libraries.keys():
        readspertr.append(data.loc[i]['Nreads'].median()/data.loc[i]['Ninsertions'].median())
        
    return readspertr
    
def transposon_density(data,names_libraries):
    density=[]
    for i in names_libraries.keys():
        
        density.append(data.loc[i]['Ninsertions']/data.loc[i]['Nbasepairs'])        
    return density
def median_feature(data,names_libraries,feature):
    median_feature=[]
    for i in names_libraries.keys():
        
        median_feature.append(data.loc[i][feature].median())        
    return median_feature
    
def median_feature_essentials(data,names_libraries,feature):
    median_feature=[]
    for i in names_libraries.keys():
        
        median_feature.append(data.loc[i][data.loc[i]['Essentiality']==1][feature].median())        
    return median_feature
    
def median_feature_nonessentials(data,names_libraries,feature):
    median_feature=[]
    for i in names_libraries.keys():
        
        median_feature.append(data.loc[i][data.loc[i]['Essentiality']==0][feature].median())        
    return median_feature