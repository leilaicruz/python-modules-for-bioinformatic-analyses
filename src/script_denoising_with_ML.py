# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 15:23:33 2021

@author: linigodelacruz
"""
from sklearn.neighbors import KNeighborsRegressor
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd

df=data_wt
clf = KNeighborsRegressor(n_neighbors=100, weights='uniform')
y=df.loc[:, 'Nreads']
clf.fit(df.index.values[:, np.newaxis], y)
y_pred = clf.predict(df.index.values[:, np.newaxis])
ax = pd.Series(y).plot(color='lightgray')
pd.Series(y_pred).plot(color='black', ax=ax, figsize=(12, 8))

#ax.set_ylim(0,40000)
