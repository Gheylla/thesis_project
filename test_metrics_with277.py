# -*- coding: utf-8 -*-
"""
Created on Tue May  9 15:13:48 2023

@author: LENOVO
"""

import numpy
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import scipy as sc
from scipy.signal import argrelmin 
import datetime as dt
import inspect 
import cloudmetrics
import plotly.graph_objects as go
import seaborn as sn
import os

#data = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_2000.h5')
observations = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\complete_observations_280290.h5')
data2 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_1195.h5')
data3 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_2000.h5')
data4 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_3950.h5')
data5 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_4000.h5')
data6 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_6000.h5')


plt.figure(figsize = (15,10))
plt.plot(data2['Time'],data2['cloud_fraction'])
plt.plot(data3['Time'],data3['cloud_fraction'])
#plt.plot(data4['Time'],data4['cloud_fraction'])
plt.plot(data5['Time'],data5['cloud_fraction'])
plt.plot(data6['Time'],data6['cloud_fraction'])
plt.plot(observations['Time'][1:8000] , observations['cloud_fraction'][1:8000], color = 'black')


# =============================================================================
# 
# keys = data2.keys()
# 
# for i in range(len(keys)-1):
#     print(keys[i])
#     
#     plt.figure(figsize= (25,10))
#     plt.plot(observations[keys[i]][1:8879])
#     #plt.plot(observations['Time'][1:1967], data[keys[i]])
#     plt.plot(data2['Time'], data2[keys[i]])
#     plt.plot(data3['Time'], data3[keys[i]])
#     plt.plot(data4['Time'].values, data4)
# 
# 
#     plt.title(keys[i])
# =============================================================================

