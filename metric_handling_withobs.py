# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 09:33:01 2023

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
import datetime as dt
import h5py
import os


observations = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\complete_observations.h5')
obs_new = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_293_286_test.h5')
df_metrics_noHGTQS43 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\HARMONIE\df_metrics_noHGTQS.h5')
df_metrics_noHGTQS_noSHAL = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\HARMONIE\df_metrics_noHGTQS_noSHAL.h5')
df_metrics_noHGTQS_noUVmix = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\HARMONIE\df_metrics_noHGTQS_noUVmix.h5')


metrics = observations.keys()

for i in range(len(metrics) -1):
    plt.figure(figsize = (25,10))
    #plt.plot(df_metrics_noHGTQS43.index, df_metrics_noHGTQS43[metrics[i]], label = 'HARMONIE noHGTQS', color = 'red')
    #plt.plot(df_metrics_noHGTQS_noSHAL.index, df_metrics_noHGTQS_noSHAL[metrics[i]], label = 'HARMONIE noHGTQS noSHAL', color = 'green')
    #plt.plot(df_metrics_noHGTQS_noUVmix.index, df_metrics_noHGTQS_noUVmix[metrics[i]], label = 'HARMONIE noHGTQS noUVmix', color = 'blue')
    plt.plot(observations.Time[3938:3965], observations[metrics[i]][3938:3965], label = '290-280K', color = 'black', linewidth = 3.5)
    plt.plot(obs_new.Time[1:], obs_new[metrics[i]][1:], label = '293-286K', color = 'red', linewidth = 3.5)
    plt.legend(loc='upper right')
    plt.grid()
    plt.title(f' Timeseries of {metrics[i]}')
    plt.xlabel('Time')
    plt.ylabel(f'{metrics[i]} [-]')
    
# =============================================================================
# hours = np.linspace(0,23,24)
# for i in range(len(metrics) -1):
#     noHGTQS_diur = df_metrics_noHGTQS43[metrics[i]].groupby(df_metrics_noHGTQS43.index.hour).mean()
#     noHGTQS_noSHAL_diur = df_metrics_noHGTQS_noSHAL[metrics[i]].groupby(df_metrics_noHGTQS_noSHAL.index.hour).mean()
#     noHGTQS_noUVmix_diur = df_metrics_noHGTQS_noUVmix[metrics[i]].groupby(df_metrics_noHGTQS_noSHAL.index.hour).mean()    
#     obs_diur = observations.groupby([observations.index.hour])[metrics[i]].mean()
#     plt.figure(figsize =(15,10))
#     plt.plot(noHGTQS_diur.index, noHGTQS_diur.values, color = 'red', label = 'HARMONIE noHGTQS')
#     plt.plot(noHGTQS_noSHAL_diur.index, noHGTQS_noSHAL_diur.values, color = 'green', label = 'HARMONIE noHGTQS noSHAL')
#     plt.plot(noHGTQS_noUVmix_diur.index, noHGTQS_noUVmix_diur.values, color = 'blue', label = 'HARMONIE noHGTQS noUVmix')
#     plt.plot(obs_diur.index, obs_diur.values, label = 'Observations', color = 'black', linewidth=3.5)
#     plt.legend(loc='upper right')
#     plt.grid()
#     plt.title(f'Composite diurnal cycle of {metrics[i]}')
#     plt.xlabel('Time [hr]')
#     plt.ylabel(f'{metrics[i]} [-]')
#     plt.xticks(hours)
# =============================================================================
    