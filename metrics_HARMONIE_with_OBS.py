# -*- coding: utf-8 -*-
"""
Created on Mon May 15 14:12:08 2023

@author: Gheylla Liberia
This code tries to analyze the various HARMONIE metric outputs
"""

#%%
'''Import packages'''
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


'''Import data'''
harm_nohgtqs = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\HARMONIE\df_metrics_noHGTQS.h5')
harm_nohgtqs_noshal = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\HARMONIE\df_metrics_noHGTQS_noSHAL.h5')
harm_nohgtqs_nouvmix = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\HARMONIE\df_metrics_noHGTQS_noUVmix.h5')


#import observation data
goes_metrics277 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\complete_observations_277290.h5')
goes_metrics280 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\complete_observations_280290.h5')
#goes_metrics277 = goes_metrics277.dropna(axis = 0)

# %%
'''Plot all data '''
plt.figure(figsize = (25,10))
plt.plot(harm_nohgtqs.index.values, harm_nohgtqs['cloud_fraction'], color = 'red', label = 'Metrics cll HARMONIE noHGTQS')
plt.plot(harm_nohgtqs_noshal.index.values, harm_nohgtqs_noshal['cloud_fraction'], color = 'blue', label = 'Metrics cll HARMONIE noHGTQS noSHAL')
plt.plot(harm_nohgtqs_nouvmix.index.values, harm_nohgtqs_nouvmix['cloud_fraction'], color = 'green', label = 'Metrics cll HARMONIE noHGTQS noUVmix')

plt.plot(goes_metrics277.index.values[1:8879], goes_metrics277['cloud_fraction'][1:8879], label = 'Observations mask 280-290K', color = 'black', linestyle = 'dashed')
plt.plot(goes_metrics280.index.values[1:8879], goes_metrics280['cloud_fraction'][1:8879], label = 'Observations mask 277-290K', color = 'black')


plt.grid()
plt.title('Timeseries of low cloud cover')
plt.xlabel('Time')
plt.ylabel('CC [-]')
plt.legend(loc = 'upper right')



# %%
'''Make the timeseries plots for all other metrics'''
keys = harm_nohgtqs.keys()

for i in range(len(keys) -1):
    plt.figure(figsize = (25,10))
    plt.plot(harm_nohgtqs.index.values, harm_nohgtqs[keys[i]].values, color = 'red', label = 'HARMONIE noHGTQS')
    plt.plot(harm_nohgtqs_noshal.index.values, harm_nohgtqs_noshal[keys[i]].values, color = 'blue', label = 'HARMONIE noHGTQS noSHAL')
    plt.plot(harm_nohgtqs_nouvmix.index.values, harm_nohgtqs_nouvmix[keys[i]].values, color = 'green', label = 'HARMONIE noHGTQS noUVmix')
    
    plt.grid()
    plt.title(f'Timeseries of {keys[i]}')
    plt.xlabel('Time')
    plt.ylabel(f'{keys[i]} [-]')
    plt.legend(loc = 'upper right')
    
    
# %%
'''Make diurnal cycles of all metrics'''
keys = harm_nohgtqs.keys()
hours = np.linspace(00,23, 24)
lt = [20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]


for ii in range(len(keys)-1):
    print(keys[ii])
    diur_goes_277 = goes_metrics277[keys[ii]].groupby(goes_metrics277.index.hour).mean()
    diur_goes_280 = goes_metrics280[keys[ii]].groupby(goes_metrics280.index.hour).mean()
    diur_nohgtqs = harm_nohgtqs[keys[ii]].groupby(harm_nohgtqs.index.hour).mean()
    diur_nohgtqs_noshal = harm_nohgtqs_noshal[keys[ii]].groupby(harm_nohgtqs_noshal.index.hour).mean()
    diur_nohgtqs_nouvmix = harm_nohgtqs_nouvmix[keys[ii]].groupby(harm_nohgtqs_nouvmix.index.hour).mean()
    
    fig, ax = plt.subplots(layout='constrained', figsize = (15,10))
    

    ax.plot(diur_nohgtqs.index.values, diur_nohgtqs, color = 'red', label = 'HARMONIE noHGTQS' )
    ax.plot(diur_nohgtqs_noshal.index.values, diur_nohgtqs_noshal, color = 'blue', label = 'HARMONIE noHGTQS noSHAL' )
    ax.plot(diur_nohgtqs_nouvmix.index.values, diur_nohgtqs_nouvmix, color = 'green', label = 'HARMONIE noHGTQS noUVmix' )
    ax.plot(diur_goes_277.index.values, diur_goes_277, color = 'black', label = 'Observations mask 277-290K')
    ax.plot(diur_goes_280.index.values, diur_goes_280, color = 'black', linestyle = 'dashed', label = 'Observations mask 280-290K')
    
    ax.grid()
    ax.set_xticks(hours)
    ax.set_xlabel('UTC Time [hr]')
    secax = ax.secondary_xaxis('top')
    secax.set_xticks(hours, lt)
    secax.set_xlabel('Local Time [hr]')
    plt.ylabel(f'{keys[ii]} [-]')
    plt.title(f'Composite Diurnal cycle of {keys[ii]}')
    plt.legend()
    plt.show()
    

# %%
'''Make correlation plot of all metrics'''
keys = harm_nohgtqs.keys()

for jj in range(len(keys) - 1):
    plt.figure()
    plt.scatter(goes_metrics277[keys[jj]][1:8879], goes_metrics277[keys[jj + 1 ]][1:8879], s = 1, color = 'black')
    plt.grid()
    plt.xlabel(f'{keys[jj]}')
    plt.ylabel(f'{keys[jj + 1]}')
    plt.title('Correlation for GOES 16 metrics')


plt.figure()
plt.scatter(goes_metrics277['mean_length_scale'][1:8879], goes_metrics277['open_sky'][1:8879], s = 1, color = 'black')
plt.grid()
plt.xlabel('mean length scale')
plt.ylabel('open sky')
plt.title('Correlation for GOES 16 metrics')

    