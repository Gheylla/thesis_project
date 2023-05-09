# -*- coding: utf-8 -*-
"""
Created on Thu May  4 13:55:01 2023

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


# %%
'''Import the data'''
metrics_mart = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\Observations\MetricsRef1_MARTINpaper.h5')
metrics_hauk = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\Observations\GOES16_IR_nc_Iorg_EUREC4A_10-20_-58--48.nc')
metrics_gheyl = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\complete_observations.h5')


#%%
'''Compute the diurnal cycle of the metrics'''
diur_gheyl_iorg = metrics_gheyl.groupby([metrics_gheyl.index.hour])['iorg'].mean()
diur_hauk_iorg = metrics_hauk.Iorg.groupby(metrics_hauk.time.dt.hour).mean()

diur_gheyl_num = metrics_gheyl.groupby([metrics_gheyl.index.hour])['num_objects'].mean()
diur_hauk_num = metrics_hauk.N.groupby(metrics_hauk.time.dt.hour).mean()

diur_gheyl_cf = metrics_gheyl.groupby([metrics_gheyl.index.hour])['cloud_fraction'].mean()
diur_hauk_cf = metrics_hauk.cloud_fraction.groupby(metrics_hauk.time.dt.hour).mean()

diur_hauk_size = metrics_hauk.cluster_size_mean.groupby(metrics_hauk.time.dt.hour).mean()

# %%
'''Plot graphs of time series'''
plt.figure(figsize = (25,10))
plt.plot(metrics_hauk.time.values, metrics_hauk.Iorg.values, label = 'Hauke CONT3XT data')
plt.plot(metrics_gheyl.Time.values[1:8879], metrics_gheyl.iorg.values[1:8879], label = 'Gheylla GOES-16 data')
plt.title('Timeseries of Iorg')
plt.xlabel('Time')
plt.ylabel('Iorg [-]')
plt.grid()
plt.legend()


plt.figure(figsize = (25,10))
plt.plot(metrics_hauk.time.values, metrics_hauk.N.values, label = 'Hauke CONT3XT data')
plt.plot(metrics_gheyl.Time.values[1:8879], metrics_gheyl.num_objects.values[1:8879], label = 'Gheylla GOES-16 data')
plt.title('Timeseries of number of objects')
plt.xlabel('Time')
plt.ylabel('Num Objects [-]')
plt.grid()
plt.legend()

plt.figure(figsize = (25,10))
plt.plot(metrics_hauk.time.values, metrics_hauk.cloud_fraction.values, label = 'Hauke CONT3XT data')
plt.plot(metrics_gheyl.Time.values[1:8879], metrics_gheyl.cloud_fraction.values[1:8879], label = 'Gheylla GOES-16 data')
plt.title('Timeseries of cloud fraction')
plt.xlabel('Time')
plt.ylabel('CF[-]')
plt.grid()
plt.legend()


# %%
'''Plot graphs of diurnal cycle'''
plt.figure(figsize = (25,10))
plt.plot(diur_hauk_iorg.hour.values, diur_hauk_iorg, label = 'Hauke CONT3XT data')
plt.plot(diur_gheyl_iorg.index.values, diur_gheyl_iorg, label = 'Gheylla GOES-16 data')
plt.title('Diurnal cycle of Iorg')
plt.xlabel('Time')
plt.ylabel('Iorg [-]')
plt.grid()
plt.legend()

plt.figure(figsize = (25,10))
plt.plot(diur_hauk_iorg.hour.values, diur_hauk_num, label = 'Hauke CONT3XT data')
plt.plot(diur_gheyl_iorg.index.values, diur_gheyl_num, label = 'Gheylla GOES-16 data')
plt.title('Diurnal cycle of number of objects')
plt.xlabel('Time')
plt.ylabel('Number Objects [-]')
plt.grid()
plt.legend()

plt.figure(figsize = (25,10))
plt.plot(diur_hauk_iorg.hour.values, diur_hauk_cf, label = 'Hauke CONT3XT data')
plt.plot(diur_gheyl_iorg.index.values, diur_gheyl_cf, label = 'Gheylla GOES-16 data')
plt.title('Diurnal cycle of cloud fraction')
plt.xlabel('Time')
plt.ylabel('CF [-]')
plt.grid()
plt.legend()



# %%
'''Compute the mean cluster size'''
S_gheyl = metrics_gheyl['fractal_dimension'][1:8879] / metrics_gheyl['num_objects'][1:8879]

plt.figure(figsize = (25,10))
plt.plot(metrics_hauk.time.values, metrics_hauk.cluster_size_mean.values, label = 'Hauke CONT3XT data')
plt.plot(metrics_gheyl.Time.values[1:8879], S_gheyl, label = 'Gheylla GOES-16 data')
plt.title('Timeseries of mean cluster size')
plt.xlabel('Time')
plt.ylabel('Mean size[-]')
plt.grid()
plt.legend()

plt.figure(figsize = (25,10))
plt.plot(diur_hauk_iorg.hour.values, diur_hauk_size, label = 'Hauke CONT3XT data')
plt.title('Diurnal cycle of mean size')
plt.xlabel('Time')
plt.ylabel('Cluster size [-]')
plt.grid()
plt.legend()