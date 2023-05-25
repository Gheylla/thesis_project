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
#metrics_mart = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\Observations\MetricsRef1_MARTINpaper.h5')
metrics_hauk = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\Observations\GOES16_IR_nc_Iorg_EUREC4A_10-20_-58--48.nc') #https://www.zenodo.org/record/5979718#.ZGHwp3ZBxD8
metrics_gheyl_280290 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\complete_observations_280290.h5')
metrics_gheyl_277290 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\complete_observations_277290.h5')


#%%
'''Compute the diurnal cycle of the metrics'''


diur_gheyl_iorg_280290 = metrics_gheyl_280290.groupby([metrics_gheyl_280290.index.hour])['iorg'].mean()
diur_gheyl_iorg_277290 = metrics_gheyl_277290.groupby([metrics_gheyl_277290.index.hour])['iorg'].mean()
diur_hauk_iorg = metrics_hauk.Iorg.groupby(metrics_hauk.time.dt.hour).mean()

diur_gheyl_num_280290 = metrics_gheyl_280290.groupby([metrics_gheyl_280290.index.hour])['num_objects'].mean()
diur_gheyl_num_277290 = metrics_gheyl_277290.groupby([metrics_gheyl_277290.index.hour])['num_objects'].mean()
diur_hauk_num = metrics_hauk.N.groupby(metrics_hauk.time.dt.hour).mean()

diur_gheyl_cf_280290 = metrics_gheyl_280290.groupby([metrics_gheyl_280290.index.hour])['cloud_fraction'].mean()
diur_gheyl_cf_277290 = metrics_gheyl_277290.groupby([metrics_gheyl_277290.index.hour])['cloud_fraction'].mean()
diur_hauk_cf = metrics_hauk.cloud_fraction.groupby(metrics_hauk.time.dt.hour).mean()

diur_hauk_size = metrics_hauk.cluster_size_mean.groupby(metrics_hauk.time.dt.hour).mean()

# %%
'''Plot graphs of time series'''
plt.figure(figsize = (25,10))
plt.plot(metrics_hauk.time.values, metrics_hauk.Iorg.values, label = 'Hauke CONT3XT data', color = 'red')



plt.plot(metrics_gheyl_280290.Time.values[1:8879], metrics_gheyl_280290.iorg.values[1:8879], label = 'Gheylla GOES-16 data mask 280-290K', color = 'black')
plt.plot(metrics_gheyl_277290.Time.values[1:8879], metrics_gheyl_277290.iorg.values[1:8879], label = 'Gheylla GOES-16 data mask 277-290K', color = 'black', linestyle = 'dashed')
plt.title('Timeseries of Iorg')
plt.xlabel('Time')
plt.ylabel('Iorg [-]')
plt.grid()
plt.legend()

#%%
plt.figure(figsize = (25,10))
plt.plot(metrics_hauk.time.values, metrics_hauk.N.values, label = 'Hauke CONT3XT data', color = 'red')

plt.plot(metrics_gheyl_280290.Time.values[1:8879], metrics_gheyl_280290.num_objects.values[1:8879], label = 'Gheylla GOES-16 data mask 280-290K', color = 'black')
plt.plot(metrics_gheyl_277290.Time.values[1:8879], metrics_gheyl_277290.num_objects.values[1:8879], label = 'Gheylla GOES-16 data mask 277-290K', color = 'black',  linestyle = 'dashed')
plt.title('Timeseries of number of objects')
plt.xlabel('Time')
plt.ylabel('Num Objects [-]')
plt.grid()
plt.legend()

#%%

plt.figure(figsize = (25,10))
plt.plot(metrics_hauk.time.values, metrics_hauk.cloud_fraction.values, label = 'Hauke CONT3XT data', color = 'red')


plt.plot(metrics_gheyl_280290.Time.values[1:8879], metrics_gheyl_280290.cloud_fraction.values[1:8879], label = 'Gheylla GOES-16 data mask 280-290K', color = 'black')
plt.plot(metrics_gheyl_277290.Time.values[1:8879], metrics_gheyl_277290.cloud_fraction.values[1:8879], label = 'Gheylla GOES-16 data mask 277-290K', color = 'black', linestyle = 'dashed')
plt.title('Timeseries of cloud fraction')
plt.xlabel('Time')
plt.ylabel('CF[-]')
plt.grid()
plt.legend()


# %%
'''Plot graphs of diurnal cycle'''
hours = np.linspace(00,23, 24)
lt = [20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

fig, ax = plt.subplots(layout='constrained', figsize = (15,10))


ax.plot(diur_hauk_iorg.hour.values, diur_hauk_iorg, label = 'Hauke CONT3XT data', color = 'red')
ax.plot(diur_gheyl_iorg_280290.index.values, diur_gheyl_iorg_280290, label = 'Gheylla GOES-16 data mask 280-290K', color = 'black')
ax.plot(diur_gheyl_iorg_277290.index.values, diur_gheyl_iorg_277290, label = 'Gheylla GOES-16 data mask 277-290K', color = 'black', linestyle = 'dashed')

ax.grid()
ax.set_xticks(hours)
ax.set_xlabel('UTC Time [hr]')
secax = ax.secondary_xaxis('top')
secax.set_xticks(hours, lt)
secax.set_xlabel('Local Time [hr]')

plt.title('Diurnal cycle of Iorg')
plt.ylabel('Iorg [-]')
plt.legend()

#%%
hours = np.linspace(00,23, 24)
lt = [20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

fig, ax = plt.subplots(layout='constrained', figsize = (25,10))


ax.plot(diur_hauk_iorg.hour.values, diur_hauk_num, label = 'Hauke CONT3XT data', color = 'red')
ax.plot(diur_gheyl_iorg_280290.index.values, diur_gheyl_num_280290, label = 'Gheylla GOES-16 data mask 280-290K', color = 'black')
ax.plot(diur_gheyl_iorg_277290.index.values, diur_gheyl_num_277290, label = 'Gheylla GOES-16 data mask 277-290K', color = 'black', linestyle = 'dashed')

ax.grid()
ax.set_xticks(hours)
ax.set_xlabel('UTC Time [hr]')
secax = ax.secondary_xaxis('top')
secax.set_xticks(hours, lt)
secax.set_xlabel('Local Time [hr]')

plt.title('Diurnal cycle of number of objects')
plt.ylabel('Number Objects [-]')
plt.legend()

#%%
hours = np.linspace(00,23, 24)
lt = [20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

fig, ax = plt.subplots(layout='constrained', figsize = (15,7))

ax.plot(diur_hauk_iorg.hour.values, diur_hauk_cf, label = 'Hauke CONT3XT data', color = 'red')
ax.plot(diur_gheyl_iorg_280290.index.values, diur_gheyl_cf_280290, label = 'Gheylla GOES-16 data mask 280-290K', color = 'black')
ax.plot(diur_gheyl_iorg_277290.index.values, diur_gheyl_cf_277290, label = 'Gheylla GOES-16 data mask 277-290K', color = 'black', linestyle = 'dashed')

ax.grid()
ax.set_xticks(hours)
ax.set_xlabel('UTC Time [hr]')
secax = ax.secondary_xaxis('top')
secax.set_xticks(hours, lt)
secax.set_xlabel('Local Time [hr]')

plt.title('Diurnal cycle of cloud fraction')
plt.ylabel('CF [-]')
plt.legend()



