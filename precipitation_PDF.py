# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 11:03:35 2023

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
import os 
import datetime as dt
from cut_subgrid_module import cut_sub_grid
#%%

'''Import data'''

data_pr_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\pr_data\pr_noHGTQS_cut.nc')
data_pr_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\pr_data\pr_noHGTQS_noSHAL_cut.nc')
data_pr_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\pr_data\pr_noHGTQS_noUVmix_cut.nc')

#%%
'''Make mean of the area'''
mean_pr_noHGTQS = data_pr_noHGTQS.pr.mean(dim = ('x', 'y'))
mean_pr_noHGTQS_noSHAL = data_pr_noHGTQS_noSHAL.pr.mean(dim = ('x', 'y'))
mean_pr_noHGTQS_noUVmix = data_pr_noHGTQS_noUVmix.pr.mean(dim = ('x', 'y'))

#%%

#plt.plot(mean_pr_noHGTQS.values)

# =============================================================================
# min_pr_noHGTQS = np.min(mean_pr_noHGTQS.values)
# min_pr_noHGTQS_noSHAL = np.min(mean_pr_noHGTQS_noSHAL.values)
# min_pr_noHGTQS_noUVmix = np.min(mean_pr_noHGTQS_noUVmix.values)
# 
# loc_min_noHGTQS = np.where(mean_pr_noHGTQS.values == min_pr_noHGTQS )
# loc_min_noHGTQS_noSHAL = np.where(mean_pr_noHGTQS_noSHAL.values == min_pr_noHGTQS_noSHAL )
# loc_min_noHGTQS_noSHAL = np.where(mean_pr_noHGTQS_noUVmix.values == min_pr_noHGTQS_noUVmix )
# =============================================================================

mean_pr_noHGTQS = mean_pr_noHGTQS.values
mean_pr_noHGTQS_noSHAL = mean_pr_noHGTQS_noSHAL.values
mean_pr_noHGTQS_noUVmix = mean_pr_noHGTQS_noUVmix.values

#%%
'''determine the hourly rain fall'''

hr_pr_noHGTQS = np.zeros(len(data_pr_noHGTQS.time.values))
hr_pr_noHGTQS_noSHAL = np.zeros(len(data_pr_noHGTQS_noSHAL.time.values))
hr_pr_noHGTQS_noUVmix = np.zeros(len(data_pr_noHGTQS_noUVmix.time.values))


for i in range((len(data_pr_noHGTQS.time.values)) - 1):
    print(i,data_pr_noHGTQS.time.values[i] )
    if i == 0: 
        hr_pr_noHGTQS[i] = mean_pr_noHGTQS[i]
        hr_pr_noHGTQS_noSHAL[i] = mean_pr_noHGTQS_noSHAL[i]
        hr_pr_noHGTQS_noUVmix[i] = mean_pr_noHGTQS_noUVmix[i]
    else:
        pr_noHGTQS_before = mean_pr_noHGTQS[i - 1]
        pr_noHGTQS_noSHAL_before = mean_pr_noHGTQS_noSHAL[i - 1]
        pr_noHGTQS_noUVmix_before = mean_pr_noHGTQS_noUVmix[i - 1]
        pr_noHGTQS = mean_pr_noHGTQS[i]
        pr_noHGTQS_noSHAL = mean_pr_noHGTQS_noSHAL[i]
        pr_noHGTQS_noUVmix = mean_pr_noHGTQS_noUVmix[i]
        if i == 744:  #location where the new data starts
            #print(data_pr_noHGTQS.time.values[744])
            hr_pr_noHGTQS[i] = 0
            hr_pr_noHGTQS_noSHAL[i] = 0
            hr_pr_noHGTQS_noUVmix[i] = 0
        else: 
            hr_pr_noHGTQS[i] = pr_noHGTQS - pr_noHGTQS_before
            hr_pr_noHGTQS_noSHAL[i] = pr_noHGTQS_noSHAL - pr_noHGTQS_noSHAL_before
            hr_pr_noHGTQS_noUVmix[i] = pr_noHGTQS_noUVmix - pr_noHGTQS_noUVmix_before


#%% 
'''Make new dataset'''

xr_precip_noHGTQS = xr.Dataset({'precip': (('time'), hr_pr_noHGTQS )},
                  {'time': ( data_pr_noHGTQS.time.values)}, 
                  {'units': ('mm/hr'), 
                   'long_name': ('Precipitation mm/hr noHGTQS')})


xr_precip_noHGTQS_noSHAL = xr.Dataset({'precip': (('time'), hr_pr_noHGTQS_noSHAL )},
                  {'time': ( data_pr_noHGTQS.time.values)}, 
                  {'units': ('mm/hr'), 
                   'long_name': ('Precipitation mm/hr noHGTQS noSHAL')})


xr_precip_noHGTQS_noUVmix = xr.Dataset({'precip': (('time'), hr_pr_noHGTQS_noUVmix )},
                  {'time': ( data_pr_noHGTQS.time.values)}, 
                  {'units': ('mm/hr'), 
                   'long_name': ('Precipitation mm/hr noHGTQS noUVmix')})


#%%
xr_precip_noHGTQS.to_netcdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\pr_data\pr_nohgtqs_hour.nc')
xr_precip_noHGTQS_noSHAL.to_netcdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\pr_data\pr_nohgtqsnoshal_hour.nc')
xr_precip_noHGTQS_noUVmix.to_netcdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\pr_data\pr_nohgtqsnouvmix_hour.nc')

# %%
'''Composite diurnal precipitation'''
precip_noHGTQS = xr_precip_noHGTQS.groupby(xr_precip_noHGTQS.time.dt.hour).mean()
precip_noHGTQS_noSHAL = xr_precip_noHGTQS_noSHAL.groupby(xr_precip_noHGTQS_noSHAL.time.dt.hour).mean()
precip_noHGTQS_noUVmix = xr_precip_noHGTQS_noUVmix.groupby(xr_precip_noHGTQS_noUVmix.time.dt.hour).mean()



#%%
'''Plot the data'''
hours = np.linspace(0,23,24)


plt.figure(figsize = (15,10))
plt.plot(precip_noHGTQS.hour, (precip_noHGTQS.precip.values * (1/997) * 1000), label = 'HARMONIE noHGTQS')
plt.plot(precip_noHGTQS_noSHAL.hour, (precip_noHGTQS_noSHAL.precip.values * (1/997) * 1000), label = 'HARMONIE noHGTQS noSHAL')
plt.plot(precip_noHGTQS_noUVmix.hour, (precip_noHGTQS_noUVmix.precip.values * (1/997) * 1000), label = 'HARMONIE noHGTQS noUVmix')
plt.title('Composite Diurnal Cycle of Precipitation')
plt.grid()
plt.xlabel('Time [hr]')
plt.ylabel('Precipitation [mm]')
plt.legend()
plt.xticks(precip_noHGTQS_noUVmix.hour, hours)


plt.figure(figsize = (15,10))
plt.plot(data_pr_noHGTQS.time.values, ((mean_pr_noHGTQS * (1/997)) * 1000), label = 'HARMONIE noHGTQS')
plt.plot(data_pr_noHGTQS_noSHAL.time.values, ((mean_pr_noHGTQS_noSHAL * (1/997)) * 1000), label = 'HARMONIE noHGTQS noSHAL')
plt.plot(data_pr_noHGTQS_noUVmix.time.values, ((mean_pr_noHGTQS_noUVmix * (1/997)) * 1000), label = 'HARMONIE noHGTQS noUVmix')
plt.title('Accumulated precipitation')
plt.xlabel('Time')
plt.ylabel('Accumulated Precipitation [mm]')
plt.grid()
plt.legend()


plt.figure(figsize = (15,10))
plt.plot(data_pr_noHGTQS.time.values[0:1440], hr_pr_noHGTQS, label = 'HARMONIE noHGTQS')
plt.plot(data_pr_noHGTQS_noSHAL.time.values[0:1440], hr_pr_noHGTQS_noSHAL, label = 'HARMONIE noHGTQS noSHAL')
plt.plot(data_pr_noHGTQS_noUVmix.time.values[0:1440], hr_pr_noHGTQS_noUVmix, label = 'HARMONIE noHGTQS noUVmix')
plt.title('Precipitation')
plt.xlabel('Time')
plt.ylabel('Precipitation [mm/hr]')
plt.grid()
plt.legend()
