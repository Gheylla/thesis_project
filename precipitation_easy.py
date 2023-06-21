# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 08:43:17 2023

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

xr_precip_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\pr_data\pr_nohgtqs_hour.nc')
xr_precip_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\pr_data\pr_nohgtqsnoshal_hour.nc')
xr_precip_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\pr_data\pr_nohgtqsnouvmix_hour.nc')


# %%
'''Composite diurnal precipitation'''
precip_noHGTQS = xr_precip_noHGTQS.groupby(xr_precip_noHGTQS.time.dt.hour).mean()
precip_noHGTQS_noSHAL = xr_precip_noHGTQS_noSHAL.groupby(xr_precip_noHGTQS_noSHAL.time.dt.hour).mean()
precip_noHGTQS_noUVmix = xr_precip_noHGTQS_noUVmix.groupby(xr_precip_noHGTQS_noUVmix.time.dt.hour).mean()


#%%
'''Plot the data'''
hours = np.linspace(0,23,24)
lt = [20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

fig, ax = plt.subplots(layout='constrained', figsize = (15,5))
ax.plot(precip_noHGTQS.hour, (precip_noHGTQS.precip.values * (1/997) * 1000), label = 'HARMONIE noHGTQS', color = 'red')
ax.plot(precip_noHGTQS_noSHAL.hour, (precip_noHGTQS_noSHAL.precip.values * (1/997) * 1000), label = 'HARMONIE noHGTQS noSHAL', color = 'blue')
ax.plot(precip_noHGTQS_noUVmix.hour, (precip_noHGTQS_noUVmix.precip.values * (1/997) * 1000), label = 'HARMONIE noHGTQS noUVmix', color = 'green')



plt.title('Composite Diurnal Cycle of Precipitation')
ax.grid()




plt.ylabel('Precipitation [mm/hr]')
plt.legend()
ax.set_xticks(hours)
ax.set_xlabel('UTC Time [hr]')
secax = ax.secondary_xaxis('top')
secax.set_xticks(hours, lt)
secax.set_xlabel('Local Time [hr]')



# =============================================================================
# plt.figure(figsize = (15,10))
# plt.plot(data_pr_noHGTQS.time.values, ((mean_pr_noHGTQS * (1/997)) * 1000), label = 'HARMONIE noHGTQS')
# plt.plot(data_pr_noHGTQS_noSHAL.time.values, ((mean_pr_noHGTQS_noSHAL * (1/997)) * 1000), label = 'HARMONIE noHGTQS noSHAL')
# plt.plot(data_pr_noHGTQS_noUVmix.time.values, ((mean_pr_noHGTQS_noUVmix * (1/997)) * 1000), label = 'HARMONIE noHGTQS noUVmix')
# plt.title('Accumulated precipitation')
# plt.xlabel('Time')
# plt.ylabel('Accumulated Precipitation [mm]')
# plt.grid()
# plt.legend()
# 
# 
# plt.figure(figsize = (15,10))
# plt.plot(data_pr_noHGTQS.time.values[0:1440], hr_pr_noHGTQS, label = 'HARMONIE noHGTQS')
# plt.plot(data_pr_noHGTQS_noSHAL.time.values[0:1440], hr_pr_noHGTQS_noSHAL, label = 'HARMONIE noHGTQS noSHAL')
# plt.plot(data_pr_noHGTQS_noUVmix.time.values[0:1440], hr_pr_noHGTQS_noUVmix, label = 'HARMONIE noHGTQS noUVmix')
# plt.title('Precipitation')
# plt.xlabel('Time')
# plt.ylabel('Precipitation [mm/hr]')
# plt.grid()
# plt.legend()
# =============================================================================
