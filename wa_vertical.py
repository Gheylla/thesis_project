# -*- coding: utf-8 -*-
"""
Created on Wed May 17 13:43:00 2023

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

# %%
'''Import data'''
#remember that because this is dependend on the level the area is different. The area is smaller than what we are actually analysing. 

w_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\w_data\wa_noHTQS_total.nc')
w_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\w_data\wa_noHTQS_noSHAL_total.nc')
w_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\w_data\wa_noHTQS_noUVmix_total.nc')

# %%
'''Determine the mean of the area'''
w_noHGTQS_mean = w_noHGTQS.mean(dim = ('x', 'y'))
w_noHGTQS_noSHAL_mean = w_noHGTQS_noSHAL.mean(dim = ('x', 'y'))
w_noHGTQS_noUVmix_mean = w_noHGTQS_noUVmix.mean(dim = ('x', 'y'))

# %%
'''Groupby hours to get diurnal cycle'''
w_noHGTQS_diur = w_noHGTQS_mean.groupby(w_noHGTQS_mean.time.dt.hour)
w_noHGTQS_noSHAL_diur = w_noHGTQS_noSHAL_mean.groupby(w_noHGTQS_noSHAL_mean.time.dt.hour)
w_noHGTQS_noUVmix_diur = w_noHGTQS_noUVmix_mean.groupby(w_noHGTQS_noUVmix_mean.time.dt.hour)

# %%
'''Select the last level in order to get the surface vertical winds'''
#to get surface we first have to flip of select the last level

w_noHGTQS_surf = w_noHGTQS.sel(lev = 64)
w_noHGTQS_noSHAL_surf = w_noHGTQS_noSHAL.sel(lev = 64)
w_noHGTQS_noUVmix_surf = w_noHGTQS_noUVmix.sel(lev = 64)


w_noHGTQS_surf_mean = w_noHGTQS_surf.mean(dim = 'time')
w_noHGTQS_noSHAL_surf_mean = w_noHGTQS_noSHAL_surf.mean(dim = 'time')
w_noHGTQS_noUVmix_surf_mean = w_noHGTQS_noUVmix_surf.mean(dim = 'time')

#%%
'''Plot the graphs'''

levels = np.linspace(-0.4, 0.4, 9)

plt.figure()
plt.contourf(w_noHGTQS_surf_mean.lon.values, w_noHGTQS_surf_mean.lat.values, w_noHGTQS_surf_mean.wa.values, vmin = -0.4, vmax = 0.4, levels = levels)
color_bar = plt.colorbar(label = 'W wind speed [m/s]')
plt.title('Mean vertical wind speed (W) at the surface for HARMONIE noHGTQS')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
color_bar.minorticks_on()


plt.figure()
plt.contourf(w_noHGTQS_noSHAL_surf_mean.lon.values, w_noHGTQS_noSHAL_surf_mean.lat.values, w_noHGTQS_noSHAL_surf_mean.wa.values, vmin = -0.4, vmax = 0.4, levels = levels)
color_bar = plt.colorbar(label = 'W wind speed [m/s]')
plt.title('Mean vertical wind speed (W) at the surface for HARMONIE noHGTQS noSHAL')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
color_bar.minorticks_on()

plt.figure()
plt.contourf(w_noHGTQS_noUVmix_surf_mean.lon.values, w_noHGTQS_noUVmix_surf_mean.lat.values, w_noHGTQS_noUVmix_surf_mean.wa.values, vmin = -0.4, vmax = 0.4, levels = levels)
color_bar = plt.colorbar(label = 'W wind speed [m/s]')
plt.title('Mean vertical wind speed (W) at the surface for HARMONIE noHGTQS noUVmix')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
color_bar.minorticks_on()


#%%
'''Make cross section of vertical winds'''
# =============================================================================
# 
# w_noHGTQS_cross = w_noHGTQS.mean(dim = ('y', 'time'))
# w_noHGTQS_cross = w_noHGTQS_cross.transpose()
# plt.contourf(w_noHGTQS_cross.x.values, w_noHGTQS_cross.lev.values, w_noHGTQS_cross.wa.values)
# =============================================================================
