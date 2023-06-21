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

heights = pd.read_csv(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\z_43_fullheights.csv')

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
'''We want to get the snapshot of one level so determine the mean of time'''
w_noHGTQS_meant = w_noHGTQS.mean(dim = ('time'))
w_noHGTQS_noSHAL_meant = w_noHGTQS_noSHAL.mean(dim = ('time'))
w_noHGTQS_noUVmix_meant = w_noHGTQS_noUVmix.mean(dim = ('time'))

#%%
'''Get a snapshot of w winds at around 2500m'''
#loc_2500 is around index 29 which is equal to evel 36
w_noHGTQS_2500 = w_noHGTQS.isel(lev = 36)
w_noHGTQS_noSHAL_2500 = w_noHGTQS_noSHAL.isel(lev = 36)
w_noHGTQS_noUVmix_2500 = w_noHGTQS_noUVmix.isel(lev = 36)

#%%
'''Plot the w wind at around 2500'''
levels = np.linspace(-0.5, 2, 11)

for i in range(len(w_noHGTQS_noUVmix_2500.time.values)):
    w_noHGTQS_noUVmix_2500_timesel = w_noHGTQS_noUVmix_2500.isel(time = i)
    plt.contourf(w_noHGTQS_noUVmix_2500_timesel.lon.values, w_noHGTQS_noUVmix_2500_timesel.lat.values, w_noHGTQS_noUVmix_2500_timesel.wa.values, vmin = -0.3, vmax = 1.5, levels = levels)
    plt.title(f'W wind at height 2413m for HARMONIE noHGTQS noUVmix for time {w_noHGTQS_noUVmix_2500_timesel.time.values}')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.colorbar(label = 'W windspeed [m/s]')
    plt.show()
    

#%%
'''Get the distribution of vertical velocity'''
w_nohgtqs_arr = np.array(w_noHGTQS.wa.values).reshape(-1)
w_nohgtqs_noshal_arr = np.array(w_noHGTQS_noSHAL.wa.values).reshape(-1)
w_nohgtqs_nouvmix_arr = np.array(w_noHGTQS_noUVmix.wa.values).reshape(-1)

#%%
'''Make histogram'''
bin_bounds = [-1.5, -1.25, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5]
plt.figure()
plt.hist(w_nohgtqs_arr ,   bins =bin_bounds , color = 'red', label = 'HARMONIE noHGTQS', density = True, histtype = 'step', align = 'mid')
plt.hist(w_nohgtqs_noshal_arr ,   bins = bin_bounds, color = 'blue', label = 'HARMONIE noHGTQS noSHAL', density = True, histtype = 'step', align = 'mid')
plt.hist(w_nohgtqs_nouvmix_arr ,   bins = bin_bounds, color = 'green', label = 'HARMONIE noHGTQS noUVmix', density = True, histtype = 'step', align = 'mid')
plt.legend()
plt.title("Histogram of vertical wind speed")
plt.xlabel('W windspeed [m/s]')
plt.ylabel('Density')
plt.grid()
plt.show()






# =============================================================================
# plt.figure()
# plt.contourf(w_noHGTQS_2500.lon.values, w_noHGTQS_2500.lat.values, w_noHGTQS_2500.wa.values, vmin = -0.090, vmax = 0.09, levels = levels)
# plt.title('W wind at height 2413m for HARMONIE noHGTQS')
# plt.xlabel('Longitude')
# plt.ylabel('Latitude')
# plt.colorbar(label = 'W windspeed [m/s]')
# plt.show()
# 
# plt.figure()
# plt.contourf(w_noHGTQS_noSHAL_2500.lon.values, w_noHGTQS_noSHAL_2500.lat.values, w_noHGTQS_noSHAL_2500.wa.values, vmin = -0.090, vmax = 0.09, levels = levels)
# plt.title('W wind at height 2413m for HARMONIE noHGTQS noSHAL')
# plt.xlabel('Longitude')
# plt.ylabel('Latitude')
# plt.colorbar(label = 'W windspeed [m/s]')
# plt.show()
# 
# plt.figure()
# plt.contourf(w_noHGTQS_noUVmix_2500.lon.values, w_noHGTQS_noUVmix_2500.lat.values, w_noHGTQS_noUVmix_2500.wa.values, vmin = -0.090, vmax = 0.09, levels = levels)
# plt.title('W wind at height 2413m for HARMONIE noHGTQS noUVmix')
# plt.xlabel('Longitude')
# plt.ylabel('Latitude')
# plt.colorbar(label = 'W windspeed [m/s]')
# plt.show()
# =============================================================================
