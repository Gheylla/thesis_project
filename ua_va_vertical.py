# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 11:13:55 2023

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

ua_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\ua_data\ua_noHGTQS_cut.nc')
ua_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\ua_data\ua_noHGTQS_noSHAL_cut.nc')
ua_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\ua_data\ua_noHGTQS_noUVmix_cut.nc')

va_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\va_data\va_noHGTQS_cut.nc')
va_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\va_data\va_noHGTQS_noSHAL_cut.nc')
va_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\va_data\va_noHGTQS_noUVmix_cut.nc')

heights = pd.read_csv(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\z_43_fullheights.csv')


# %%
'''Determine the mean for the area'''
ua_noHGTQS = ua_noHGTQS.mean(dim = ('x', 'y'))
ua_noHGTQS_noSHAL = ua_noHGTQS_noSHAL.mean(dim = ('x', 'y'))
ua_noHGTQS_noUVmix = ua_noHGTQS_noUVmix.mean(dim = ('x', 'y'))

va_noHGTQS = va_noHGTQS.mean(dim = ('x', 'y'))
va_noHGTQS_noSHAL = va_noHGTQS_noSHAL.mean(dim = ('x', 'y'))
va_noHGTQS_noUVmix = va_noHGTQS_noUVmix.mean(dim = ('x', 'y'))


# %%
'''Determine the mean of time (only left is lev)'''
ua_noHGTQS_totmean = ua_noHGTQS.mean(dim = 'time')
ua_noHGTQS_noSHAL_totmean = ua_noHGTQS_noSHAL.mean(dim = 'time')
ua_noHGTQS_noUVmix_totmean = ua_noHGTQS_noUVmix.mean(dim = 'time')


# %%
'''Flip the data in order to get the right level at the bottom'''
ua_noHGTQS_totmean_arr = np.array(ua_noHGTQS_totmean.ua.values)
ua_noHGTQS_noSHAL_totmean_arr = np.array(ua_noHGTQS_noSHAL_totmean.ua.values)
ua_noHGTQS_noUVmix_totmean_arr = np.array(ua_noHGTQS_noUVmix_totmean.ua.values)

ua_noHGTQS_totmean_arr = np.flip(ua_noHGTQS_totmean_arr)
ua_noHGTQS_noSHAL_totmean_arr = np.flip(ua_noHGTQS_noSHAL_totmean_arr)
ua_noHGTQS_noUVmix_totmean_arr = np.flip(ua_noHGTQS_noUVmix_totmean_arr)


# %%
'''Groupby hour in order to get the diurnal cycle'''
ua_noHGTQS_diur = ua_noHGTQS.groupby(ua_noHGTQS.time.dt.hour).mean()
ua_noHGTQS_noSHAL_diur = ua_noHGTQS_noSHAL.groupby(ua_noHGTQS_noSHAL.time.dt.hour).mean()
ua_noHGTQS_noUVmix_diur = ua_noHGTQS_noUVmix.groupby(ua_noHGTQS_noUVmix.time.dt.hour).mean()

va_noHGTQS_diur = va_noHGTQS.groupby(va_noHGTQS.time.dt.hour).mean()
va_noHGTQS_noSHAL_diur = va_noHGTQS_noSHAL.groupby(va_noHGTQS_noSHAL.time.dt.hour).mean()
va_noHGTQS_noUVmix_diur = va_noHGTQS_noUVmix.groupby(va_noHGTQS_noUVmix.time.dt.hour).mean()

# %%
'''Get the mean of only the hours that are the main difference in cloud cover (5-7 and 16-18)'''
ua_noHGTQS_diur_sel = ua_noHGTQS_diur.isel(hour = slice(5,8))
ua_noHGTQS_noSHAL_diur_sel = ua_noHGTQS_noSHAL_diur.isel(hour = slice(5,8))
ua_noHGTQS_noUVmix_diur_sel = ua_noHGTQS_noUVmix_diur.isel(hour = slice(5,8))

# %%
ua_noHGTQS_diur_sel_mean = ua_noHGTQS_diur_sel.mean(dim = 'hour')
ua_noHGTQS_noSHAL_diur_sel_mean = ua_noHGTQS_noSHAL_diur_sel.mean(dim = 'hour')
ua_noHGTQS_noUVmix_diur_sel_mean = ua_noHGTQS_noUVmix_diur_sel.mean(dim = 'hour')


# %%
ua_noHGTQS_diur_sel_mean_arr = np.array(ua_noHGTQS_diur_sel_mean.ua.values)
ua_noHGTQS_noSHAL_diur_sel_mean_arr = np.array(ua_noHGTQS_noSHAL_diur_sel_mean.ua.values)
ua_noHGTQS_noUVmix_diur_sel_mean_arr = np.array(ua_noHGTQS_noUVmix_diur_sel_mean.ua.values)

ua_noHGTQS_diur_sel_mean_arr = np.flip(ua_noHGTQS_diur_sel_mean_arr)
ua_noHGTQS_noSHAL_diur_sel_mean_arr = np.flip(ua_noHGTQS_noSHAL_diur_sel_mean_arr)
ua_noHGTQS_noUVmix_diur_sel_mean_arr = np.flip(ua_noHGTQS_noUVmix_diur_sel_mean_arr)


# %%
'''Plot the graphs'''
#this is the vertical profile of the mean of the area and the whole eurec4a period
plt.figure(figsize = (10,10))
plt.plot(ua_noHGTQS_totmean_arr, heights['Full heights'], label = 'HARMONIE noHGTQS')
plt.plot(ua_noHGTQS_noSHAL_totmean_arr, heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL')
plt.plot(ua_noHGTQS_noUVmix_totmean_arr, heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix')
plt.grid()
plt.legend()
plt.title('Vertical profile of Eastward Wind Velocity')
plt.xlabel('Eastward Wind Velocity [m/s]')
plt.ylabel('Height [m]')
plt.ylim([0, 5000])
plt.xlim([-10, -2])


#this plots the graphs for the specific hours that have the most deviations
plt.figure(figsize = (10,10))
plt.plot(ua_noHGTQS_diur_sel_mean_arr, heights['Full heights'], label = 'HARMONIE noHGTQS')
plt.plot(ua_noHGTQS_noSHAL_diur_sel_mean_arr, heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL')
plt.plot(ua_noHGTQS_noUVmix_diur_sel_mean_arr, heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix')
plt.grid()
plt.legend()
plt.title('Vertical profile of Eastward Wind Velocity')
plt.xlabel('Eastward Wind Velocity [m/s]')
plt.ylabel('Height [m]')
plt.ylim([0, 3000])
plt.xlim([-10, -2])


