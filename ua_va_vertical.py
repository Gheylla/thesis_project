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


#%% 
'''Get the mean vertical profile of the whole period'''

ua_noHGTQS_mean_tot = ua_noHGTQS.mean(dim = ('time'))
ua_noHGTQS_noSHAL_mean_tot = ua_noHGTQS_noSHAL.mean(dim = ('time'))
ua_noHGTQS_noUVmix_mean_tot = ua_noHGTQS_noUVmix.mean(dim = ('time'))

va_noHGTQS_mean_tot = va_noHGTQS.mean(dim = ('time'))
va_noHGTQS_noSHAL_mean_tot = va_noHGTQS_noSHAL.mean(dim = ('time'))
va_noHGTQS_noUVmix_mean_tot = va_noHGTQS_noUVmix.mean(dim = ('time'))

#%%
'''Plot the mean vertical profile of the winds'''
plt.figure(figsize = (10,10))
plt.plot(np.flip(ua_noHGTQS_mean_tot.ua.values), heights['Full heights'], color = 'red', label = 'HARMONIE noHGTQS')
plt.plot(np.flip(ua_noHGTQS_noSHAL_mean_tot.ua.values), heights['Full heights'], color = 'blue', label = 'HARMONIE noHGTQS noSHAL' )
plt.plot(np.flip(ua_noHGTQS_noUVmix_mean_tot.ua.values), heights['Full heights'], color = 'green' , label = 'HARMONIE noHGTQS noUVmix')
plt.ylim([0,3000])
plt.xlim([-11, -2.5])
plt.legend()
plt.grid()
plt.title('Mean vertical profile of east-ward winds')
plt.xlabel('U windspeed [m/s]')
plt.ylabel('Height [m]')

#%%

'''Plot the mean vertical profile of the winds'''
plt.figure(figsize = (10,10))
plt.plot(np.flip(va_noHGTQS_mean_tot.va.values), heights['Full heights'], color = 'red', label = 'HARMONIE noHGTQS')
plt.plot(np.flip(va_noHGTQS_noSHAL_mean_tot.va.values), heights['Full heights'], color = 'blue', label = 'HARMONIE noHGTQS noSHAL' )
plt.plot(np.flip(va_noHGTQS_noUVmix_mean_tot.va.values), heights['Full heights'], color = 'green' , label = 'HARMONIE noHGTQS noUVmix')
plt.ylim([0,3000])
plt.xlim([-11, -2.5])
plt.legend()
plt.grid()
plt.title('Mean vertical profile of north-ward winds')
plt.xlabel('V windspeed [m/s]')
plt.ylabel('Height [m]')


# %%
'''Determine the diurnal cycle of the winds'''
ua_noHGTQS_diur = ua_noHGTQS.groupby(ua_noHGTQS.time.dt.hour).mean()
ua_noHGTQS_noSHAL_diur = ua_noHGTQS_noSHAL.groupby(ua_noHGTQS_noSHAL.time.dt.hour).mean()
ua_noHGTQS_noUVmix_diur = ua_noHGTQS_noUVmix.groupby(ua_noHGTQS_noUVmix.time.dt.hour).mean()

va_noHGTQS_diur = va_noHGTQS.groupby(va_noHGTQS.time.dt.hour).mean()
va_noHGTQS_noSHAL_diur = va_noHGTQS_noSHAL.groupby(va_noHGTQS.time.dt.hour).mean()
va_noHGTQS_noUVmix_diur = va_noHGTQS_noUVmix.groupby(va_noHGTQS.time.dt.hour).mean()

#transpose the table so we can make contourplot
ua_noHGTQS_diur = ua_noHGTQS_diur.transpose()
ua_noHGTQS_noSHAL_diur = ua_noHGTQS_noSHAL_diur.transpose()
ua_noHGTQS_noUVmix_diur = ua_noHGTQS_noUVmix_diur.transpose()

va_noHGTQS_diur = va_noHGTQS_diur.transpose()
va_noHGTQS_noSHAL_diur = va_noHGTQS_noSHAL_diur.transpose()
va_noHGTQS_noUVmix_diur = va_noHGTQS_noUVmix_diur.transpose()


#%%
'''Plot the diurnal cycle of u winds at all levels in order to make the video (inspired by Stephan presentation)'''
for i in range(len(ua_noHGTQS_diur.hour.values)):
    print(i)
    ua_hour_noHGTQS = ua_noHGTQS_diur.sel(hour = i)
    ua_hour_noHGTQS_noSHAL = ua_noHGTQS_noSHAL_diur.sel(hour = i)
    ua_hour_noHGTQS_noUVmix = ua_noHGTQS_noUVmix_diur.sel(hour = i)
    plt.figure(figsize = (10,10))
    plt.plot(np.flip(ua_hour_noHGTQS.ua.values), heights['Full heights'], color = 'red', label = 'HARMONIE noHGTQS')
    plt.plot(np.flip(ua_hour_noHGTQS_noSHAL.ua.values), heights['Full heights'], color = 'blue', label = 'HARMONIE noHGTQS noSHAL')
    plt.plot(np.flip(ua_hour_noHGTQS_noUVmix.ua.values), heights['Full heights'], color = 'green', label = 'HARMONIE noHGTQS noUVmix')
    plt.ylim([0,3000])
    plt.xlim([-12, 2.5])
    plt.title(f'Vertical profile of U winds for hour {ua_noHGTQS_diur.hour.values[i]}   [UTC] ')
    plt.xlabel('U windspeed [m/s]')
    plt.ylabel('Height [m]')
    plt.grid()
    plt.legend()
    plt.show()
    
# %%
'''Plot the diurnal cycle of v winds at all levels in order to make the video (inspired by Stephan presentation)'''
for i in range(len(ua_noHGTQS_diur.hour.values)):
    print(i)
    va_hour_noHGTQS = va_noHGTQS_diur.sel(hour = i)
    va_hour_noHGTQS_noSHAL = va_noHGTQS_noSHAL_diur.sel(hour = i)
    va_hour_noHGTQS_noUVmix = va_noHGTQS_noUVmix_diur.sel(hour = i)
    plt.figure(figsize = (10,10))
    plt.plot(np.flip(va_hour_noHGTQS.va.values), heights['Full heights'], color = 'red', label = 'HARMONIE noHGTQS', linestyle = 'dashed')
    plt.plot(np.flip(va_hour_noHGTQS_noSHAL.va.values), heights['Full heights'], color = 'blue', label = 'HARMONIE noHGTQS noSHAL', linestyle = 'dashed')
    plt.plot(np.flip(va_hour_noHGTQS_noUVmix.va.values), heights['Full heights'], color = 'green', label = 'HARMONIE noHGTQS noUVmix', linestyle = 'dashed')
    plt.ylim([0,3000])
    plt.xlim([-4, 2])
    plt.title(f'Vertical profile of V winds for hour {va_noHGTQS_diur.hour.values[i]}   [UTC] ')
    plt.xlabel('V windspeed [m/s]')
    plt.ylabel('Height [m]')
    plt.grid()
    plt.legend()    
    plt.show()
    
    
    
    
    
    
    
    
    

# %% 
'''Plot contour of u winds'''
plt.figure(figsize = (15,10))
plt.contourf(ua_noHGTQS_diur.hour.values, heights['Full heights'], ua_noHGTQS_diur.ua.values, colormap = 'jet')
plt.title('Composite diurnal cycle of Eastward winds HARMONIE noHGTQS')
plt.xlabel('UTC Time [hr]')
plt.ylabel('Height [m]')
plt.ylim([0,3000])
plt.colorbar(label = 'u windspeed [m/s]')

plt.figure(figsize = (15,10))
plt.contourf(ua_noHGTQS_noSHAL_diur.hour.values, heights['Full heights'], ua_noHGTQS_noSHAL_diur.ua.values, colormap = 'jet')
plt.title('Composite diurnal cycle of Eastward winds HARMONIE noHGTQS noSHAL')
plt.xlabel('UTC Time [hr]')
plt.ylabel('Height [m]')
plt.ylim([0,3000])
plt.colorbar(label = 'u windspeed [m/s]')

plt.figure(figsize = (15,10))
plt.contourf(ua_noHGTQS_noUVmix_diur.hour.values, heights['Full heights'], ua_noHGTQS_noUVmix_diur.ua.values, colormap = 'jet')
plt.title('Composite diurnal cycle of Eastward winds HARMONIE noHGTQS noUVmix')
plt.xlabel('UTC Time [hr]')
plt.ylabel('Height [m]')
plt.ylim([0,3000])
plt.colorbar(label = 'u windspeed [m/s]')

# %%
'''Plot contour of v winds'''
plt.figure(figsize = (15,10))
plt.contourf(va_noHGTQS_diur.hour.values, heights['Full heights'], va_noHGTQS_diur.va.values, colormap = 'jet')
plt.title('Composite diurnal cycle of Eastward winds HARMONIE noHGTQS')
plt.xlabel('UTC Time [hr]')
plt.ylabel('Height [m]')
plt.ylim([0,3000])
plt.colorbar(label = 'v windspeed [m/s]')

plt.figure(figsize = (15,10))
plt.contourf(va_noHGTQS_noSHAL_diur.hour.values, heights['Full heights'], va_noHGTQS_noSHAL_diur.va.values, colormap = 'jet')
plt.title('Composite diurnal cycle of Eastward winds HARMONIE noHGTQS noSHAL')
plt.xlabel('UTC Time [hr]')
plt.ylabel('Height [m]')
plt.ylim([0,3000])
plt.colorbar(label = 'v windspeed [m/s]')

plt.figure(figsize = (15,10))
plt.contourf(va_noHGTQS_noUVmix_diur.hour.values, heights['Full heights'], va_noHGTQS_noUVmix_diur.va.values, colormap = 'jet')
plt.title('Composite diurnal cycle of Eastward winds HARMONIE noHGTQS noUVmix')
plt.xlabel('UTC Time [hr]')
plt.ylabel('Height [m]')
plt.ylim([0,3000])
plt.colorbar(label = 'v windspeed [m/s]')

# %% 
'''Select hours that deviate the most'''
#first for ua
ua_noHGTQS_567 = ua_noHGTQS_diur.sel(hour = slice(5,7))
ua_noHGTQS_noSHAL_567 = ua_noHGTQS_noSHAL_diur.sel(hour = slice(5,7))
ua_noHGTQS_noUVmix_567 = ua_noHGTQS_noUVmix_diur.sel(hour = slice(5,7))


ua_noHGTQS_567mean = ua_noHGTQS_567.mean(dim = 'hour')
ua_noHGTQS_noSHAL_567mean = ua_noHGTQS_noSHAL_567.mean(dim = 'hour')
ua_noHGTQS_noUVmix_567mean = ua_noHGTQS_noUVmix_567.mean(dim = 'hour')


#for va
va_noHGTQS_567 = va_noHGTQS_diur.sel(hour = slice(5,7))
va_noHGTQS_noSHAL_567 = va_noHGTQS_noSHAL_diur.sel(hour = slice(5,7))
va_noHGTQS_noUVmix_567 = va_noHGTQS_noUVmix_diur.sel(hour = slice(5,7))


va_noHGTQS_567mean = va_noHGTQS_567.mean(dim = 'hour')
va_noHGTQS_noSHAL_567mean = va_noHGTQS_noSHAL_567.mean(dim = 'hour')
va_noHGTQS_noUVmix_567mean = va_noHGTQS_noUVmix_567.mean(dim = 'hour')


#now for hours 17-19
ua_noHGTQS_1719 = ua_noHGTQS_diur.sel(hour = slice(17,19))
ua_noHGTQS_noSHAL_1719 = ua_noHGTQS_noSHAL_diur.sel(hour = slice(17,19))
ua_noHGTQS_noUVmix_1719 = ua_noHGTQS_noUVmix_diur.sel(hour = slice(17,19))


ua_noHGTQS_1719mean = ua_noHGTQS_1719.mean(dim = 'hour')
ua_noHGTQS_noSHAL_1719mean = ua_noHGTQS_noSHAL_1719.mean(dim = 'hour')
ua_noHGTQS_noUVmix_1719mean = ua_noHGTQS_noUVmix_1719.mean(dim = 'hour')


#for va
va_noHGTQS_1719 = va_noHGTQS_diur.sel(hour = slice(17,19))
va_noHGTQS_noSHAL_1719 = va_noHGTQS_noSHAL_diur.sel(hour = slice(17,19))
va_noHGTQS_noUVmix_1719 = va_noHGTQS_noUVmix_diur.sel(hour = slice(17,19))


va_noHGTQS_1719mean = va_noHGTQS_1719.mean(dim = 'hour')
va_noHGTQS_noSHAL_1719mean = va_noHGTQS_noSHAL_1719.mean(dim = 'hour')
va_noHGTQS_noUVmix_1719mean = va_noHGTQS_noUVmix_1719.mean(dim = 'hour')


# %%
'''Plot the vertical profiles of the winds at the hours where the cloud cover deviates the most'''
plt.figure(figsize = (10,10))
plt.plot(np.flip(ua_noHGTQS_567mean.ua.values), heights['Full heights'], label = 'HARMONIE noHGTQS', color = 'red')
plt.plot(np.flip(ua_noHGTQS_noSHAL_567mean.ua.values), heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL', color = 'blue')
plt.plot(np.flip(ua_noHGTQS_noUVmix_567mean.ua.values), heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix', color = 'green')

plt.title('Mean Eastward vertical wind profiles for hours 5, 6 and 7 [UTC]')
plt.ylabel('Height [m]')
plt.xlabel('Windspeed [m/s]')
plt.grid()
plt.ylim([0,3000])
plt.xlim([-11, -2.5])
plt.legend()

# =============================================================================
# plt.figure(figsize = (10,10))
# plt.plot(np.flip(va_noHGTQS_567mean.va.values), heights['Full heights'], label = 'HARMONIE noHGTQS', color = 'red')
# plt.plot(np.flip(va_noHGTQS_noSHAL_567mean.va.values), heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL', color = 'blue')
# plt.plot(np.flip(va_noHGTQS_noUVmix_567mean.va.values), heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix', color = 'green')
# 
# plt.title('Mean Northward vertical wind profiles for hours 5, 6 and 7 [UTC]')
# plt.ylabel('Height [m]')
# plt.xlabel('Windspeed [m/s]')
# plt.grid()
# plt.ylim([0,3000])
# plt.legend()
# =============================================================================


plt.figure(figsize = (10,10))
plt.plot(np.flip(ua_noHGTQS_1719mean.ua.values), heights['Full heights'], label = 'HARMONIE noHGTQS', color = 'red')
plt.plot(np.flip(ua_noHGTQS_noSHAL_1719mean.ua.values), heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL', color = 'blue')
plt.plot(np.flip(ua_noHGTQS_noUVmix_1719mean.ua.values), heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix', color = 'green')

plt.title('Mean Eastward vertical wind profiles for hours 17, 18 and 19 [UTC]')
plt.ylabel('Height [m]')
plt.xlabel('Windspeed [m/s]')
plt.grid()
plt.ylim([0,3000])
plt.xlim([-11, -2.5])
plt.legend()

# =============================================================================
# plt.figure(figsize = (10,10))
# plt.plot(np.flip(va_noHGTQS_1719mean.va.values), heights['Full heights'], label = 'HARMONIE noHGTQS', color = 'red')
# plt.plot(np.flip(va_noHGTQS_noSHAL_1719mean.va.values), heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL', color = 'blue')
# plt.plot(np.flip(va_noHGTQS_noUVmix_1719mean.va.values), heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix', color = 'green')
# 
# plt.title('Mean Northward vertical wind profiles for hours 17, 18 and 19 [UTC]')
# plt.ylabel('Height [m]')
# plt.xlabel('Windspeed [m/s]')
# plt.grid()
# plt.ylim([0,3000])
# plt.legend()
# 
# =============================================================================




























