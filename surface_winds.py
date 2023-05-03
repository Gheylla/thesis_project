# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 08:26:42 2023

@author: gheylla
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
from cut_subgrid_module import cut_sub_grid
import datetime as dt


uas_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\uas_data\uas_noHGTQS_cut.nc')
uas_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\uas_data\uas_noHGTQS_noSHAL_cut.nc')
uas_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\uas_data\uas_noHGTQS_noUVmix_cut.nc')


vas_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\vas_data\vas_noHGTQS_cut.nc')
vas_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\vas_data\vas_noHGTQS_noSHAL_cut.nc')
vas_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\vas_data\vas_noHGTQS_noUVmix_cut.nc')


'''Calculate the total wind speed '''
# u and v are orthogonal to each other. The total wind speed becomes sqrt u2 + v2 

tot_speed_noHGTQS = np.sqrt(uas_noHGTQS.uas ** 2 + vas_noHGTQS.vas ** 2)
tot_speed_noHGTQS_noSHAL = np.sqrt(uas_noHGTQS_noSHAL.uas ** 2 + vas_noHGTQS_noSHAL.vas** 2)
tot_speed_noHGTQS_noUVmix = np.sqrt(uas_noHGTQS_noUVmix.uas ** 2 + vas_noHGTQS_noUVmix.vas ** 2)

mean_tot_speed_noHGTQS = tot_speed_noHGTQS.mean(dim = ('x', 'y'))
mean_tot_speed_noHGTQS_noSHAL = tot_speed_noHGTQS_noSHAL.mean(dim = ('x', 'y'))
mean_tot_speed_noHGTQS_noUVmix = tot_speed_noHGTQS_noUVmix.mean(dim = ('x', 'y'))


diur_noHGTQS = mean_tot_speed_noHGTQS.groupby(mean_tot_speed_noHGTQS.time.dt.hour).mean()
diur_noHGTQS_noSHAL = mean_tot_speed_noHGTQS_noSHAL.groupby(mean_tot_speed_noHGTQS_noSHAL.time.dt.hour).mean()
diur_noHGTQS_noUVmix = mean_tot_speed_noHGTQS_noUVmix.groupby(mean_tot_speed_noHGTQS_noUVmix.time.dt.hour).mean()


plt.figure(figsize = (15,10))
plt.plot(diur_noHGTQS.hour.values, diur_noHGTQS.values, label = 'HARMONIE noHGTQS')
plt.plot(diur_noHGTQS_noSHAL.hour.values, diur_noHGTQS_noSHAL.values, label = 'HARMONIE noHGTQS noSHAL')
plt.plot(diur_noHGTQS_noUVmix.hour.values, diur_noHGTQS_noUVmix.values, label = 'HARMONIE noHGTQS noUVmix')
plt.title('Composite diurnal cycle of total wind speed')
plt.xlabel('Time [hr]')
plt.ylabel('Wind speed [m/s]')
plt.grid()
plt.legend()

hours = np.linspace(0, 23, 24)
plt.xticks(diur_noHGTQS_noUVmix.hour.values, hours)

# =============================================================================
# mean_uas_noHGTQS = uas_noHGTQS.mean(dim = ('x', 'y'))
# mean_uas_noHGTQS = mean_uas_noHGTQS.groupby(uas_noHGTQS.time.dt.hour).mean()
# 
# mean_uas_noHGTQS_noSHAL = uas_noHGTQS_noSHAL.mean(dim = ('x', 'y'))
# mean_uas_noHGTQS_noSHAL = mean_uas_noHGTQS_noSHAL.groupby(uas_noHGTQS.time.dt.hour).mean()
# 
# mean_uas_noHGTQS_noUVmix = uas_noHGTQS_noUVmix.mean(dim = ('x', 'y'))
# mean_uas_noHGTQS_noUVmix = mean_uas_noHGTQS_noUVmix.groupby(uas_noHGTQS.time.dt.hour).mean()
# 
# 
# plt.figure()
# plt.plot(mean_uas_noHGTQS.hour.values, mean_uas_noHGTQS.uas.values, label = 'HARMONIE noHGTQS')
# plt.plot(mean_uas_noHGTQS_noSHAL.hour.values, mean_uas_noHGTQS_noSHAL.uas.values, label = 'HARMONIE noHGTQS noSHAL')
# plt.plot(mean_uas_noHGTQS_noUVmix.hour.values, mean_uas_noHGTQS_noUVmix.uas.values, label = 'HARMONIE noHGTQS noUVmix')
# 
# 
# plt.ylabel('Windspeed [m/s]')
# plt.title('East-ward near surface wind speed')
# plt.xlabel('Time [hr]')
# plt.grid()
# plt.legend()
# 
# 
# 
# 
# 
# mean_vas_noHGTQS = vas_noHGTQS.mean(dim = ('x', 'y'))
# mean_vas_noHGTQS = mean_vas_noHGTQS.groupby(vas_noHGTQS.time.dt.hour).mean()
# 
# mean_vas_noHGTQS_noSHAL = vas_noHGTQS_noSHAL.mean(dim = ('x', 'y'))
# mean_vas_noHGTQS_noSHAL = mean_vas_noHGTQS_noSHAL.groupby(vas_noHGTQS.time.dt.hour).mean()
# 
# mean_vas_noHGTQS_noUVmix = vas_noHGTQS_noUVmix.mean(dim = ('x', 'y'))
# mean_vas_noHGTQS_noUVmix = mean_vas_noHGTQS_noUVmix.groupby(vas_noHGTQS.time.dt.hour).mean()
# 
# 
# plt.figure()
# plt.plot(mean_vas_noHGTQS.hour.values, mean_vas_noHGTQS.vas.values, label = 'HARMONIE noHGTQS')
# plt.plot(mean_vas_noHGTQS_noSHAL.hour.values, mean_vas_noHGTQS_noSHAL.vas.values, label = 'HARMONIE noHGTQS noSHAL')
# plt.plot(mean_vas_noHGTQS_noUVmix.hour.values, mean_vas_noHGTQS_noUVmix.vas.values, label = 'HARMONIE noHGTQS UVmix')
# 
# 
# plt.ylabel('Windspeed [m/s]')
# plt.title('North-ward near surface wind speed')
# plt.xlabel('Time [hr]')
# plt.grid()
# plt.legend()
# =============================================================================
