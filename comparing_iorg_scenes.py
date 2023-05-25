# -*- coding: utf-8 -*-
"""
Created on Thu May 25 11:01:54 2023

@author: Gheylla Liberia
this code tries plots 2 scenes from HARMONIE that have different Iorg values to see what is happening with the Iorg
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import cloudmetrics
import inspect
from skimage.transform import resize
import seaborn as sn

#%% 
'''Import data'''
harm_nohgtqs = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cll_data\cll_noHGTQS_cut.nc')

#%%
#from the timeseries of iorg we check which days have the most different iorg valeues and plot these scenes to get a feeling of what is happening
#min iorg at time 2020-02-07 08:00:00
#max iorg at time '2020-02-01 07:00:00


harm_nohgtqs_min = harm_nohgtqs.sel(time = '2020-02-07 08:00:00')
harm_nohgtqs_max = harm_nohgtqs.sel(time = '2020-01-28 21:00:00')

#%%
'''plot the data'''
plt.figure(figsize = (10,10))
plt.contourf(harm_nohgtqs_min.lon.values, harm_nohgtqs_min.lat.values, harm_nohgtqs_min.cll.values, cmap = 'Blues_r')
plt.title(f'Low cloud cover for time {harm_nohgtqs_min.time.values} with iorg of 0.414')
plt.xlabel('Longitude [degrees]')
plt.ylabel('Latitude [degrees]')
plt.colorbar(label = 'CF [-]')


plt.figure(figsize = (10,10))
plt.contourf(harm_nohgtqs_max.lon.values, harm_nohgtqs_max.lat.values, harm_nohgtqs_max.cll.values, cmap = 'Blues_r')
plt.title(f'Low cloud cover for time {harm_nohgtqs_max.time.values} with iorg of 0.875')
plt.xlabel('Longitude [degrees]')
plt.ylabel('Latitude [degrees]')
plt.colorbar(label = 'CF [-]')