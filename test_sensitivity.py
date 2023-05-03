# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 15:16:17 2023

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
import pyhdf
#import GOES
from calculating_latlon import calculate_degrees
from cut_subgrid_module import cut_sub_grid
import os


# =============================================================================
# data_dir = r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cut_sat\OR_ABI-L2-CMIPF-M6C13_G16_s20200530630188_e20200530639508_c20200530639595.nc'
# #entries = os.listdir(directory)
# 
# 
# data = xr.open_mfdataset(data_dir, engine = 'netcdf4')
# data = data.isel(time = 0)
# plt.contourf(data.lon.values, data.lat.values, data.C13_Brightness_Temp.values, cmap = 'Blues')
# plt.colorbar()
# plt.legend()
# 
# plt.title(f'Cloud cover for day {data.time}')
# plt.xlabel('Longitude')
# plt.ylabel('Latitude')
# =============================================================================


metrics_290280 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\metrics_obs_test290280.h5')
metrics_295280 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\metrics_obs_test295280.h5')
metrics_300280 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\metrics_obs_test300280.h5')
metrics_290270 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\metrics_obs_test290270.h5')
metrics_290260 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\metrics_obs_test290260.h5')


plt.figure(figsize = (10,5))
plt.plot(metrics_290280.Time.values ,metrics_290280['mean_length_scale'], marker = 'o', label = '290-280')
#plt.plot(metrics_295280.Time.values ,metrics_295280['mean_length_scale'], marker = 'o', label = '295-280')
#plt.plot(metrics_300280.Time.values ,metrics_300280['mean_length_scale'], marker = 'o', label = '300-280')
plt.plot(metrics_290270.Time.values ,metrics_290270['mean_length_scale'], marker = 'o', label = '290-270')
#plt.plot(metrics_290260.Time.values ,metrics_290260['cloud_fraction'], marker = 'o', label = '290-260')

plt.legend()
plt.grid()
plt.title('Mean length scale')
plt.xlabel('Date')
plt.ylabel('Mean length scale[-]')
