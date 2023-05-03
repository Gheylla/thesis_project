# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 13:35:52 2023

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
from cut_subgrid_module import cut_sub_grid


cll = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cll_data\cll_*.nc', combine = 'by_coords')

cll = cut_sub_grid(cll)


x = len(cll.x)
y = len(cll.y)
imin = np.min(cll.cll.shape)
#imin = cll.cll.shape[0]
cll = cll.isel(x = slice(0, (imin)), y = slice(0, (imin)))

cll.to_netcdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cll_data\cll_noHGTQS_noUVmix_cut.nc')

y = 13.1939
x = -59.5432

cll_cut = cll.isel(time = 1)
plt.figure(figsize = (15,10))
levels = np.linspace(0,1,10)
plt.contourf(cll_cut.lon.values, cll_cut.lat.values, (cll_cut.cll.values), cmap = 'Blues_r', vmin = 0, vmax = 1, levels = levels)
plt.plot(x, y, marker="o", markersize=20, markeredgecolor="red", markerfacecolor="red")
plt.text((x + 0.2), (y + 0.2), 'BCO')
plt.title(f'Low Cloud Cover for HARMONIE noHGTQS at time = {cll_cut.time.values} ')
plt.xlabel('Longitude [degrees]')
plt.ylabel('Latitude [degrees]')
plt.colorbar(label = 'Low Cloud Cover [-]')

# =============================================================================
# 
# for i in range(len(cll.time.values)):
#     cll_cut = cll.isel(time = i)
#     plt.figure(figsize = (15,10))
#     levels = np.linspace(0,1,10)
#     plt.contourf(cll_cut.lon.values, cll_cut.lat.values, (cll_cut.cll.values), cmap = 'Blues_r', vmin = 0, vmax = 1, levels = levels)
#     plt.plot(x, y, marker="o", markersize=20, markeredgecolor="red", markerfacecolor="red")
#     plt.text((x + 0.2), (y + 0.2), 'BCO')
#     plt.title(f'Low Cloud Cover for HARMONIE noHGTQS noUVmix at time = {cll_cut.time.values} ')
#     plt.xlabel('Longitude [degrees]')
#     plt.ylabel('Latitude [degrees]')
#     plt.colorbar(label = 'Low Cloud Cover [-]')
#     
# 
# =============================================================================

