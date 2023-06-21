# -*- coding: utf-8 -*-
"""
Created on Tue May 30 16:24:11 2023

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
import seaborn as sn


cll_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cll_data\cll_noHGTQS_cut.nc')
cll_noHGTQSnoSHAL =  xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cll_data\cll_noHGTQS_noSHAL_cut.nc')
cll_noHGTQSnoUVmix =  xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cll_data\cll_noHGTQS_noUVmix_cut.nc')

# %%
'''Plot the surface cape'''


cll_noHGTQS_arr = np.array(cll_noHGTQS.cll.values).reshape(-1)
cll_noHGTQS_noSHAL_arr = np.array(cll_noHGTQSnoSHAL.cll.values).reshape(-1)
cll_noHGTQS_noUVmix_arr = np.array(cll_noHGTQSnoUVmix.cll.values).reshape(-1)

#%%
plt.figure()
plt.hist(cll_noHGTQS_arr,  bins = 10, color = 'red', label = 'HARMONIE noHGTQS')
plt.legend()
plt.title("Histogram of CC")
plt.xlabel('CC [-]')
plt.ylabel('Count')
plt.grid()
plt.show()


plt.figure()
plt.hist(cll_noHGTQS_noSHAL_arr,   bins = 10, color = 'blue', label = 'HARMONIE noHGTQS noSHAL')
plt.legend()
plt.title("Histogram of CC")
plt.xlabel('CC [-]')
plt.ylabel('Count')
plt.grid()
plt.show()

plt.figure()
plt.hist(cll_noHGTQS_noUVmix_arr,   bins = 10, color = 'green', label = 'HARMONIE noHGTQS noUVmix')
plt.legend()
plt.title("Histogram of CC")
plt.xlabel('CC [-]')
plt.ylabel('Count')
plt.grid()
plt.show()
#%%
time = '2020-02-03T02:00:00.000000000'



cll_noHGTQS = cll_noHGTQS.sel(time = time)
cll_noHGTQS_noSHAL = cll_noHGTQSnoSHAL.sel(time = time)
cll_noHGTQS_noUVmix = cll_noHGTQSnoUVmix.sel(time = time)

cmap = 'Blues_r'
levels = np.linspace(0,1,11)


plt.figure(figsize = (10,10))
plt.contourf(cll_noHGTQS.lon.values, cll_noHGTQS.lat.values, cll_noHGTQS.cll.values, cmap = cmap, vmin = 0, vmax = 1, levels = levels)
plt.title(f'Cloud cover HARMONIE noHGTQS for time {time} [UTC]')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.colorbar(label = 'Cloud cover [-]')


plt.figure(figsize = (10,10))
plt.contourf(cll_noHGTQS_noSHAL.lon.values, cll_noHGTQS_noSHAL.lat.values, cll_noHGTQS_noSHAL.cll.values, cmap = cmap, vmin = 0, vmax = 1, levels = levels)
plt.title(f'Cloud cover HARMONIE noHGTQS noSHAL for time {time} [UTC]')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.colorbar(label = 'Cloud cover [-]')



plt.figure(figsize = (10,10))
plt.contourf(cll_noHGTQS_noUVmix.lon.values, cll_noHGTQS_noUVmix.lat.values, cll_noHGTQS_noUVmix.cll.values, cmap = cmap, vmin = 0, vmax = 1, levels = levels)
plt.title(f'Cloud cover HARMONIE noHGTQS noUVmix for time {time} [UTC]')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.colorbar(label = 'Cloud cover [-]')
