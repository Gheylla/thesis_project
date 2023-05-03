# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 09:43:21 2023

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
from skimage.transform import resize
from PIL import Image 
import netCDF4 as nc

'''Step 1: Import data'''
print('Step 1: Import data')
'''For this specific case we are using low cloud cover fields from HARMONIE outputs'''

path_dataset = r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cut_sat\OR_ABI-L2-CMIPF-M6C13_G16_s20200181110194_e20200181119514_c20200181120008.nc'
dataset = xr.open_mfdataset(path_dataset, combine = 'by_coords', engine = 'netcdf4')



#with x and y you get the whole space filled with values, with lat and lon you get nan values
plt.figure()
dataset = dataset.isel(time = 0)

scene = dataset.C13_Brightness_Temp
plt.contourf(dataset.lon.values, dataset.lat.values, dataset.C13_Brightness_Temp.values)

imin = np.argmin(scene.shape)
scene = resize(scene, (scene.shape[imin],scene.shape[imin]), anti_aliasing=True)

assert(np.isnan(scene).any() == False)


vmin = 270
vmax = 300
cmap = "RdBu_r"
plt.figure()
plt.imshow(scene, cmap=cmap, vmin=vmin, vmax=vmax)
plt.show()
# =============================================================================
# 
# plt.figure()
# plt.contourf(dataset.y.values, dataset.x.values, dataset.C13_Brightness_Temp.values)
# #dataset.C13_Brightness_Temp.plot(x = dataset.lon.values, y = dataset.lat.values)
# =============================================================================

# =============================================================================
# plt.scatter()
# =============================================================================
