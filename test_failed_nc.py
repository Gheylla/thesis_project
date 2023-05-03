# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 16:26:54 2023

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
# data = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cut_sat\OR_ABI-L2-CMIPF-M6C13_G16_s20200101600212_e20200101609532_c20200101610017.nc')
# data = data.isel(time = 0)
# 
# =============================================================================

directory = r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cut_sat'
entries = os.listdir(directory)

for entry in entries:
    file_location_name = os.path.join(directory,entry)
    data = xr.open_mfdataset(file_location_name)
    data = data.isel(time = 0)
    plt.figure()
    plt.contourf(data.lon.values, data.lat.values, data.C13_Brightness_Temp.values, cmap = 'Blues_r')


# =============================================================================
# for i in range(len(data.C13_Brightness_Temp.values[0,:])):
#     print(i)
#     for ii in range(len(data.C13_Brightness_Temp.values[:,0])):
#         print(ii)
#         if np.isnan(data.C13_Brightness_Temp.values[i,ii]) == False: 
#             data.C13_Brightness_Temp.values[i,ii] = 0
# 
# =============================================================================
