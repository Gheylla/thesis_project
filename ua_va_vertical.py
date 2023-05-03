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

'''Plot the vertical distibution of the horizontal winds'''
#remember that because this is dependend on the level the area is different. The area is smaller than what we are actually analysing. 

ua_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\ua_data\ua_noHGTQS_cut.nc')
ua_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\ua_data\ua_noHGTQS_noSHAL_cut.nc')
ua_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\ua_data\ua_noHGTQS_noUVmix_cut.nc')

va_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\va_data\va_noHGTQS_cut.nc')
va_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\va_data\va_noHGTQS_noSHAL_cut.nc')
va_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\va_data\va_noHGTQS_noUVmix_cut.nc')

ua_noHGTQS = ua_noHGTQS.mean(dim = ('x', 'y'))
ua_noHGTQS_noSHAL = ua_noHGTQS_noSHAL.mean(dim = ('x', 'y'))
ua_noHGTQS_noUVmix = ua_noHGTQS_noUVmix.mean(dim = ('x', 'y'))

va_noHGTQS = va_noHGTQS.mean(dim = ('x', 'y'))
va_noHGTQS_noSHAL = va_noHGTQS_noSHAL.mean(dim = ('x', 'y'))
va_noHGTQS_noUVmix = va_noHGTQS_noUVmix.mean(dim = ('x', 'y'))

ua_noHGTQS = ua_noHGTQS.groupby(ua_noHGTQS.time.dt.hour).mean()
ua_noHGTQS_noSHAL = ua_noHGTQS_noSHAL.groupby(ua_noHGTQS_noSHAL.time.dt.hour).mean()
ua_noHGTQS_noUVmix = ua_noHGTQS_noUVmix.groupby(ua_noHGTQS_noUVmix.time.dt.hour).mean()

va_noHGTQS = va_noHGTQS.groupby(va_noHGTQS.time.dt.hour).mean()
va_noHGTQS_noSHAL = va_noHGTQS_noSHAL.groupby(va_noHGTQS_noSHAL.time.dt.hour).mean()
va_noHGTQS_noUVmix = va_noHGTQS_noUVmix.groupby(va_noHGTQS_noUVmix.time.dt.hour).mean()

plt.plot(va_noHGTQS.va.values[5], va_noHGTQS.lev.values)
plt.plot(va_noHGTQS_noSHAL.va.values[5], va_noHGTQS_noSHAL.lev.values)
plt.plot(va_noHGTQS_noUVmix.va.values[5], va_noHGTQS_noUVmix.lev.values)
# =============================================================================
# ua_hr4 = ua_noHGTQS.ua.values[4]
# ua_hr4 = ua_hr4[::-1]
# 
# plt.plot(ua_hr4[0:20])
# =============================================================================
