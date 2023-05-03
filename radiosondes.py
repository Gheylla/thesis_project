# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 17:49:29 2023

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
#import pyhdf
#import GOES
import os
from numpy import nan

data = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\BCO_data\EUREC4A_BCO_Vaisala-RS_L2_v3.0.0.nc')


# =============================================================================
# rel_hum = np.linspace(0, 1, 11)
# temp = np.linspace(0,3500, 101)
# loc_cb = []
# 
# for i in range(len(data.sounding) - 330):
#     data_sel = data.isel(sounding = i)
#     plt.figure(figsize = (15, 35))
#     plt.plot(data_sel.dp.values[0:350], data_sel.alt.values[0:350], label = 'Dew temperature')
#     plt.plot(data_sel.ta.values[0:350], data_sel.alt.values[0:350], label = 'Dry bulb temperature')
#     plt.legend()
#     plt.yticks(temp)
#     plt.grid()
#     plt.axhline(y = 600, color = 'r', linestyle = '-')
#     
# =============================================================================
    
    #print(data.launch_time.values)
    
    
loc_base = np.where(data.alt.values == 2000)[0][0]
temp_base = []

for i in range(len(data.sounding)):
    data_new = data.isel(sounding = i)
    data_new = data_new.ta.values.astype(str)
    if data_new[loc_base] != 'nan':
        print('not nan')
        data_new = data_new.astype(float)
        temp_base.append(data_new[loc_base])
    else:
        print('it has nan')
    
    

#temp_base= temp_base[~np.isnan(temp_base)]


