# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 10:06:11 2023

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
import os
import eurec4a
from skimage.transform import resize

'''Cloud base height from mesonet data'''
data = pd.read_csv(r'C:\Users\LENOVO\Desktop\TU_Delft\TBPB.csv')
#this data was found in the paper satellite retrieval of cloudbase height and geometric thickness of low-level cloud based on CALIPSO. 
#the data is a cloud base level 1 data and is retrieved from barbados and can be found on the website https://mesonet.agron.iastate.edu/request/download.phtml?network=BB__ASOS

for i in range(len(data['skyl1'])):
    if data['skyl1'][i] == 'M':
        data = data.drop(i)
        
data['skyl1'] = data['skyl1'].astype(float)
conversion = 3.281
data['skyl1_meters'] = data['skyl1'] / conversion
mean_cbh = data['skyl1_meters'].mean()


'''Cloud base height from eurec4a'''

from intake import open_catalog
cat = open_catalog("https://raw.githubusercontent.com/eurec4a/eurec4a-intake/master/catalog.yml")

cat = eurec4a.get_intake_catalog()

#data = xr.concat([cat.barbados.bco.ceilometer['latlongrid'](date=d).to_dask().chunk() for d in pd.date_range('2020-01-12','2020-01-20')], dim='time')
data_ceilo = cat.barbados.bco.ceilometer.CBH.to_dask()
ceilo_val = data_ceilo.cbh_1.values
ceilo_val = ceilo_val[~np.isnan(ceilo_val)]
mean_cbh_ceilo = np.mean(ceilo_val)


'''Plot the different cloud bases'''
# =============================================================================
# plt.figure(figsize = (15,10))
# plt.plot( ceilo_val)
# 
# plt.figure(figsize = (15,10))
# plt.plot(data['skyl1_meters'])
# =============================================================================
