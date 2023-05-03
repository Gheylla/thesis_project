# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 11:00:28 2023

@author: LENOVO
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 14:45:25 2023

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

t0 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\metrics_obs_292286_new.h5')
t1 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\metrics_obs_293280_new.h5')
t2 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\metrics_obs_293281_new.h5')
t3 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\metrics_obs_293282_new.h5')
t4 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\metrics_obs_293283_new.h5')
t5 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\metrics_obs_293284_new.h5')
t6 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\metrics_obs_293285_new.h5')
t7 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\metrics_obs_293286_new.h5')


# =============================================================================
# 
# data = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cut_sat\OR_ABI-L2-CMIPF-M6C13_G16_s20200181110194_e20200181119514_c20200181120008.nc')
# data = data.isel(time =0)
# 
# data.C13_Brightness_Temp.plot(cmap = 'Blues')
# =============================================================================
