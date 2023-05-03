# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 11:26:43 2023

@author: Gheylla Liberia
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

lwp = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\clw_Slev_his_EUREC4Acircle_HA43h22tg3_clim_noHGTQS_1hr_202001210000-202002010000.nc')

#plot all of the levels separately but make them transclucent and place them on top of eachother


lwp = lwp.isel(time = 3)
lwp = lwp.isel(lev = 60)
levels = [-3.38e-21, 1]
plt.contourf(lwp.lon.values, lwp.lat.values, (lwp.clw.values * 1000), vmin = np.min(lwp.clw.values), vmax = np.max(lwp.clw.values))

for i in range(len(lwp.time.values)):
    #print(i)
    lwp_cut = lwp.isel(time = i)
    #print(lwp_cut)
    for ii in range(len(lwp_cut.lev.values)):
        lwp_cut_lev = lwp_cut.isel(lev = ii)
        plt.figure()
        plt.contour(lwp_cut_lev.lat.values, lwp_cut_lev.lon.values, lwp_cut_lev.clw.values)
    #print(lwp_cut.time.values)
    #print(lwp_cut.lev.values)
    #plt.figure()
    #plt.plot()
    #print(lwp1.time.values)

