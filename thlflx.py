# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 17:58:40 2023

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


#%%
'''Import data'''
conv_dry_nohgtqs = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\thflx_data\thflx_conv_dry_noHGTQS.nc')
conv_dry_nohgtqs_noshal = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\thflx_data\thflx_conv_dry_noHGTQS_noSHAL.nc')
conv_moist_nohgtqs = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\thflx_data\thflx_conv_moist_noHGTQS.nc')
conv_moist_nohgtqs_noshal = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\thflx_data\thflx_conv_moist_noHGTQS_noSHAL.nc')
turb_nohgtqs = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\thflx_data\thflx_turb_noHGTQS.nc')
turb_nohgtqs_noshal = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\thflx_data\thflx_turb_noHGTQS_noSHAL.nc')

heights = pd.read_csv(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\z_43_fullheights.csv')

#%%
'''Get the mean of the area'''
conv_dry_nohgtqs_mean = conv_dry_nohgtqs.mean(dim = ('x', 'y'))
conv_dry_nohgtqs_noshal_mean = conv_dry_nohgtqs_noshal.mean(dim = ('x', 'y'))
conv_moist_nohgtqs_mean = conv_moist_nohgtqs.mean(dim = ('x', 'y'))
conv_moist_nohgtqs_noshal_mean = conv_moist_nohgtqs_noshal.mean(dim = ('x', 'y'))
turb_nohgtqs_mean = turb_nohgtqs.mean(dim = ('x', 'y'))
turb_nohgtqs_noshal_mean = turb_nohgtqs_noshal.mean(dim = ('x', 'y'))


#%%
'''Tranpose the data'''
conv_dry_nohgtqs_mean_trans = conv_dry_nohgtqs_mean.transpose()
conv_dry_nohgtqs_noshal_trans = conv_dry_nohgtqs_noshal_mean.transpose()
conv_moist_nohgtqs_trans = conv_moist_nohgtqs_mean.transpose()
conv_moist_nohgtqs_noshal_trans = conv_moist_nohgtqs_noshal_mean.transpose()
turb_nohgtqs_trans = turb_nohgtqs_mean.transpose()
turb_nohgtqs_noshal_trans = turb_nohgtqs_noshal_mean.transpose()


#%%
'''Plot the data'''
plt.figure(figsize = (15,10))
plt.contourf(conv_dry_nohgtqs_mean_trans.time.values, heights['Full heights'], np.flip(conv_dry_nohgtqs_mean_trans.Thlflx_conv_dry.values, axis = 0))
plt.colorbar(label = 'theta_l Flux [K/ms]')
plt.title('Accumulated theta_l Flux Dry Convection for HARMONIE noHGTQS')
plt.ylabel('Height [m]')
plt.xlabel('Time')
plt.ylim([0,3000])
plt.show()

plt.figure(figsize = (15,10))
plt.contourf(conv_dry_nohgtqs_noshal_trans.time.values, heights['Full heights'], np.flip(conv_dry_nohgtqs_noshal_trans.Thlflx_conv_dry.values, axis = 0))
plt.colorbar(label = 'theta_l Flux [K/ms]')
plt.title('Accumulated theta_l Flux Dry Convection for HARMONIE noHGTQS noSHAL')
plt.ylabel('Height [m]')
plt.xlabel('Time')
plt.ylim([0,3000])
plt.show()


plt.figure(figsize = (15,10))
plt.contourf(conv_moist_nohgtqs_trans.time.values, heights['Full heights'], np.flip(conv_moist_nohgtqs_trans.Thlflx_conv_mois.values, axis = 0))
plt.colorbar(label = 'theta_l Flux [K/ms]')
plt.title('Accumulated theta_l Flux Moist Convection for HARMONIE noHGTQS')
plt.ylabel('Height [m]')
plt.xlabel('Time')
plt.ylim([0,3000])
plt.show()


plt.figure(figsize = (15,10))
plt.contourf(conv_moist_nohgtqs_noshal_trans.time.values, heights['Full heights'], np.flip(conv_moist_nohgtqs_noshal_trans.Thlflx_conv_mois.values, axis = 0))
plt.colorbar(label = 'theta_l Flux [K/ms]')
plt.title('Accumulated theta_l Flux moist Convection for HARMONIE noHGTQS noSHAL')
plt.ylabel('Height [m]')
plt.xlabel('Time')
plt.ylim([0,3000])
plt.show()


plt.figure(figsize = (15,10))
plt.contourf(turb_nohgtqs_trans.time.values, heights['Full heights'], np.flip(turb_nohgtqs_trans.Thlflx_turb.values, axis = 0))
plt.colorbar(label = 'theta_l Flux [K/ms]')
plt.title('Accumulated theta_l Flux turbulent for HARMONIE noHGTQS')
plt.ylabel('Height [m]')
plt.xlabel('Time')
plt.ylim([0,3000])
plt.show()


plt.figure(figsize = (15,10))
plt.contourf(turb_nohgtqs_noshal_trans.time.values, heights['Full heights'], np.flip(turb_nohgtqs_noshal_trans.Thlflx_turb.values, axis = 0))
plt.colorbar(label = 'theta_l Flux [K/ms]')
plt.title('Accumulated theta_l Flux turbulent HARMONIE noHGTQS noSHAL')
plt.ylabel('Height [m]')
plt.xlabel('Time')
plt.ylim([0,3000])
plt.show()


#%%
'''Remove the accumulation'''
conv_dry_nohgtqs_nonaccum = []

for i in range(len(conv_dry_nohgtqs_mean_trans.time.values)):
    for ii in range(len(conv_dry_nohgtqs_mean_trans.lev.values)):
        print(i)
        if i == 0 :
            conv_dry_nohgtqs_nonaccum.append(conv_moist_nohgtqs_trans.Thlflx_conv_mois.values[i,ii])
        
        
