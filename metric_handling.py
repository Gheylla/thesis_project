# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 15:59:09 2023

@author: User
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
import datetime as dt

'''Import data'''
# =============================================================================
# data = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cut_sat\OR_ABI-L2-CMIPF-M6C13_G16_s20200340320124_e20200340329444_c20200340329541.nc')
# 
# data = data.isel(time = 0)
# plt.contourf(data.lon.values, data.lat.values, data.C13_Brightness_Temp.values)
# =============================================================================
'''Metrics without resize function'''
df_metrics_noHGTQS43 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\HARMONIE\df_metrics_noHGTQS.h5')
df_metrics_noHGTQS_noSHAL = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\HARMONIE\df_metrics_noHGTQS_noSHAL.h5')
df_metrics_noHGTQS_noUVmix = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\HARMONIE\df_metrics_noHGTQS_noUVmix.h5')

cll_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cll_data\cll_noHGTQS_cut.nc')
cll_noHGTQSnoSHAL =  xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cll_data\cll_noHGTQS_noSHAL_cut.nc')
cll_noHGTQSnoUVmix =  xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cll_data\cll_noHGTQS_noUVmix_cut.nc')


metrics_obs = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs.h5')
metrics_obs_1 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_1.h5')
metrics_obs_2 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_2.h5')
metrics_obs_3 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_3.h5')
metrics_obs_4 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_4.h5')
metrics_obs_5 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_5.h5')
metrics_obs_6 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_6.h5')
metrics_obs_7 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_7.h5')
metrics_obs_8 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_8.h5')
metrics_obs_9 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_9.h5')
metrics_obs_10 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_10.h5')
metrics_obs_11 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_11.h5')
metrics_obs_12 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_12.h5')
metrics_obs_13 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_13.h5')
metrics_obs_14 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_14.h5')
metrics_obs_15 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_15.h5')
#metrics_obs_16 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_16.h5')
#number 16 is wrong it was done using the wrong threshold
metrics_obs_17 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_17.h5')
metrics_obs_18 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_18.h5')
metrics_obs_19 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_19.h5')
metrics_obs_20 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\observations\metrics_obs_20.h5')



plt.figure(figsize = (20,10))

'''Plotting observations'''
plt.plot(metrics_obs['Time'], metrics_obs['mean_length_scale'], color = 'black', linewidth=2.5, label = 'Observations')
plt.plot(metrics_obs_1['Time'], metrics_obs_1['mean_length_scale'], color = 'black', linewidth=2.5)
plt.plot(metrics_obs_2['Time'], metrics_obs_2['mean_length_scale'], color = 'black', linewidth=2.5)
plt.plot(metrics_obs_3['Time'], metrics_obs_3['mean_length_scale'], color = 'black', linewidth=2.5)
plt.plot(metrics_obs_4['Time'], metrics_obs_4['mean_length_scale'], color = 'black', linewidth=2.5)
plt.plot(metrics_obs_5['Time'], metrics_obs_5['mean_length_scale'], color = 'black', linewidth=2.5)
plt.plot(metrics_obs_6['Time'], metrics_obs_6['mean_length_scale'], color = 'black', linewidth=2.5)
plt.plot(metrics_obs_7['Time'], metrics_obs_7['mean_length_scale'], color = 'black', linewidth=2.5)
plt.plot(metrics_obs_8['Time'], metrics_obs_8['mean_length_scale'], color = 'black', linewidth=2.5)
plt.plot(metrics_obs_9['Time'], metrics_obs_9['mean_length_scale'], color = 'black', linewidth=2.5)
plt.plot(metrics_obs_10['Time'], metrics_obs_10['mean_length_scale'], color = 'black', linewidth=2.5)
plt.plot(metrics_obs_11['Time'], metrics_obs_11['mean_length_scale'], color = 'black', linewidth=2.5)
plt.plot(metrics_obs_12['Time'], metrics_obs_12['mean_length_scale'], color = 'black', linewidth=2.5)
plt.plot(metrics_obs_13['Time'], metrics_obs_13['mean_length_scale'], color = 'black', linewidth=2.5)
plt.plot(metrics_obs_14['Time'], metrics_obs_14['mean_length_scale'], color = 'black', linewidth=2.5)
plt.plot(metrics_obs_15['Time'], metrics_obs_15['mean_length_scale'], color = 'black', linewidth=2.5)
#plt.plot(metrics_obs_16['Time'], metrics_obs_16['mean_length_scale'], color = 'black', linewidth=2.5)
plt.plot(metrics_obs_17['Time'], metrics_obs_17['mean_length_scale'], color = 'black', linewidth=2.5)
#plt.plot(metrics_obs_18['Time'][1:], metrics_obs_18['mean_length_scale'][1:], color = 'black', linewidth=2.5)
plt.plot(metrics_obs_19['Time'], metrics_obs_19['mean_length_scale'], color = 'black', linewidth=2.5)
plt.plot(metrics_obs_20['Time'], metrics_obs_20['mean_length_scale'], color = 'black', linewidth=2.5)



'''Plotting experiments'''
plt.plot(df_metrics_noHGTQS43.index, df_metrics_noHGTQS43['mean_length_scale'], label = 'HARMONIE noHGTQS')
plt.plot(df_metrics_noHGTQS_noSHAL.index, df_metrics_noHGTQS_noSHAL['mean_length_scale'], label = 'HARMONIE noHGTQS noSHAL')
plt.plot(df_metrics_noHGTQS_noUVmix.index, df_metrics_noHGTQS_noUVmix['mean_length_scale'], label = 'HARMONIE noHGTQS noUVmix')


plt.title('Time Series of mean length scale')
plt.ylabel('Mean length scale[-]')
plt.xlabel('Time')
plt.legend()
plt.grid()


# =============================================================================
# '
# '''Plotting observations'''
# plt.plot(metrics_obs['Time'], metrics_obs['num_objects'], color = 'black', linewidth=2.5, linestyle = 'dashed', label = 'Observations')
# plt.plot(metrics_obs_1['Time'], metrics_obs_1['num_objects'], color = 'black', linewidth=2.5, linestyle = 'dashed')
# plt.plot(metrics_obs_2['Time'], metrics_obs_2['num_objects'], color = 'black', linewidth=2.5, linestyle = 'dashed')
# plt.plot(metrics_obs_3['Time'], metrics_obs_3['num_objects'], color = 'black', linewidth=2.5, linestyle = 'dashed')
# plt.plot(metrics_obs_4['Time'], metrics_obs_4['num_objects'], color = 'black', linewidth=2.5, linestyle = 'dashed')
# plt.plot(metrics_obs_5['Time'], metrics_obs_5['num_objects'], color = 'black', linewidth=2.5, linestyle = 'dashed')
# plt.plot(metrics_obs_6['Time'], metrics_obs_6['num_objects'], color = 'black', linewidth=2.5, linestyle = 'dashed')
# plt.plot(metrics_obs_7['Time'], metrics_obs_7['num_objects'], color = 'black', linewidth=2.5, linestyle = 'dashed')'
# 
# =============================================================================


# =============================================================================
# '''Plotting experiments'''
# plt.plot(df_metrics_noHGTQS43.index, df_metrics_noHGTQS43['num_objects'], label = 'HARMONIE noHGTQS')
# plt.plot(df_metrics_noHGTQS_noSHAL.index, df_metrics_noHGTQS_noSHAL['num_objects'], label = 'HARMONIE noHGTQS noSHAL')
# plt.plot(df_metrics_noHGTQS_noUVmix.index, df_metrics_noHGTQS_noUVmix['num_objects'], label = 'HARMONIE noHGTQS noUVmix')
# 
# =============================================================================



plt.title('Number of objects')
plt.ylabel('Number of objects[-]')
plt.xlabel('Time')
plt.legend()
plt.grid()


# =============================================================================
# 
# '''Compute the diurnal cycle of the low cloud cover'''
# cll_noHGTQS = cll_noHGTQS.mean(dim= ('x', 'y'))
# cll_noHGTQS = cll_noHGTQS.groupby(cll_noHGTQS.time.dt.hour).mean()
# 
# cll_noHGTQSnoSHAL = cll_noHGTQSnoSHAL.mean(dim= ('x', 'y'))
# cll_noHGTQSnoSHAL = cll_noHGTQSnoSHAL.groupby(cll_noHGTQSnoSHAL.time.dt.hour).mean()
# 
# cll_noHGTQSnoUVmix = cll_noHGTQSnoUVmix.mean(dim= ('x', 'y'))
# cll_noHGTQSnoUVmix = cll_noHGTQSnoUVmix.groupby(cll_noHGTQSnoUVmix.time.dt.hour).mean()
# =============================================================================

# =============================================================================
# 
# '''Compute the diurnal cycle of the low cloud cover'''
# open_noHGTQS = df_metrics_noHGTQS43['num_objects'].groupby(df_metrics_noHGTQS43.index.hour).mean()
# open_noHGTQS_noSHAL = df_metrics_noHGTQS_noSHAL['num_objects'].groupby(df_metrics_noHGTQS_noSHAL.index.hour).mean()
# open_noHGTQS_noUVmix = df_metrics_noHGTQS_noUVmix['num_objects'].groupby(df_metrics_noHGTQS_noUVmix.index.hour).mean()
# 
# 
# plt.figure(figsize= (15,10))
# plt.plot(open_noHGTQS.index, open_noHGTQS, label = 'HARMONIE noHGTQS')
# plt.plot(open_noHGTQS_noSHAL.index, open_noHGTQS_noSHAL, label = 'HARMONIE noHGTQS noSHAL')
# plt.plot(open_noHGTQS_noUVmix.index, open_noHGTQS_noUVmix, label = 'HARMONIE noHGTQS noUVmix')
# plt.title('Composite diurnal cycle of number of objects')
# plt.grid()
# plt.xlabel('Time [hr]')
# plt.ylabel('Number of objects')
# plt.legend()
# 
# 
# hours = np.linspace(0, 23, 24)
# plt.xticks(open_noHGTQS_noUVmix.index, hours)
# '''Plot the graphs'''
# =============================================================================
# =============================================================================
# plt.figure(figsize = (15,5))
# plt.plot(cll_noHGTQS.hour.values, cll_noHGTQS.cll.values, label = 'HARMONIE noHGTQS')
# plt.plot(cll_noHGTQSnoSHAL.hour.values, cll_noHGTQSnoSHAL.cll.values, label = 'HARMONIE noHGTQS noSHAL')
# plt.plot(cll_noHGTQSnoUVmix.hour.values, cll_noHGTQSnoUVmix.cll.values, label = 'HARMONIE noHGTQS noUVmix')
# plt.grid()
# plt.title('Composite diurnal cycle of low cloud cover')
# plt.xlabel('Time [hour]')
# plt.ylabel('Low cloud cover [-]')
# plt.legend()
# 
# hours = np.linspace(0,23,24)
# plt.xticks(cll_noHGTQSnoUVmix.hour.values, hours)
# =============================================================================

# =============================================================================
# numobj_noHGTQS = df_metrics_noHGTQS43.num_objects.groupby(df_metrics_noHGTQS43.index.hour).mean()
# numobj_noHGTQS_noSHAL = df_metrics_noHGTQS_noSHAL.num_objects.groupby(df_metrics_noHGTQS_noSHAL.index.hour).mean()
# numobj_noHGTQS_noUVmix = df_metrics_noHGTQS_noUVmix.num_objects.groupby(df_metrics_noHGTQS_noUVmix.index.hour).mean()
# 
# 
# hours = np.linspace(0,23,24)
# plt.figure(figsize = (15,5))
# plt.plot(hours, numobj_noHGTQS, label = 'HARMONIE cy43 no HGTQS')
# plt.plot(hours, numobj_noHGTQS_noSHAL, label = 'HARMONIE cy43 noHGTQS noSHAL', color = 'orange')
# plt.plot(hours, numobj_noHGTQS_noUVmix, label = 'HARMONIE cy43 noHGTQS noUVmix', color = 'green')
# plt.legend()
# plt.xlabel('Hour')
# plt.ylabel('Number of objects [-]')
# plt.title('Composite diurnal cycle of number of objects')
# plt.grid()
# plt.xticks(hours, hours)
# =============================================================================


# =============================================================================
# 
# #df_sat = pd.read_csv('C:/Users/LENOVO/Desktop/TU_Delft/thesis/csv_combined.csv', parse_dates= True , index_col=0)
# #df_sat.resample('H').mean()
# metrics = df_metrics_noHGTQS43.columns
# 
# 
# num_objects_noHGTQS = df_metrics_noHGTQS43.num_objects.groupby(df_metrics_noHGTQS43.index.hour).mean()
# num_objects_noHGTQS_noSHAL = df_metrics_noHGTQS_noSHAL.num_objects.groupby(df_metrics_noHGTQS_noSHAL.index.hour).mean()
# num_objects_noHGTQS_noUVmix = df_metrics_noHGTQS_noUVmix.num_objects.groupby(df_metrics_noHGTQS_noUVmix.index.hour).mean()
# 
# hours = np.linspace(0,23,24)
# =============================================================================

# =============================================================================
# plt.figure(figsize = (15,5))
# plt.plot(hours, num_objects_noHGTQS, label = 'HARMONIE cy43 no HGTQS')
# plt.plot(hours, num_objects_noHGTQS_noSHAL, label = 'HARMONIE cy43 noHGTQS noSHAL', color = 'orange')
# plt.plot(hours, num_objects_noHGTQS_noUVmix, label = 'HARMONIE cy43 noHGTQS noUVmix', color = 'green')
# plt.legend()
# plt.xlabel('Hour')
# plt.ylabel('I_org')
# plt.title('Composite diurnal cycle of I_org')
# plt.grid()
# plt.xticks(hours, hours)
# 
# =============================================================================

# =============================================================================
# for i in range(len(metrics)):
#     print(metrics[i])
#     '''Plot all of the metrics'''
#     plt.figure(figsize = (25,10))
#     plt.plot(df_metrics_noHGTQS43.index, df_metrics_noHGTQS43[metrics[i]], label = 'HARMONIE cy43 no HGTQS')
#     plt.plot(df_metrics_noHGTQS_noSHAL.index, df_metrics_noHGTQS_noSHAL[metrics[i]], label = 'HARMONIE cy43 noHGTQS noSHAL', color = 'orange')
#     plt.plot(df_metrics_noHGTQS_noUVmix.index, df_metrics_noHGTQS_noUVmix[metrics[i]], label = 'HARMONIE cy43 noHGTQS noUVmix', color = 'green')
#     
#     plt.title(f'Timeseries {metrics[i]}')
#     plt.grid()
#     plt.legend()
#     plt.xlabel('Time')
#     plt.ylabel(metrics[i])
#     
# =============================================================================

# =============================================================================
# plt.figure(figsize = (15,5))    
# plt.plot(df_metrics_noHGTQS43_noresize['num_objects'])
# plt.plot(df_metrics_noHGTQS_noSHAL_noresize['num_objects'])
# 
# plt.figure()
# plt.plot(df_sat.index, df_sat['num_objects'])
# 
# 
# =============================================================================
















