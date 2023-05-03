# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 13:07:02 2023

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
from cut_subgrid_module import cut_sub_grid


evs_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\evs_data\evspsbl_noHGTQS_cut.nc', combine ='by_coords') 
evs_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\evs_data\evspsbl_noHGTQS_noSHAL_cut.nc', combine ='by_coords')
evs_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\evs_data\evspsbl_noHGTQS_noUVmix_cut.nc', combine ='by_coords')

hfls_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\hfls_data\hfls_eva_noHGTQS_cut.nc', combine ='by_coords') #wrong
hfls_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\hfls_data\hfls_eva_noHGTQS_noSHAL_cut.nc', combine ='by_coords')
hfls_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\hfls_data\hfls_eva_noHGTQS_noUVmix_cut.nc', combine ='by_coords')
       
# =============================================================================
# evs_noHGTQS = cut_sub_grid(evs_noHGTQS)
# hfls_noHGTQS = cut_sub_grid(hfls_noHGTQS)
# 
# evs_noHGTQS.to_netcdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\evspsbl_noHGTQS_cut')
# hfls_noHGTQS.to_netcdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\hfls_eva_noHGTQS_cut')
# =============================================================================

non_accum_evs_noHGTQS = np.load(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\evs_data\hr_evs_noHGTQS.npy')
non_accum_evs_noHGTQS_noSHAL = np.load(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\evs_data\hr_evs_noHGTQS_noSHAL.npy')
non_accum_evs_noHGTQS_noUVmix = np.load(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\evs_data\hr_evs_noHGTQS_noUVmix.npy')

non_accum_hfls_noHGTQS = np.load(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\hfls_data\hr_hfls_noHGTQS.npy')
non_accum_hfls_noHGTQS_noSHAL = np.load(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\hfls_data\hr_hfls_noHGTQS_noSHAL.npy')
non_accum_hfls_noHGTQS_noUVmix = np.load(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\hfls_data\hr_hfls_noHGTQS_noUVmix.npy')



plt.figure(figsize = (15,10))
plt.plot(evs_noHGTQS.time.values, non_accum_evs_noHGTQS, label = 'HARMONIE noHGTQS')
plt.plot(evs_noHGTQS_noSHAL.time.values, non_accum_evs_noHGTQS_noSHAL, label = 'HARMONIE noHGTQS noSHAL')
plt.plot(evs_noHGTQS_noUVmix.time.values, non_accum_evs_noHGTQS_noUVmix, label = 'HARMONIE noHGTQS noUVmix')
plt.grid()
plt.title('Surface Water Evaporation Flux')
plt.xlabel('Time [hr]')
plt.ylabel('Surface Water Evaporation Flux [kg m-2]')
plt.legend()

plt.figure(figsize = (15,10))
plt.plot(hfls_noHGTQS.time.values, non_accum_hfls_noHGTQS, label = 'HARMONIE noHGTQS')
plt.plot(hfls_noHGTQS_noSHAL.time.values, non_accum_hfls_noHGTQS_noSHAL, label = 'HARMONIE noHGTQS noSHAL')
plt.plot(hfls_noHGTQS_noUVmix.time.values, non_accum_hfls_noHGTQS_noUVmix, label = 'HARMONIE noHGTQS noUVmix')
plt.grid()
plt.title('Surface Upward Latent Heat Flux Due To Evavaporation')
plt.xlabel('Time [hr]')
plt.ylabel('Surface Upward Latent Heat Flux [J m-2]')
plt.legend()


xr_evs_noHGTQS = xr.Dataset({'evs': (('time'), non_accum_evs_noHGTQS )},
                  {'time': ( evs_noHGTQS.time.values)}
                  )

xr_evs_noHGTQS_noSHAL = xr.Dataset({'evs': (('time'), non_accum_evs_noHGTQS_noSHAL )},
                  {'time': ( evs_noHGTQS.time.values)}
                  )

xr_evs_noHGTQS_noUVmix = xr.Dataset({'evs': (('time'), non_accum_evs_noHGTQS_noUVmix )},
                  {'time': ( evs_noHGTQS.time.values)}
                  )



xr_hfls_noHGTQS = xr.Dataset({'hfls': (('time'), non_accum_hfls_noHGTQS )},
                  {'time': ( hfls_noHGTQS.time.values)}
                  )

xr_hfls_noHGTQS_noSHAL = xr.Dataset({'hfls': (('time'), non_accum_hfls_noHGTQS_noSHAL )},
                  {'time': ( hfls_noHGTQS_noSHAL.time.values)}
                  )

xr_hfls_noHGTQS_noUVmix = xr.Dataset({'hfls': (('time'), non_accum_hfls_noHGTQS_noUVmix )},
                  {'time': ( hfls_noHGTQS_noUVmix.time.values)}
                  )


evs_gr_noHGTQS = xr_evs_noHGTQS.groupby(xr_evs_noHGTQS.time.dt.hour).mean()
evs_gr_noHGTQS_noSHAL =xr_evs_noHGTQS_noSHAL.groupby(xr_evs_noHGTQS_noSHAL.time.dt.hour).mean()
evs_gr_noHGTQS_noUVmix = xr_evs_noHGTQS_noUVmix.groupby(xr_evs_noHGTQS_noUVmix.time.dt.hour).mean()

hfls_gr_noHGTQS = xr_hfls_noHGTQS.groupby(xr_hfls_noHGTQS.time.dt.hour).mean()
hfls_gr_noHGTQS_noSHAL =xr_hfls_noHGTQS_noSHAL.groupby(xr_hfls_noHGTQS_noSHAL.time.dt.hour).mean()
hfls_gr_noHGTQS_noUVmix = xr_hfls_noHGTQS_noUVmix.groupby(xr_hfls_noHGTQS_noUVmix.time.dt.hour).mean()


hours = np.linspace(0,23,24)

plt.figure(figsize = (15,10))
plt.plot(evs_gr_noHGTQS.hour, evs_gr_noHGTQS.evs.values, label = 'HARMONIE noHGTQS')
plt.plot(evs_gr_noHGTQS_noSHAL.hour, evs_gr_noHGTQS_noSHAL.evs.values, label = 'HARMONIE noHGTQS noSHAL')
plt.plot(evs_gr_noHGTQS_noUVmix.hour, evs_gr_noHGTQS_noUVmix.evs.values, label = 'HARMONIE noHGTQS noUVmix')
plt.title('Composite Diurnal Cycle of Surface Water Evaporation Flux')
plt.grid()
plt.xlabel('Time [hr]')
plt.ylabel('Surface Water Evaporation Flux [kg m-2]')
plt.legend()
plt.xticks(evs_gr_noHGTQS_noUVmix.hour, hours)


plt.figure(figsize = (15,10))
plt.plot(hfls_gr_noHGTQS.hour, hfls_gr_noHGTQS.hfls.values, label = 'HARMONIE noHGTQS')
plt.plot(hfls_gr_noHGTQS_noSHAL.hour, hfls_gr_noHGTQS_noSHAL.hfls.values, label = 'HARMONIE noHGTQS noSHAL')
plt.plot(hfls_gr_noHGTQS_noUVmix.hour, hfls_gr_noHGTQS_noUVmix.hfls.values, label = 'HARMONIE noHGTQS noUVmix')
plt.title('Composite Diurnal Cycle of Surface Upward Latent Heat Flux')
plt.grid()
plt.xlabel('Time [hr]')
plt.ylabel('Surface Upward Latent Heat Flux [J m-2 hr -1]')
plt.legend()
plt.xticks(evs_gr_noHGTQS_noUVmix.hour, hours)

