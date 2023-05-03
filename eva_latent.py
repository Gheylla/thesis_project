# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 09:24:09 2023

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

hfls_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\hfls_data\hfls_eva_noHGTQS_cut.nc', combine ='by_coords')
hfls_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\hfls_data\hfls_eva_noHGTQS_noSHAL_cut.nc', combine ='by_coords')
hfls_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\hfls_data\hfls_eva_noHGTQS_noUVmix_cut.nc', combine ='by_coords')



evs_noHGTQS = evs_noHGTQS.evspsbl.mean(dim = ('x', 'y'))
evs_noHGTQS_noSHAL = evs_noHGTQS_noSHAL.evspsbl.mean(dim = ('x', 'y'))
evs_noHGTQS_noUVmix = evs_noHGTQS_noUVmix.evspsbl.mean(dim = ('x', 'y'))

min_evs_noHGTQS = np.min(evs_noHGTQS.values)
min_evs_noHGTQS_noSHAL = np.min(evs_noHGTQS_noSHAL.values)
min_evs_noHGTQS_noUVmix = np.min(evs_noHGTQS_noUVmix.values)

loc_min_noHGTQS = np.where(evs_noHGTQS.values == min_evs_noHGTQS )
loc_min_noHGTQS_noSHAL = np.where(evs_noHGTQS_noSHAL.values == min_evs_noHGTQS_noSHAL )
loc_min_noHGTQS_noSHAL = np.where(evs_noHGTQS_noUVmix.values == min_evs_noHGTQS_noUVmix )

hfls_noHGTQS = hfls_noHGTQS.mean(dim = ('x', 'y'))
hfls_noHGTQS_noSHAL = hfls_noHGTQS_noSHAL.mean(dim = ('x', 'y'))
hfls_noHGTQS_noUVmix = hfls_noHGTQS_noUVmix.mean(dim = ('x', 'y'))

plt.figure(figsize = (15,10))
plt.plot(evs_noHGTQS.time.values, evs_noHGTQS.values, label = 'HARMONIE noHGTQS')
plt.plot(evs_noHGTQS_noSHAL.time.values, evs_noHGTQS_noSHAL.values, label = 'HARMONIE noHGTQS noSHAL')
plt.plot(evs_noHGTQS_noUVmix.time.values, evs_noHGTQS_noUVmix.values, label = 'HARMONIE noHGTQS noUVmix')
plt.grid()
plt.title('Accumulated Surface Water Evaporation Flux')
plt.xlabel('Time [hr]')
plt.ylabel('Surface Water Evaporation Flux [kg m-2]')
plt.legend()

plt.figure(figsize = (15,10))
plt.plot(hfls_noHGTQS.time.values, hfls_noHGTQS.hfls_eva.values, label = 'HARMONIE noHGTQS' )
plt.plot(hfls_noHGTQS_noSHAL.time.values, hfls_noHGTQS_noSHAL.hfls_eva.values, label = 'HARMONIE noHGTQS noSHAL')
plt.plot(hfls_noHGTQS_noUVmix.time.values, hfls_noHGTQS_noUVmix.hfls_eva.values, label = 'HARMONIE noHGTQS noUVmix')
plt.grid()
plt.title('Accumulated Surface Upward Latent Heat Flux Due To Evavaporation')
plt.xlabel('Time [hr]')
plt.ylabel('Surface Upward Latent Heat Flux [J m-2]')
plt.legend()


hr_evs_noHGTQS = np.zeros(len(evs_noHGTQS.time.values))
hr_evs_noHGTQS_noSHAL = np.zeros(len(evs_noHGTQS_noSHAL.time.values))
hr_evs_noHGTQS_noUVmix = np.zeros(len(evs_noHGTQS_noUVmix.time.values))

for i in range((len(evs_noHGTQS.time.values)) - 1):
    print(i)
    if i == 0: 
        hr_evs_noHGTQS[i] = evs_noHGTQS.values[i]
        hr_evs_noHGTQS_noSHAL[i] = evs_noHGTQS_noSHAL.values[i]
        hr_evs_noHGTQS_noUVmix[i] = evs_noHGTQS_noUVmix.values[i]
    else:
        evs_noHGTQS_before = evs_noHGTQS.values[i - 1]
        evs_noHGTQS_noSHAL_before = evs_noHGTQS_noSHAL.values[i - 1]
        evs_noHGTQS_noUVmix_before = evs_noHGTQS_noUVmix.values[i - 1]
        evs_noHGTQS_now = evs_noHGTQS.values[i]
        evs_noHGTQS_noSHAL_now = evs_noHGTQS_noSHAL.values[i]
        evs_noHGTQS_noUVmix_now = evs_noHGTQS_noUVmix.values[i]
        hr_evs_noHGTQS[i] =  evs_noHGTQS_now -  evs_noHGTQS_before
        hr_evs_noHGTQS_noSHAL[i] =  evs_noHGTQS_noSHAL_now -  evs_noHGTQS_noSHAL_before
        hr_evs_noHGTQS_noUVmix[i] =  evs_noHGTQS_noUVmix_now -  evs_noHGTQS_noUVmix_before
        if i == loc_min_noHGTQS[0][0]:
            hr_evs_noHGTQS[i] = evs_noHGTQS.values[i]
            hr_evs_noHGTQS_noSHAL[i] = evs_noHGTQS_noSHAL.values[i]
            hr_evs_noHGTQS_noUVmix[i] = evs_noHGTQS_noUVmix.values[i]


plt.figure(figsize = (15,10))
plt.plot(evs_noHGTQS.time.values, hr_evs_noHGTQS, label = 'HARMONIE noHGTQS ' )          
plt.plot(evs_noHGTQS.time.values, hr_evs_noHGTQS_noSHAL, label = 'HARMONIE noHGTQS noSHAL' )   
plt.plot(evs_noHGTQS.time.values, hr_evs_noHGTQS_noUVmix, label = 'HARMONIE noHGTQS noUVmix')     
plt.grid()
plt.title('Surface Water Evaporation Flux')
plt.xlabel('Time [hr]')
plt.ylabel('Surface Water Evaporation Flux [kg m-2]')
plt.legend()



hr_hfls_noHGTQS = np.zeros(len(hfls_noHGTQS.time.values))
hr_hfls_noHGTQS_noSHAL = np.zeros(len(hfls_noHGTQS_noSHAL.time.values))
hr_hfls_noHGTQS_noUVmix = np.zeros(len(hfls_noHGTQS_noUVmix.time.values))

for i in range((len(hfls_noHGTQS.time.values)) - 1):
    print(i)
    if i == 0: 
        hr_hfls_noHGTQS[i] = hfls_noHGTQS.hfls_eva.values[i]
        hr_hfls_noHGTQS_noSHAL[i] = hfls_noHGTQS_noSHAL.hfls_eva.values[i]
        hr_hfls_noHGTQS_noUVmix[i] = hfls_noHGTQS_noUVmix.hfls_eva.values[i]
    else:
        hfls_noHGTQS_before = hfls_noHGTQS.hfls_eva.values[i - 1]
        hfls_noHGTQS_noSHAL_before = hfls_noHGTQS_noSHAL.hfls_eva.values[i - 1]
        hfls_noHGTQS_noUVmix_before = hfls_noHGTQS_noUVmix.hfls_eva.values[i - 1]
        hfls_noHGTQS_now = hfls_noHGTQS.hfls_eva.values[i]
        hfls_noHGTQS_noSHAL_now = hfls_noHGTQS_noSHAL.hfls_eva.values[i]
        hfls_noHGTQS_noUVmix_now = hfls_noHGTQS_noUVmix.hfls_eva.values[i]
        hr_hfls_noHGTQS[i] =  hfls_noHGTQS_now -  hfls_noHGTQS_before
        hr_hfls_noHGTQS_noSHAL[i] = hfls_noHGTQS_noSHAL_now - hfls_noHGTQS_noSHAL_before
        hr_hfls_noHGTQS_noUVmix[i] = hfls_noHGTQS_noUVmix_now - hfls_noHGTQS_noUVmix_before
        if i == loc_min_noHGTQS[0][0]:
            hr_hfls_noHGTQS[i] = hfls_noHGTQS.hfls_eva.values[i]
            hr_hfls_noHGTQS_noSHAL[i] = hfls_noHGTQS_noSHAL.hfls_eva.values[i]
            hr_hfls_noHGTQS_noUVmix[i] = hfls_noHGTQS_noUVmix.hfls_eva.values[i]
  
            
            
plt.figure(figsize = (15,10))
plt.plot(hfls_noHGTQS.time.values, hr_hfls_noHGTQS, label = 'HARMONIE noHGTQS')          
plt.plot(hfls_noHGTQS_noSHAL.time.values, hr_hfls_noHGTQS_noSHAL, label = 'HARMONIE noHGTQS noSHAL')   
plt.plot(hfls_noHGTQS_noUVmix.time.values, hr_hfls_noHGTQS_noUVmix, label = 'HARMONIE noHGTQS noUVmix')     
plt.grid()
plt.title('Surface Upward Latent Heat Flux Due To Evavaporation')
plt.xlabel('Time [hr]')
plt.ylabel('Surface Upward Latent Heat Flux [J m-2]')
plt.legend()
