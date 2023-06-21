# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 15:30:34 2023

@author: Gheylla Liberia
Make snapshot of resolved mass flux
by getting all the areas with cloud cover >0.5 and multiplying this with the vertical windspeed.
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
from cut_subgrid_module import cut_sub_grid
#%%
'''Import data'''
w_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\w_data\wa_noHTQS_total.nc')
w_noHGTQS_nouvmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\w_data\wa_noHTQS_noSHAL_total.nc')
w_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\w_data\wa_noHTQS_noUVmix_total.nc')

heights = pd.read_csv(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\z_43_fullheights.csv')

cf_nohgtqs = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cf_data\cl_Slev_noHGTQS.nc')
cf_nohgtqs_noshal = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cf_data\cl_Slev_noHGTQS_noSHAL.nc')
cf_nohgtqs_nouvmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cf_data\cl_Slev_noHGTQS_noUVmix.nc')

#%%
'''Get the mask'''

cll_threshold = 0.5
cloud_mask_nohgtqs = np.zeros(cf_nohgtqs.cl.shape,dtype=int)
cloud_mask_nohgtqs_noshal = np.zeros(cf_nohgtqs_noshal.cl.shape,dtype=int)
cloud_mask_nohgtqs_nouvmix = np.zeros(cf_nohgtqs_nouvmix.cl.shape,dtype=int)



cloud_mask_nohgtqs[(cf_nohgtqs.cl > cll_threshold)] = 1
cloud_mask_nohgtqs_noshal[(cf_nohgtqs_noshal.cl > cll_threshold)] = 1
cloud_mask_nohgtqs_nouvmix[(cf_nohgtqs_nouvmix.cl > cll_threshold)] = 1


#%%
'''Save in data array to make easier'''
xr_mask_nohgtqs = xr.Dataset({'cll_mask': (('time', 'lev', 'y', 'x'), cloud_mask_nohgtqs )},
                  {'time': (cf_nohgtqs.time.values)}, 
                  {'units': ('-'), 
                   'long_name': ('cloud mask for HARMONIE noHGTQS')})

xr_mask_nohgtqs_noshal = xr.Dataset({'cll_mask': (('time', 'lev', 'y', 'x'), cloud_mask_nohgtqs_noshal )},
                  {'time': (cf_nohgtqs_noshal.time.values)}, 
                  {'units': ('-'), 
                   'long_name': ('cloud mask for HARMONIE noHGTQS noSHAL')})


xr_mask_nohgtqs_nouvmix = xr.Dataset({'cll_mask': (('time', 'lev', 'y', 'x'), cloud_mask_nohgtqs_nouvmix )},
                  {'time': (cf_nohgtqs_nouvmix.time.values)}, 
                  {'units': ('-'), 
                   'long_name': ('cloud mask for HARMONIE noHGTQS noUVmix')})

xr_mask_nohgtqs.to_netcdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\mask_cf_nohgtqs.nc')
xr_mask_nohgtqs_noshal.to_netcdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\mask_cf_nohgtqs_noshal.nc')
xr_mask_nohgtqs_nouvmix.to_netcdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\mask_cf_nohgtqs_nouvmix.nc')

#%%
'''Import data again but the mask immediately'''
w_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\w_data\wa_noHTQS_total.nc')
w_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\w_data\wa_noHTQS_noSHAL_total.nc')
w_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\w_data\wa_noHTQS_noUVmix_total.nc')

heights = pd.read_csv(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\z_43_fullheights.csv')

mask_nohgtqs = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\mask_cf_nohgtqs.nc')
mask_nohgtqs_noshal = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\mask_cf_nohgtqs_noshal.nc')
mask_nohgtqs_nouvmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\mask_cf_nohgtqs_nouvmix.nc')


#%%
'''Calculate the mask flux'''
w_noHGTQS_sel = w_noHGTQS.isel(lev = 36)
mask_nohgtqs_sel = mask_nohgtqs.isel(lev = 36)
levels = np.linspace(-0.5, 2 , 11)

for i in range(len(w_noHGTQS_sel.time.values)):
    print(i)
    w_noHGTQS_sel_time = w_noHGTQS_sel.isel(time = i)
    mask_nohgtqs_sel_time = mask_nohgtqs_sel.isel(time = i)
    plt.figure(figsize = (15,15))
    mass_flux = w_noHGTQS_sel_time.wa.values * mask_nohgtqs_sel_time.cll_mask.values
    plt.contourf(w_noHGTQS_sel.lon.values, w_noHGTQS_sel.lat.values, mass_flux, vmin = -0.5, vmax = 1.8, levels = levels)
    plt.title(f'HARMONIE noHGTQS at height 2413m at time {w_noHGTQS_sel.time.values[i]}')
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    plt.colorbar(label = 'Mass flux [m/s]')
    plt.show()
    
    
#%%
'''For noshal'''
w_noHGTQS_noSHAL_sel = w_noHGTQS_noSHAL.isel(lev = 36)
mask_nohgtqs_noshal_sel = mask_nohgtqs_noshal.isel(lev = 36)
levels = np.linspace(-0.5, 2 , 11)

for i in range(len(w_noHGTQS_noSHAL_sel.time.values)):
    print(i)
    w_noHGTQS_noSHAL_sel_time = w_noHGTQS_noSHAL_sel.isel(time = i)
    mask_nohgtqs_noshal_sel_time = mask_nohgtqs_noshal_sel.isel(time = i)
    plt.figure(figsize = (15,15))
    mass_flux_noshal = w_noHGTQS_noSHAL_sel_time.wa.values * mask_nohgtqs_noshal_sel_time.cll_mask.values
    plt.contourf(w_noHGTQS_noSHAL_sel.lon.values, w_noHGTQS_noSHAL_sel.lat.values, mass_flux_noshal, vmin = -0.5, vmax = 1.8, levels = levels)
    plt.title(f'HARMONIE noHGTQS noSHAL at height 2413m at time {w_noHGTQS_noSHAL_sel.time.values[i]}')
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    plt.colorbar(label = 'Mass flux [m/s]')
    plt.show()




#%%
'''For no uvmix'''
w_noHGTQS_noUVmix_sel = w_noHGTQS_noUVmix.isel(lev = 36)
mask_nohgtqs_nouvmix_sel = mask_nohgtqs_nouvmix.isel(lev = 36)
levels = np.linspace(-0.5, 2 , 11)

for i in range(len(w_noHGTQS_noUVmix_sel.time.values)):
    print(i)
    w_noHGTQS_noUVmix_sel_time = w_noHGTQS_noUVmix_sel.isel(time = i)
    mask_nohgtqs_nouvmix_sel_time = mask_nohgtqs_nouvmix_sel.isel(time = i)
    plt.figure(figsize = (15,15))
    mass_flux_nouvmix = w_noHGTQS_noUVmix_sel_time.wa.values * mask_nohgtqs_nouvmix_sel_time.cll_mask.values
    plt.contourf(w_noHGTQS_noUVmix_sel.lon.values, w_noHGTQS_noUVmix_sel.lat.values, mass_flux_nouvmix, vmin = -0.5, vmax = 1.8, levels = levels)
    plt.title(f'HARMONIE noHGTQS noUVmix at height 2413m at time {w_noHGTQS_noUVmix_sel.time.values[i]}')
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    plt.colorbar(label = 'Mass flux [m/s]')
    plt.show()