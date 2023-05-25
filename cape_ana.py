# -*- coding: utf-8 -*-
"""
Created on Wed May  3 11:44:40 2023

@author: LENOVO
"""


import eurec4a
from cut_subgrid_module import cut_sub_grid
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import cloudmetrics
import inspect
from skimage.transform import resize
import seaborn as sn


'''Import data'''
#CAPE is not dependend on level so is it on the surface? 
cape_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cape_data\cape_noHGTQS.nc')
cape_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cape_data\cape_noHGTQS_noSHAL.nc')
cape_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cape_data\cape_noHGTQS_noUVmix.nc')

# %%
'''Cut data for desired grid'''
cape_noHGTQS = cut_sub_grid(cape_noHGTQS)
cape_noHGTQS_noSHAL = cut_sub_grid(cape_noHGTQS_noSHAL)
cape_noHGTQS_noUVmix = cut_sub_grid(cape_noHGTQS_noUVmix)

#%%
'''Testing'''
levels = np.linspace(0, 3000, 10)

for i in range(len(cape_noHGTQS.time.values)):
    cape_noHGTQS_sel = cape_noHGTQS.isel(time = i)
    plt.contourf(cape_noHGTQS_sel.lon.values, cape_noHGTQS_sel.lat.values, cape_noHGTQS_sel.cape.values, vmin = 0, vmax = 3000, levels = levels)
    plt.title(f'CAPE for HARMONIE noHGTQS at time {cape_noHGTQS_sel.time.values}')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.colorbar(label = 'CAPE [J kg-1]')
    plt.show()
    
    
#%%    
for i in range(len(cape_noHGTQS_noSHAL.time.values)):
    cape_noHGTQS_noSHAL_sel = cape_noHGTQS_noSHAL.isel(time = i)
    plt.contourf(cape_noHGTQS_noSHAL_sel.lon.values, cape_noHGTQS_noSHAL_sel.lat.values, cape_noHGTQS_noSHAL_sel.cape.values, vmin = 0, vmax = 3000, levels = levels)
    plt.title(f'CAPE for HARMONIE noHGTQS noSHAL at time {cape_noHGTQS_noSHAL_sel.time.values}')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.colorbar(label = 'CAPE [J kg-1]')
    plt.show()
    
    
    
#%%
for i in range(len(cape_noHGTQS_noUVmix.time.values)):
    cape_noHGTQS_noUVmix_sel = cape_noHGTQS_noUVmix.isel(time = i)
    plt.contourf(cape_noHGTQS_noUVmix_sel.lon.values, cape_noHGTQS_noUVmix_sel.lat.values, cape_noHGTQS_noUVmix_sel.cape.values, vmin = 0, vmax = 3000, levels = levels)
    plt.title(f'CAPE for HARMONIE noHGTQS noUVmix at time {cape_noHGTQS_noUVmix_sel.time.values}')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.colorbar(label = 'CAPE [J kg-1]')
    plt.show()    

# %%
'''Plot the surface cape'''
plt.figure()
sn.distplot(cape_noHGTQS.cape.values, label = 'noHGTQS', hist = False)
sn.distplot(cape_noHGTQS_noSHAL.cape.values, label = 'noHGTQS noSHAL', hist = False)
sn.distplot(cape_noHGTQS_noUVmix.cape.values, label = 'noHGTQS noUVmix', hist = False)
plt.legend()
plt.title("Density plot of CAPE")
plt.xlabel('Count [-]')
plt.grid()