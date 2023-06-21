# -*- coding: utf-8 -*-
"""
Created on Wed May  3 11:44:40 2023

@author: LENOVO
"""


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
'''Get the mean of the area'''
cape_noHGTQS_mean = cape_noHGTQS.mean(dim = ('x', 'y'))
cape_noHGTQS_noSHAL_mean = cape_noHGTQS_noSHAL.mean(dim = ('x', 'y'))
cape_noHGTQS_noUVmix_mean = cape_noHGTQS_noUVmix.mean(dim = ('x', 'y'))


#%%
'''Determine the diurnal cycle of cape'''
cape_noHGTQS_diur = cape_noHGTQS_mean.groupby(cape_noHGTQS_mean.time.dt.hour).mean()
cape_noHGTQS_noSHAL_diur = cape_noHGTQS_noSHAL_mean.groupby(cape_noHGTQS_noSHAL_mean.time.dt.hour).mean()
cape_noHGTQS_noUVmix_diur = cape_noHGTQS_noUVmix_mean.groupby(cape_noHGTQS_noUVmix_mean.time.dt.hour).mean()


#%%
'''Transpose the data for plotting'''
cape_noHGTQS_diur_trans  = cape_noHGTQS_diur.transpose()
cape_noHGTQS_noSHAL_diur_trans = cape_noHGTQS_noSHAL_diur.transpose()
cape_noHGTQS_noUVmix_diur_trans = cape_noHGTQS_noUVmix_diur.transpose()


#%%
'''Get contour plot of diurnalcycle'''
hours = np.linspace(00,23, 24)
lt = [20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

  
fig, ax = plt.subplots(layout='constrained', figsize = (15,5))
ax.plot(cape_noHGTQS_diur_trans.hour.values, cape_noHGTQS_diur_trans.cape.values, label = 'HARMONIE noHGQTS', color = 'red')
ax.plot(cape_noHGTQS_diur_trans.hour.values, cape_noHGTQS_noSHAL_diur_trans.cape.values, label = 'HARMONIE noHGQTS noSHAL', color = 'blue')
ax.plot(cape_noHGTQS_diur_trans.hour.values, cape_noHGTQS_noUVmix_diur_trans.cape.values, label = 'HARMONIE noHGQTS noUVmix', color = 'green')



ax.grid()
ax.set_xticks(hours)
ax.set_xlabel('UTC Time [hr]')
secax = ax.secondary_xaxis('top')
secax.set_xticks(hours, lt)
secax.set_xlabel('Local Time [hr]')
plt.ylabel('CAPE [J/kg]')
plt.title('Composite diurnal cycle of CAPE')
plt.legend()
plt.show()
    
#%%
'''plot cape for video'''
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
'''Plot the surface cape distribution'''
plt.figure()
sn.distplot(cape_noHGTQS.cape.values,  hist = False, color = 'red', label = 'HARMONIE noHGTQS')
sn.distplot(cape_noHGTQS_noSHAL.cape.values,  hist = False, color = 'blue', label = 'HARMONIE noHGTQS noSHAL')
sn.distplot(cape_noHGTQS_noUVmix.cape.values,  hist = False, color = 'green', label = 'HARMONIE noHGTQS noUVmix')
plt.legend()
plt.title("Density plot of CAPE")
plt.xlabel('CAPE ')
plt.grid()