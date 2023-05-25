# -*- coding: utf-8 -*-
"""
Created on Mon May  1 13:35:29 2023

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
import datetime as dt

'''Import data'''
data_hus_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\hus_data\hus_Slev_noHGTQS.nc')
data_hus_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\hus_data\hus_Slev_noHGTQS_noSHAL.nc')
data_hus_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\hus_data\hus_Slev_noHGTQS_noUVmix.nc')

heights = pd.read_csv(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\z_43_fullheights.csv')

# %%
'''Take mean of the area (remember that this is not the complete area we are analyzing but only the halo circle'''
hus_mean_noHGTQS = data_hus_noHGTQS.mean(dim = ('x', 'y'))
hus_mean_noHGTQS_noSHAL = data_hus_noHGTQS_noSHAL.mean(dim = ('x', 'y'))
hus_mean_noHGTQS_noUVmix = data_hus_noHGTQS_noUVmix.mean(dim = ('x', 'y'))

# %%
'''Make a diurnal overview of the cloud fraction'''
hus_diur_noHGTQS = hus_mean_noHGTQS.groupby(hus_mean_noHGTQS.time.dt.hour).mean()
hus_diur_noHGTQS_noSHAL = hus_mean_noHGTQS_noSHAL.groupby(hus_mean_noHGTQS_noSHAL.time.dt.hour).mean()
hus_diur_noHGTQS_noUVmix = hus_mean_noHGTQS_noUVmix.groupby(hus_mean_noHGTQS_noUVmix.time.dt.hour).mean()

# %%
'''Transpose the data and flip the values (because levels are flipped)'''
hus_diur_noHGTQS = hus_diur_noHGTQS.transpose()
hus_diur_noHGTQS_noSHAL = hus_diur_noHGTQS_noSHAL.transpose()
hus_diur_noHGTQS_noUVmix = hus_diur_noHGTQS_noUVmix.transpose()

hus_diur_noHGTQS_arr = np.array(hus_diur_noHGTQS.hus.values)
hus_diur_noHGTQS_arr = np.flip(hus_diur_noHGTQS_arr, axis = 0)

hus_diur_noHGTQS_noSHAL_arr = np.array(hus_diur_noHGTQS_noSHAL.hus.values)
hus_diur_noHGTQS_noSHAL_arr = np.flip(hus_diur_noHGTQS_noSHAL_arr, axis = 0)

hus_diur_noHGTQS_noUVmix_arr = np.array(hus_diur_noHGTQS_noUVmix.hus.values)
hus_diur_noHGTQS_noUVmix_arr = np.flip(hus_diur_noHGTQS_noUVmix_arr, axis = 0)

# %%
'''Plot the data for levels and diurnal cycle'''
hours = np.linspace(0,23,24)
#we still need to change the limit of the colorbar to set all of them at the same colorbar limit

plt.figure(figsize = (15,10))
plt.contourf(hus_diur_noHGTQS.hour.values, heights['Full heights'], hus_diur_noHGTQS_arr * 1000, cmap = 'BuGn')
plt.title('Composite diurnal cycle of specific humidity for HARMONIE noHGTQS')
plt.xlabel('Time [hr]')
plt.ylabel('Height [m]')
plt.grid()
plt.xticks(hours)
plt.ylim([0, 3000])
#plt.clim(vmin = 0, vmax = 100)
plt.colorbar(label= '$q_v$ [g/kg]')


plt.figure(figsize = (15,10))
plt.contourf(hus_diur_noHGTQS_noSHAL.hour.values, heights['Full heights'], hus_diur_noHGTQS_noSHAL_arr * 1000, cmap = 'BuGn')
plt.title('Composite diurnal cycle of specific humidity for HARMONIE noHGTQS noSHAL')
plt.xlabel('Time [hr]')
plt.ylabel('Height [m]')
plt.grid()
plt.xticks(hours)
plt.ylim([0, 3000])
plt.colorbar(label= '$q_v$ [g/kg]')
#plt.clim(vmin = 0, vmax = 1)



plt.figure(figsize = (15,10))
plt.contourf(hus_diur_noHGTQS_noUVmix.hour.values, heights['Full heights'], hus_diur_noHGTQS_noUVmix_arr * 1000, cmap = 'BuGn')
plt.title('Composite diurnal cycle of specific humidity for HARMONIE noHGTQS noUVmix')
plt.xlabel('Time [hr]')
plt.ylabel('Height [m]')
plt.grid()
plt.xticks(hours)
plt.ylim([0, 3000])
plt.colorbar(label= '$q_v$ [g/kg]')
#plt.clim(vmin = 0, vmax = 1)

# %%
'''Get the mean vertical profile of the hours that deviate the most'''
#hours 5, 6 and 7 
no_hgtqs_567 = hus_diur_noHGTQS.sel(hour = slice(5,7))
no_hgtqs_567_mean = no_hgtqs_567.mean(dim = 'hour')
no_hgtqs_567_mean = np.flip(no_hgtqs_567_mean.hus.values)

no_hgtqs_noshal_567 = hus_diur_noHGTQS_noSHAL.sel(hour = slice(5,7))
no_hgtqs_noshal_567_mean = no_hgtqs_noshal_567.mean(dim = 'hour')
no_hgtqs_noshal_567_mean = np.flip(no_hgtqs_noshal_567_mean.hus.values)

no_hgtqs_nouvmix_567 = hus_diur_noHGTQS_noUVmix.sel(hour = slice(5,7))
no_hgtqs_nouvmix_567_mean = no_hgtqs_nouvmix_567.mean(dim = 'hour')
no_hgtqs_nouvmix_567_mean = np.flip(no_hgtqs_nouvmix_567_mean.hus.values)

# %%
'''Get the mean vertical profile of the hours that deviate the most'''
#hours 5, 6 and 7 
no_hgtqs_1719 = hus_diur_noHGTQS.sel(hour = slice(17,19))
no_hgtqs_1719_mean = no_hgtqs_567.mean(dim = 'hour')
no_hgtqs_1719_mean = np.flip(no_hgtqs_1719_mean.hus.values)

no_hgtqs_noshal_1719 = hus_diur_noHGTQS_noSHAL.sel(hour = slice(17,19))
no_hgtqs_noshal_1719_mean = no_hgtqs_noshal_1719.mean(dim = 'hour')
no_hgtqs_noshal_1719_mean = np.flip(no_hgtqs_noshal_1719_mean.hus.values)

no_hgtqs_nouvmix_1719 = hus_diur_noHGTQS_noUVmix.sel(hour = slice(17,19))
no_hgtqs_nouvmix_1719_mean = no_hgtqs_nouvmix_1719.mean(dim = 'hour')
no_hgtqs_nouvmix_1719_mean = np.flip(no_hgtqs_nouvmix_1719_mean.hus.values)




# %%
'''Plot the mean of the hours'''
plt.figure(figsize = (10,10))
plt.plot(no_hgtqs_567_mean * 1000, heights['Full heights'], color = 'red', label = 'HARMONIE noHGTQS')
plt.plot(no_hgtqs_noshal_567_mean * 1000, heights['Full heights'], color = 'blue', label = 'HARMONIE noHGTQS noSHAL')
plt.plot(no_hgtqs_nouvmix_567_mean * 1000, heights['Full heights'], color = 'green', label = 'HARMONIE noHGTQS noUVmix')

plt.ylim([0, 3000])
plt.grid()
plt.title('Mean specific humidity for hours 5, 6 and 7 [UTC]')
plt.xlabel('$q_v$ [g/kg]')
plt.legend()
plt.ylabel('Height [m]')

plt.figure(figsize = (10,10))
plt.plot(no_hgtqs_567_mean * 1000, heights['Full heights'], color = 'red', label = 'HARMONIE noHGTQS')
plt.plot(no_hgtqs_noshal_567_mean * 1000, heights['Full heights'], color = 'blue', label = 'HARMONIE noHGTQS noSHAL')
plt.plot(no_hgtqs_nouvmix_567_mean * 1000, heights['Full heights'], color = 'green', label = 'HARMONIE noHGTQS noUVmix')

plt.ylim([0, 500])
plt.xlim([12, 18])
plt.grid()
plt.title('Mean specific humidity for hours 5, 6 and 7 [UTC]')
plt.xlabel('$q_v$ [g/kg]')
plt.legend()
plt.ylabel('Height [m]')

'''Plot the mean of the hours'''
plt.figure(figsize = (10,10))
plt.plot(no_hgtqs_1719_mean * 1000, heights['Full heights'], color = 'red', label = 'HARMONIE noHGTQS')
plt.plot(no_hgtqs_noshal_1719_mean * 1000, heights['Full heights'], color = 'blue', label = 'HARMONIE noHGTQS noSHAL')
plt.plot(no_hgtqs_nouvmix_1719_mean * 1000, heights['Full heights'], color = 'green', label = 'HARMONIE noHGTQS noUVmix')

plt.ylim([0, 3000])
plt.grid()
plt.title('Mean specific humidity for hours 17, 18 and 19 [UTC]')
plt.xlabel('$q_v$ [g/kg]')
plt.legend()
plt.ylabel('Height [m]')

plt.figure(figsize = (10,10))
plt.plot(no_hgtqs_1719_mean * 1000, heights['Full heights'], color = 'red', label = 'HARMONIE noHGTQS')
plt.plot(no_hgtqs_noshal_1719_mean * 1000, heights['Full heights'], color = 'blue', label = 'HARMONIE noHGTQS noSHAL')
plt.plot(no_hgtqs_nouvmix_1719_mean * 1000, heights['Full heights'], color = 'green', label = 'HARMONIE noHGTQS noUVmix')

plt.ylim([0, 500])
plt.xlim([12, 18])
plt.grid()
plt.title('Mean specific humidity for hours 17, 18 and 19 [UTC]')
plt.xlabel('$q_v$ [g/kg]')
plt.legend()
plt.ylabel('Height [m]')


# %%
'''Select only the hours where the cloud fraction deviates the most'''
hr5_loc = np.where(hus_diur_noHGTQS.hour.values == 5)[0][0]

hus_noHGTQS_hr5 = hus_diur_noHGTQS.hus.values[:, hr5_loc]
hus_noHGTQS_hr5 = np.array(hus_noHGTQS_hr5)
hus_noHGTQS_hr5 = np.flip(hus_noHGTQS_hr5)

hus_noHGTQS_noSHAL_hr5 = hus_diur_noHGTQS_noSHAL.hus.values[:, hr5_loc]
hus_noHGTQS_noSHAL_hr5 = np.array(hus_noHGTQS_noSHAL_hr5)
hus_noHGTQS_noSHAL_hr5 = np.flip(hus_noHGTQS_noSHAL_hr5)

hus_noHGTQS_noUVmix_hr5 = hus_diur_noHGTQS_noUVmix.hus.values[:, hr5_loc]
hus_noHGTQS_noUVmix_hr5 = np.array(hus_noHGTQS_noUVmix_hr5)
hus_noHGTQS_noUVmix_hr5 = np.flip(hus_noHGTQS_noUVmix_hr5)

plt.figure(figsize = (10,10))
plt.plot(hus_noHGTQS_hr5, heights['Full heights'], label = 'HARMONIE noHGTQS')
plt.plot(hus_noHGTQS_noSHAL_hr5, heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL' )
plt.plot(hus_noHGTQS_noUVmix_hr5, heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix' )
plt.ylim([0, 3000])
plt.legend()
plt.title('Specific Humidity at time = 5 ')
plt.ylabel('Height [m]')
plt.xlabel('$q_v$ [g/kg]')
plt.grid()

# %%

