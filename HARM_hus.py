# -*- coding: utf-8 -*-
"""
Created on Mon May  1 13:35:29 2023

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

'''Hour 6'''
hr6_loc = np.where(hus_diur_noHGTQS.hour.values == 6)[0][0]
hus_noHGTQS_hr6 = hus_diur_noHGTQS.hus.values[:, hr6_loc]
hus_noHGTQS_hr6 = np.array(hus_noHGTQS_hr6)
hus_noHGTQS_hr6 = np.flip(hus_noHGTQS_hr6)

hus_noHGTQS_noSHAL_hr6 = hus_diur_noHGTQS_noSHAL.hus.values[:, hr6_loc]
hus_noHGTQS_noSHAL_hr6 = np.array(hus_noHGTQS_noSHAL_hr6)
hus_noHGTQS_noSHAL_hr6 = np.flip(hus_noHGTQS_noSHAL_hr6)

hus_noHGTQS_noUVmix_hr6 = hus_diur_noHGTQS_noUVmix.hus.values[:, hr6_loc]
hus_noHGTQS_noUVmix_hr6 = np.array(hus_noHGTQS_noUVmix_hr6)
hus_noHGTQS_noUVmix_hr6 = np.flip(hus_noHGTQS_noUVmix_hr6)

plt.figure(figsize = (10,10))
plt.plot(hus_noHGTQS_hr6 * 1000, heights['Full heights'], label = 'HARMONIE noHGTQS')
plt.plot(hus_noHGTQS_noSHAL_hr6 * 1000, heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL' )
plt.plot(hus_noHGTQS_noUVmix_hr6 * 1000, heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix' )
plt.ylim([0, 3000])
plt.legend()
plt.title('Specific Humidity at time = 6 ')
plt.ylabel('Height [m]')
plt.xlabel('$q_v$ [g/kg]')
plt.grid()

# %%

'''Hour 7'''
hr7_loc = np.where(hus_diur_noHGTQS.hour.values == 7)[0][0]
hus_noHGTQS_hr7 = hus_diur_noHGTQS.hus.values[:, hr7_loc]
hus_noHGTQS_hr7 = np.array(hus_noHGTQS_hr7)
hus_noHGTQS_hr7 = np.flip(hus_noHGTQS_hr7)

hus_noHGTQS_noSHAL_hr7 = hus_diur_noHGTQS_noSHAL.hus.values[:, hr7_loc]
hus_noHGTQS_noSHAL_hr7 = np.array(hus_noHGTQS_noSHAL_hr7)
hus_noHGTQS_noSHAL_hr7 = np.flip(hus_noHGTQS_noSHAL_hr7)

hus_noHGTQS_noUVmix_hr7 = hus_diur_noHGTQS_noUVmix.hus.values[:, hr7_loc]
hus_noHGTQS_noUVmix_hr7 = np.array(hus_noHGTQS_noUVmix_hr7)
hus_noHGTQS_noUVmix_hr7 = np.flip(hus_noHGTQS_noUVmix_hr7)

plt.figure(figsize = (10,10))
plt.plot(hus_noHGTQS_hr7 * 1000, heights['Full heights'], label = 'HARMONIE noHGTQS')
plt.plot(hus_noHGTQS_noSHAL_hr7 * 1000, heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL' )
plt.plot(hus_noHGTQS_noUVmix_hr7 * 1000, heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix' )
plt.ylim([0, 3000])
plt.legend()
plt.title('Specific Humidity at time = 7 ')
plt.ylabel('Height [m]')
plt.xlabel('$q_v$ [g/kg]')
plt.grid()

# %%

'''Hour 17'''
hr17_loc = np.where(hus_diur_noHGTQS.hour.values == 17)[0][0]
hus_noHGTQS_hr17 = hus_diur_noHGTQS.hus.values[:, hr17_loc]
hus_noHGTQS_hr17 = np.array(hus_noHGTQS_hr17)
hus_noHGTQS_hr17 = np.flip(hus_noHGTQS_hr17)

hus_noHGTQS_noSHAL_hr17 = hus_diur_noHGTQS_noSHAL.hus.values[:, hr17_loc]
hus_noHGTQS_noSHAL_hr17 = np.array(hus_noHGTQS_noSHAL_hr17)
hus_noHGTQS_noSHAL_hr17 = np.flip(hus_noHGTQS_noSHAL_hr17)

hus_noHGTQS_noUVmix_hr17 = hus_diur_noHGTQS_noUVmix.hus.values[:, hr17_loc]
hus_noHGTQS_noUVmix_hr17 = np.array(hus_noHGTQS_noUVmix_hr17)
hus_noHGTQS_noUVmix_hr17 = np.flip(hus_noHGTQS_noUVmix_hr17)

plt.figure(figsize = (10,10))
plt.plot(hus_noHGTQS_hr17 * 1000, heights['Full heights'], label = 'HARMONIE noHGTQS')
plt.plot(hus_noHGTQS_noSHAL_hr17 * 1000, heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL' )
plt.plot(hus_noHGTQS_noUVmix_hr17 * 10000, heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix' )
plt.ylim([0, 3000])
plt.legend()
plt.title('Specific Humidity at time = 17 ')
plt.ylabel('Height [m]')
plt.xlabel('$q_v$ [g/kg]')
plt.grid()

# %%

'''Hour 18'''
hr18_loc = np.where(hus_diur_noHGTQS.hour.values == 18)[0][0]
hus_noHGTQS_hr18 = hus_diur_noHGTQS.hus.values[:, hr18_loc]
hus_noHGTQS_hr18 = np.array(hus_noHGTQS_hr18)
hus_noHGTQS_hr18 = np.flip(hus_noHGTQS_hr18)

hus_noHGTQS_noSHAL_hr18 = hus_diur_noHGTQS_noSHAL.hus.values[:, hr18_loc]
hus_noHGTQS_noSHAL_hr18 = np.array(hus_noHGTQS_noSHAL_hr18)
hus_noHGTQS_noSHAL_hr18 = np.flip(hus_noHGTQS_noSHAL_hr18)

hus_noHGTQS_noUVmix_hr18 = hus_diur_noHGTQS_noUVmix.hus.values[:, hr18_loc]
hus_noHGTQS_noUVmix_hr18 = np.array(hus_noHGTQS_noUVmix_hr18)
hus_noHGTQS_noUVmix_hr18 = np.flip(hus_noHGTQS_noUVmix_hr18)

plt.figure(figsize = (10,10))
plt.plot(hus_noHGTQS_hr18 * 1000, heights['Full heights'], label = 'HARMONIE noHGTQS')
plt.plot(hus_noHGTQS_noSHAL_hr18 * 1000, heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL' )
plt.plot(hus_noHGTQS_noUVmix_hr18 * 1000, heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix' )
plt.ylim([0, 3000])
plt.legend()
plt.title('Specific Humidity at time = 18 ')
plt.ylabel('Height [m]')
plt.xlabel('$q_v$ [g/kg]')
plt.grid()

# %%

'''Hour 19'''
hr19_loc = np.where(hus_diur_noHGTQS.hour.values == 19)[0][0]
hus_noHGTQS_hr19 = hus_diur_noHGTQS.hus.values[:, hr19_loc]
hus_noHGTQS_hr19 = np.array(hus_noHGTQS_hr19)
hus_noHGTQS_hr19 = np.flip(hus_noHGTQS_hr19)

hus_noHGTQS_noSHAL_hr19 = hus_diur_noHGTQS_noSHAL.hus.values[:, hr19_loc]
hus_noHGTQS_noSHAL_hr19 = np.array(hus_noHGTQS_noSHAL_hr19)
hus_noHGTQS_noSHAL_hr19 = np.flip(hus_noHGTQS_noSHAL_hr19)

hus_noHGTQS_noUVmix_hr19 = hus_diur_noHGTQS_noUVmix.hus.values[:, hr19_loc]
hus_noHGTQS_noUVmix_hr19 = np.array(hus_noHGTQS_noUVmix_hr19)
hus_noHGTQS_noUVmix_hr19 = np.flip(hus_noHGTQS_noUVmix_hr19)

plt.figure(figsize = (10,10))
plt.plot(hus_noHGTQS_hr19 * 1000, heights['Full heights'], label = 'HARMONIE noHGTQS')
plt.plot(hus_noHGTQS_noSHAL_hr19 * 1000, heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL' )
plt.plot(hus_noHGTQS_noUVmix_hr19 * 1000, heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix' )
plt.ylim([0, 3000])
plt.legend()
plt.title('Specific Humidity at time = 19 ')
plt.ylabel('Height [m]')
plt.xlabel('$q_v$ [g/kg]')
plt.grid()