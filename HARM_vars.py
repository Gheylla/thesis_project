# -*- coding: utf-8 -*-
"""
Created on Mon May  1 10:32:35 2023

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
data_cl_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cf_data\cl_Slev_noHGTQS.nc')
data_cl_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cf_data\cl_Slev_noHGTQS_noSHAL.nc')
data_cl_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cf_data\cl_Slev_noHGTQS_noUVmix.nc')

heights = pd.read_csv(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\z_43_fullheights.csv')

# %%
'''Take mean of the area (remember that this is not the complete area we are analyzing but only the halo circle'''
cl_mean_noHGTQS = data_cl_noHGTQS.mean(dim = ('x', 'y'))
cl_mean_noHGTQS_noSHAL = data_cl_noHGTQS_noSHAL.mean(dim = ('x', 'y'))
cl_mean_noHGTQS_noUVmix = data_cl_noHGTQS_noUVmix.mean(dim = ('x', 'y'))

# %%
'''Make a diurnal overview of the cloud fraction'''
cl_diur_noHGTQS = cl_mean_noHGTQS.groupby(cl_mean_noHGTQS.time.dt.hour).mean()
cl_diur_noHGTQS_noSHAL = cl_mean_noHGTQS_noSHAL.groupby(cl_mean_noHGTQS_noSHAL.time.dt.hour).mean()
cl_diur_noHGTQS_noUVmix = cl_mean_noHGTQS_noUVmix.groupby(cl_mean_noHGTQS_noUVmix.time.dt.hour).mean()


cl_diur_noHGTQS = cl_diur_noHGTQS.transpose()
cl_diur_noHGTQS_noSHAL = cl_diur_noHGTQS_noSHAL.transpose()
cl_diur_noHGTQS_noUVmix = cl_diur_noHGTQS_noUVmix.transpose()



cl_diur_noHGTQS_arr = np.array(cl_diur_noHGTQS.cl.values)
cl_diur_noHGTQS_arr = np.flip(cl_diur_noHGTQS_arr, axis = 0)

cl_diur_noHGTQS_noSHAL_arr = np.array(cl_diur_noHGTQS_noSHAL.cl.values)
cl_diur_noHGTQS_noSHAL_arr = np.flip(cl_diur_noHGTQS_noSHAL_arr, axis = 0)

cl_diur_noHGTQS_noUVmix_arr = np.array(cl_diur_noHGTQS_noUVmix.cl.values)
cl_diur_noHGTQS_noUVmix_arr = np.flip(cl_diur_noHGTQS_noUVmix_arr, axis = 0)


# %%
'''Plotting the data'''
hours = np.linspace(0,23,24)
#we still need to change the limit of the colorbar to set all of them at the same colorbar limit

plt.figure(figsize = (15,10))
plt.contourf(cl_diur_noHGTQS.hour.values, heights['Full heights'], cl_diur_noHGTQS_arr * 100, cmap = 'Reds', vmin = 0, vmax = 100)
plt.title('Composite diurnal cycle of cloud fraction for HARMONIE noHGTQS')
plt.xlabel('Time [hr]')
plt.ylabel('Height [m]')
plt.grid()
plt.xticks(hours)
plt.ylim([0, 3000])
plt.clim(vmin = 0, vmax = 100)
plt.colorbar(label= 'CF [%]', ticks = np.arange(0, 110, 10))


plt.figure(figsize = (15,10))
plt.contourf(cl_diur_noHGTQS_noSHAL.hour.values, heights['Full heights'], cl_diur_noHGTQS_noSHAL_arr * 100, cmap = 'Reds')
plt.title('Composite diurnal cycle of cloud fraction for HARMONIE noHGTQS noSHAL')
plt.xlabel('Time [hr]')
plt.ylabel('Height [m]')
plt.grid()
plt.xticks(hours)
plt.ylim([0, 3000])
plt.colorbar(label= 'CF [-]')
#plt.clim(vmin = 0, vmax = 1)



plt.figure(figsize = (15,10))
plt.contourf(cl_diur_noHGTQS_noUVmix.hour.values, heights['Full heights'], cl_diur_noHGTQS_noUVmix_arr * 100, cmap = 'Reds')
plt.title('Composite diurnal cycle of cloud fraction for HARMONIE noHGTQS noUVmix')
plt.xlabel('Time [hr]')
plt.ylabel('Height [m]')
plt.grid()
plt.xticks(hours)
plt.ylim([0, 3000])
plt.colorbar(label= 'CF [-]')
#plt.clim(vmin = 0, vmax = 1)

# %%
'''Plot for specific hours that are the most different'''
'''Hour 5'''
hr5_loc = np.where(cl_diur_noHGTQS.hour.values == 5)[0][0]
cl_noHGTQS_hr5 = cl_diur_noHGTQS.cl.values[:, hr5_loc]
cl_noHGTQS_hr5 = np.array(cl_noHGTQS_hr5)
cl_noHGTQS_hr5 = np.flip(cl_noHGTQS_hr5)

cl_noHGTQS_noSHAL_hr5 = cl_diur_noHGTQS_noSHAL.cl.values[:, hr5_loc]
cl_noHGTQS_noSHAL_hr5 = np.array(cl_noHGTQS_noSHAL_hr5)
cl_noHGTQS_noSHAL_hr5 = np.flip(cl_noHGTQS_noSHAL_hr5)

cl_noHGTQS_noUVmix_hr5 = cl_diur_noHGTQS_noUVmix.cl.values[:, hr5_loc]
cl_noHGTQS_noUVmix_hr5 = np.array(cl_noHGTQS_noUVmix_hr5)
cl_noHGTQS_noUVmix_hr5 = np.flip(cl_noHGTQS_noUVmix_hr5)

plt.figure(figsize = (10,10))
plt.plot(cl_noHGTQS_hr5, heights['Full heights'], label = 'HARMONIE noHGTQS')
plt.plot(cl_noHGTQS_noSHAL_hr5, heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL' )
plt.plot(cl_noHGTQS_noUVmix_hr5, heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix' )
plt.ylim([0, 3000])
plt.legend()
plt.title('Cloudfraction at time = 5 ')
plt.ylabel('Height [m]')
plt.xlabel('CF [-]')
plt.grid()

# %%
'''Hour 6'''
hr6_loc = np.where(cl_diur_noHGTQS.hour.values == 6)[0][0]
cl_noHGTQS_hr6 = cl_diur_noHGTQS.cl.values[:, hr6_loc]
cl_noHGTQS_hr6 = np.array(cl_noHGTQS_hr6)
cl_noHGTQS_hr6 = np.flip(cl_noHGTQS_hr6)

cl_noHGTQS_noSHAL_hr6 = cl_diur_noHGTQS_noSHAL.cl.values[:, hr6_loc]
cl_noHGTQS_noSHAL_hr6 = np.array(cl_noHGTQS_noSHAL_hr6)
cl_noHGTQS_noSHAL_hr6 = np.flip(cl_noHGTQS_noSHAL_hr6)

cl_noHGTQS_noUVmix_hr6 = cl_diur_noHGTQS_noUVmix.cl.values[:, hr6_loc]
cl_noHGTQS_noUVmix_hr6 = np.array(cl_noHGTQS_noUVmix_hr6)
cl_noHGTQS_noUVmix_hr6 = np.flip(cl_noHGTQS_noUVmix_hr6)

plt.figure(figsize = (10,10))
plt.plot(cl_noHGTQS_hr6, heights['Full heights'], label = 'HARMONIE noHGTQS')
plt.plot(cl_noHGTQS_noSHAL_hr6, heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL' )
plt.plot(cl_noHGTQS_noUVmix_hr6, heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix' )
plt.ylim([0, 3000])
plt.legend()
plt.title('Cloudfraction at time = 6 ')
plt.ylabel('Height [m]')
plt.xlabel('CF [-]')
plt.grid()

# %%
'''Hour 7'''
hr7_loc = np.where(cl_diur_noHGTQS.hour.values == 7)[0][0]
cl_noHGTQS_hr7 = cl_diur_noHGTQS.cl.values[:, hr7_loc]
cl_noHGTQS_hr7 = np.array(cl_noHGTQS_hr7)
cl_noHGTQS_hr7 = np.flip(cl_noHGTQS_hr7)

cl_noHGTQS_noSHAL_hr7 = cl_diur_noHGTQS_noSHAL.cl.values[:, hr7_loc]
cl_noHGTQS_noSHAL_hr7 = np.array(cl_noHGTQS_noSHAL_hr7)
cl_noHGTQS_noSHAL_hr7 = np.flip(cl_noHGTQS_noSHAL_hr7)

cl_noHGTQS_noUVmix_hr7 = cl_diur_noHGTQS_noUVmix.cl.values[:, hr7_loc]
cl_noHGTQS_noUVmix_hr7 = np.array(cl_noHGTQS_noUVmix_hr7)
cl_noHGTQS_noUVmix_hr7 = np.flip(cl_noHGTQS_noUVmix_hr7)

plt.figure(figsize = (10,10))
plt.plot(cl_noHGTQS_hr7, heights['Full heights'], label = 'HARMONIE noHGTQS')
plt.plot(cl_noHGTQS_noSHAL_hr7, heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL' )
plt.plot(cl_noHGTQS_noUVmix_hr7, heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix' )
plt.ylim([0, 3000])
plt.legend()
plt.title('Cloudfraction at time = 7 ')
plt.ylabel('Height [m]')
plt.xlabel('CF [-]')
plt.grid()

# %%
'''Hour 17'''
hr17_loc = np.where(cl_diur_noHGTQS.hour.values == 17)[0][0]
cl_noHGTQS_hr17 = cl_diur_noHGTQS.cl.values[:, hr17_loc]
cl_noHGTQS_hr17 = np.array(cl_noHGTQS_hr17)
cl_noHGTQS_hr17 = np.flip(cl_noHGTQS_hr17)

cl_noHGTQS_noSHAL_hr17 = cl_diur_noHGTQS_noSHAL.cl.values[:, hr17_loc]
cl_noHGTQS_noSHAL_hr17 = np.array(cl_noHGTQS_noSHAL_hr17)
cl_noHGTQS_noSHAL_hr17 = np.flip(cl_noHGTQS_noSHAL_hr17)

cl_noHGTQS_noUVmix_hr17 = cl_diur_noHGTQS_noUVmix.cl.values[:, hr17_loc]
cl_noHGTQS_noUVmix_hr17 = np.array(cl_noHGTQS_noUVmix_hr17)
cl_noHGTQS_noUVmix_hr17 = np.flip(cl_noHGTQS_noUVmix_hr17)

plt.figure(figsize = (10,10))
plt.plot(cl_noHGTQS_hr17, heights['Full heights'], label = 'HARMONIE noHGTQS')
plt.plot(cl_noHGTQS_noSHAL_hr17, heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL' )
plt.plot(cl_noHGTQS_noUVmix_hr17, heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix' )
plt.ylim([0, 3000])
plt.legend()
plt.title('Cloudfraction at time = 17 ')
plt.ylabel('Height [m]')
plt.xlabel('CF [-]')
plt.grid()

# %%
'''Hour 18'''
hr18_loc = np.where(cl_diur_noHGTQS.hour.values == 18)[0][0]
cl_noHGTQS_hr18 = cl_diur_noHGTQS.cl.values[:, hr18_loc]
cl_noHGTQS_hr18 = np.array(cl_noHGTQS_hr18)
cl_noHGTQS_hr18 = np.flip(cl_noHGTQS_hr18)

cl_noHGTQS_noSHAL_hr18 = cl_diur_noHGTQS_noSHAL.cl.values[:, hr18_loc]
cl_noHGTQS_noSHAL_hr18 = np.array(cl_noHGTQS_noSHAL_hr18)
cl_noHGTQS_noSHAL_hr18 = np.flip(cl_noHGTQS_noSHAL_hr18)

cl_noHGTQS_noUVmix_hr18 = cl_diur_noHGTQS_noUVmix.cl.values[:, hr18_loc]
cl_noHGTQS_noUVmix_hr18 = np.array(cl_noHGTQS_noUVmix_hr18)
cl_noHGTQS_noUVmix_hr18 = np.flip(cl_noHGTQS_noUVmix_hr18)

plt.figure(figsize = (10,10))
plt.plot(cl_noHGTQS_hr18, heights['Full heights'], label = 'HARMONIE noHGTQS')
plt.plot(cl_noHGTQS_noSHAL_hr18, heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL' )
plt.plot(cl_noHGTQS_noUVmix_hr18, heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix' )
plt.ylim([0, 3000])
plt.legend()
plt.title('Cloudfraction at time = 18 ')
plt.ylabel('Height [m]')
plt.xlabel('CF [-]')
plt.grid()


# %%
'''Hour 19'''
hr19_loc = np.where(cl_diur_noHGTQS.hour.values == 19)[0][0]
cl_noHGTQS_hr19 = cl_diur_noHGTQS.cl.values[:, hr19_loc]
cl_noHGTQS_hr19 = np.array(cl_noHGTQS_hr19)
cl_noHGTQS_hr19 = np.flip(cl_noHGTQS_hr19)

cl_noHGTQS_noSHAL_hr19 = cl_diur_noHGTQS_noSHAL.cl.values[:, hr19_loc]
cl_noHGTQS_noSHAL_hr19 = np.array(cl_noHGTQS_noSHAL_hr19)
cl_noHGTQS_noSHAL_hr19 = np.flip(cl_noHGTQS_noSHAL_hr19)

cl_noHGTQS_noUVmix_hr19 = cl_diur_noHGTQS_noUVmix.cl.values[:, hr19_loc]
cl_noHGTQS_noUVmix_hr19 = np.array(cl_noHGTQS_noUVmix_hr19)
cl_noHGTQS_noUVmix_hr19 = np.flip(cl_noHGTQS_noUVmix_hr19)

plt.figure(figsize = (10,10))
plt.plot(cl_noHGTQS_hr19, heights['Full heights'], label = 'HARMONIE noHGTQS')
plt.plot(cl_noHGTQS_noSHAL_hr19, heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL' )
plt.plot(cl_noHGTQS_noUVmix_hr19, heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix' )
plt.ylim([0, 3000])
plt.legend()
plt.title('Cloudfraction at time = 19 ')
plt.ylabel('Height [m]')
plt.xlabel('CF [-]')
plt.grid()

