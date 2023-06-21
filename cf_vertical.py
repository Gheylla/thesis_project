# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 09:20:51 2023

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

#%%
'''Import data'''
cf_nohgtqs = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cf_data\cl_Slev_noHGTQS.nc')
cf_nohgtqs_noshal = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cf_data\cl_Slev_noHGTQS_noSHAL.nc')
cf_nohgtqs_nouvmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cf_data\cl_Slev_noHGTQS_noUVmix.nc')

heights = pd.read_csv(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\z_43_fullheights.csv')

#%%
'''Determine the mean of the area'''
cf_nohgtqs_mean = cf_nohgtqs.mean(dim = ('x', 'y'))
cf_nohgtqs_noshal_mean = cf_nohgtqs_noshal.mean(dim = ('x', 'y'))
cf_nohgtqs_nouvmix_mean = cf_nohgtqs_nouvmix.mean(dim = ('x', 'y'))

#%%
'''Determine the diurnal cycle'''
cf_nohgtqs_diur = cf_nohgtqs_mean.groupby(cf_nohgtqs_mean.time.dt.hour).mean()
cf_nohgtqs_noshal_diur = cf_nohgtqs_noshal_mean.groupby(cf_nohgtqs_noshal_mean.time.dt.hour).mean()
cf_nohgtqs_nouvmix_diur = cf_nohgtqs_nouvmix_mean.groupby(cf_nohgtqs_nouvmix_mean.time.dt.hour).mean()

#%%
'''Transpose the data'''
cf_nohgtqs_diur_trans = cf_nohgtqs_diur.transpose()
cf_nohgtqs_noshal_diur_trans = cf_nohgtqs_noshal_diur.transpose()
cf_nohgtqs_nouvmix_diur_trans = cf_nohgtqs_nouvmix_diur.transpose()



#%%
'''Plot the data as contour'''
hours = np.linspace(00,23, 24)
lt = [20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
fig, ax = plt.subplots(layout='constrained', figsize = (15,5))
levels = np.linspace(0,0.3,10)

gr = ax.contourf(cf_nohgtqs_diur_trans.hour.values, heights['Full heights'], np.flip(cf_nohgtqs_diur_trans.cl.values, axis = 0), cmap = 'Blues_r', vmin = 0, vmax = 0.2, levels =levels)


plt.ylim([0,3000])
ax.grid()
ax.set_xticks(hours)
ax.set_xlabel('UTC Time [hr]')
secax = ax.secondary_xaxis('top')
secax.set_xticks(hours, lt)
secax.set_xlabel('Local Time [hr]')
fig.colorbar(gr, label = 'CF [-]')
plt.ylabel('Height [m]')
plt.title('Composite Diurnal cycle of Cloud fraction in the atmosphere for HARMONIE noHGTQS')
plt.show()

#%%
fig, ax = plt.subplots(layout='constrained', figsize = (15,5))
fr = ax.contourf(cf_nohgtqs_noshal_diur_trans.hour.values, heights['Full heights'], np.flip(cf_nohgtqs_noshal_diur_trans.cl.values, axis = 0), cmap = 'Blues_r', vmin = 0, vmax = 0.2, levels =levels)


plt.ylim([0,3000])
ax.grid()
ax.set_xticks(hours)
ax.set_xlabel('UTC Time [hr]')
secax = ax.secondary_xaxis('top')
secax.set_xticks(hours, lt)
secax.set_xlabel('Local Time [hr]')
fig.colorbar(fr, label = 'CF [-]')


plt.ylabel('Height [m]')
plt.title('Composite Diurnal cycle of Cloud fraction in the atmosphere for HARMONIE noHGTQS noSHAL')
plt.show()

#%%
fig, ax = plt.subplots(layout='constrained', figsize = (15,5))
lr = ax.contourf(cf_nohgtqs_nouvmix_diur_trans.hour.values, heights['Full heights'], np.flip(cf_nohgtqs_nouvmix_diur_trans.cl.values, axis = 0), cmap = 'Blues_r', vmin = 0, vmax = 0.2, levels =levels)


plt.ylim([0,3000])
ax.grid()
ax.set_xticks(hours)
ax.set_xlabel('UTC Time [hr]')
secax = ax.secondary_xaxis('top')
secax.set_xticks(hours, lt)
secax.set_xlabel('Local Time [hr]')
fig.colorbar(lr, label = 'CF [-]')


plt.ylabel('Height [m]')
plt.title('Composite Diurnal cycle of Cloud fraction in the atmosphere for HARMONIE noHGTQS noUVmix')
plt.show()

#%%
'''get the mean of the vertical profile'''
cf_nohgtqs_vert = cf_nohgtqs_mean.mean(dim = 'time')
cf_nohgtqs_noshal_vert = cf_nohgtqs_noshal_mean.mean(dim = 'time')
cf_nohgtqs_nouvmix_vert = cf_nohgtqs_nouvmix_mean.mean(dim = 'time')

#%%
'''Plot the mean vertical profile of cloud fraction'''
plt.figure(figsize = (10,10))

plt.plot(np.flip(cf_nohgtqs_vert.cl.values), heights['Full heights'], label = 'HARMONIE noHGTQS', color = 'red')
plt.plot(np.flip(cf_nohgtqs_noshal_vert.cl.values), heights['Full heights'], label = 'HARMONIE noHGTQS noSHAL', color = 'blue')
plt.plot(np.flip(cf_nohgtqs_nouvmix_vert.cl.values), heights['Full heights'], label = 'HARMONIE noHGTQS noUVmix', color = 'green')
plt.ylim([0,3000])
plt.grid()
plt.legend(loc = 'upper right')
plt.xlabel('CF [-]')
plt.ylabel('Height [m]')
plt.title('Mean vertical profile of cloud fraction')