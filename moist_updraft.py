# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 16:20:28 2023

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

moist_nohgtqs = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\moist_updraft_data\moist_updraft_mf_noHGTQS.nc')
moist_nohgtqs_noshal = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\moist_updraft_data\moist_updraft_mf_noHGTQS_noSHAL.nc')

#%%
heights = pd.read_csv(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\z_43_fullheights.csv')


#%%
'''Get the means of the area'''
moist_nohgtqs_mean = moist_nohgtqs.mean(dim = ('x', 'y'))
moist_nohgtqs_diur = moist_nohgtqs_mean.groupby(moist_nohgtqs_mean.time.dt.hour).mean()
moist_nohgtqs_diur = moist_nohgtqs_diur.transpose()


moist_nohgtqs_noshal_mean = moist_nohgtqs_noshal.mean(dim = ('x', 'y'))
moist_nohgtqs_noshal_diur = moist_nohgtqs_noshal_mean.groupby(moist_nohgtqs_noshal_mean.time.dt.hour).mean()
moist_nohgtqs_noshal_diur = moist_nohgtqs_noshal_diur.transpose()

#%%
'''Plot the graphs'''
hours = np.linspace(00,23, 24)
lt = [20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
fig, ax = plt.subplots(layout='constrained', figsize = (15,5))

gr = ax.contourf(moist_nohgtqs_diur.hour.values, heights['Full heights'], np.flip(moist_nohgtqs_diur.moist_updraft_mf.values, axis = 0))


#ax.colorbar(label = 'moist updraft MF [kg/m2s')
plt.ylim([0,3000])
ax.grid()
ax.set_xticks(hours)
ax.set_xlabel('UTC Time [hr]')
secax = ax.secondary_xaxis('top')
secax.set_xticks(hours, lt)
secax.set_xlabel('Local Time [hr]')
fig.colorbar(gr, label = 'moist updraft MF [kg/m2s]')


plt.ylabel('Height [m]')
plt.title('Composite Diurnal cycle of moist updraft for HARMONIE noHGTQS')
plt.show()






#%%





plt.figure()
plt.contourf( moist_nohgtqs_noshal_diur.lev.values, moist_nohgtqs_noshal_diur.hour.values, moist_nohgtqs_noshal_diur.moist_updraft_mf.values)
plt.colorbar(label = 'moist updraft MF [kg/m2s')
