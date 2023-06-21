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

dry_nohgtqs = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\dry_updraft_data\dry_updraft_mf_noHGTQS.nc')
dry_nohgtqs_noshal = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\dry_updraft_data\dry_updraft_mf_noHGTQS_noSHAL.nc')

#%%
heights = pd.read_csv(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\z_43_fullheights.csv')


#%%
'''Get the means of the area'''
dry_nohgtqs_mean = dry_nohgtqs.mean(dim = ('x', 'y'))
dry_nohgtqs_diur = dry_nohgtqs_mean.groupby(dry_nohgtqs_mean.time.dt.hour).mean()
dry_nohgtqs_diur = dry_nohgtqs_diur.transpose()


dry_nohgtqs_noshal_mean = dry_nohgtqs_noshal.mean(dim = ('x', 'y'))
dry_nohgtqs_noshal_diur = dry_nohgtqs_noshal_mean.groupby(dry_nohgtqs_noshal_mean.time.dt.hour).mean()
dry_nohgtqs_noshal_diur = dry_nohgtqs_noshal_diur.transpose()

#%%
'''Plot the graphs'''
hours = np.linspace(00,23, 24)
lt = [20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
fig, ax = plt.subplots(layout='constrained', figsize = (15,5))

gr = ax.contourf(dry_nohgtqs_diur.hour.values, heights['Full heights'], np.flip(dry_nohgtqs_diur.dry_updraft_mf.values, axis = 0))


#ax.colorbar(label = 'Dry updraft MF [kg/m2s')
plt.ylim([0,3000])
ax.grid()
ax.set_xticks(hours)
ax.set_xlabel('UTC Time [hr]')
secax = ax.secondary_xaxis('top')
secax.set_xticks(hours, lt)
secax.set_xlabel('Local Time [hr]')
fig.colorbar(gr, label = 'Dry updraft MF [kg/m2s]')


plt.ylabel('Height [m]')
plt.title('Composite Diurnal cycle of dry updraft for HARMONIE noHGTQS')
plt.show()












# =============================================================================
# plt.figure()
# plt.contourf( dry_nohgtqs_noshal_diur.lev.values, dry_nohgtqs_noshal_diur.hour.values, dry_nohgtqs_noshal_diur.dry_updraft_mf.values)
# plt.colorbar(label = 'Dry updraft MF [kg/m2s')
# =============================================================================
