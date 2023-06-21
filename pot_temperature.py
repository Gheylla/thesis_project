# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 16:57:16 2023

@author: LENOVO
"""
'''Load all required packages'''


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

temp_nohgtqs= xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\ta_data\ta_Slev_noHGTQS.nc')
temp_nohgtqs_noshal = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\ta_data\ta_Slev_noHGTQS_noSHAL.nc') 
temp_nohgtqs_nouvmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\ta_data\ta_Slev_noHGTQS_noUVmix.nc')

p_nohgtqs = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\pres_data\p_Slev_noHGTQS.nc')
p_nohgtqs_noshal = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\pres_data\p_Slev_noHGTQS_noSHAL.nc')
p_nohgtqs_nouvmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\pres_data\p_Slev_noHGTQS_noUVmix.nc')


heights = pd.read_csv(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\z_43_fullheights.csv')

#%%
'''get the mean of the area'''

temp_nohgtqs_mean = temp_nohgtqs.mean(dim = ('x', 'y'))
temp_nohgtqs_noshal_mean = temp_nohgtqs_noshal.mean(dim = ('x', 'y'))
temp_nohgtqs_nouvmix_mean = temp_nohgtqs_nouvmix.mean(dim = ('x', 'y'))

p_nohgtqs = p_nohgtqs.mean(dim = ('x', 'y'))
p_nohgtqs_noshal = p_nohgtqs_noshal.mean(dim = ('x', 'y'))
p_nohgtqs_nouvmix = p_nohgtqs_nouvmix.mean(dim = ('x', 'y'))


#%%
'''Compute potentianl temperature'''
Rd = 287.04
cp = 1004.67
p0 = 1e5


pot_temp_nohgtqs = temp_nohgtqs_mean.ta * ( p0 / p_nohgtqs.p) ** (Rd / cp)
pot_temp_nohgtqs_noshal = temp_nohgtqs_noshal_mean.ta * ( p0 / p_nohgtqs_noshal.p) ** (Rd / cp)
pot_temp_nohgtqs_nouvmix = temp_nohgtqs_nouvmix_mean.ta * ( p0 / p_nohgtqs_nouvmix.p) ** (Rd / cp)


#%%
'''Get the mean of time to only leave the levels'''

pot_temp_nohgtqs = pot_temp_nohgtqs.mean(dim = ('time'))
pot_temp_nohgtqs_noshal = pot_temp_nohgtqs_noshal.mean(dim = ('time'))
pot_temp_nohgtqs_nouvmix = pot_temp_nohgtqs_nouvmix.mean(dim = ('time'))


#%%
'''Plot the data'''
plt.figure(figsize = (10,10))
plt.plot(np.flip(pot_temp_nohgtqs.values), heights['Full heights'], color = 'red', label = 'HARMONIE noHGTQS')
plt.plot(np.flip(pot_temp_nohgtqs_noshal.values), heights['Full heights'], color = 'blue', label = 'HARMONIE noHGTQS noSHAL')
plt.plot(np.flip(pot_temp_nohgtqs_nouvmix.values), heights['Full heights'], color = 'green', label = 'HARMONIE noHGTQS noUVmix')
plt.grid()
plt.title('Mean vertical profile of potential temperature')
plt.xlabel('${\theta}$ [K]')
plt.ylabel('Height [m]')
plt.legend()
plt.ylim([0,3000])
plt.xlim([295,315])








