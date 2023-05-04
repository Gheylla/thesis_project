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
cape_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cape_data\cape_noHGTQS.nc')
cape_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cape_data\cape_noHGTQS_noSHAL.nc')
cape_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cape_data\cape_noHGTQS_noUVmix.nc')

# %%
'''Cut data for desired grid'''
cape_noHGTQS = cut_sub_grid(cape_noHGTQS)
cape_noHGTQS_noSHAL = cut_sub_grid(cape_noHGTQS_noSHAL)
cape_noHGTQS_noUVmix = cut_sub_grid(cape_noHGTQS_noUVmix)

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