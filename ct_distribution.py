# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 14:57:50 2023

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
ct_nohgtqs = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\ct_data\ct_noHGTQS.nc')
ct_nohgtqs_noshal =  xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\ct_data\ct_noHGTQS_noSHAL.nc')
ct_nohgtqs_nouvmix =  xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\ct_data\ct_noHGTQS_noUVmix.nc')

#%%
'''Cut the data to desired grid'''
ct_nohgtqs = cut_sub_grid(ct_nohgtqs)
ct_nohgtqs_noshal = cut_sub_grid(ct_nohgtqs_noshal)
ct_nohgtqs_nouvmix = cut_sub_grid(ct_nohgtqs_nouvmix)

#%%
'''Make a distribution of the cloud tops'''
#we expect no shal to have more cloud tops
ct_nohgtqs_arr = np.array(ct_nohgtqs.ct.values).reshape(-1)
ct_nohgtqs_noshal_arr = np.array(ct_nohgtqs_noshal.ct.values).reshape(-1)
ct_nohgtqs_nouvmix_arr = np.array(ct_nohgtqs_nouvmix.ct.values).reshape(-1)

#%%
'''Remove the negative values'''

for i in range(len(ct_nohgtqs_nouvmix_arr)):
    print(i) #total is 
    if ct_nohgtqs_arr[i] <0: 
        ct_nohgtqs_arr[i] = 'NaN'
    if ct_nohgtqs_noshal_arr[i] < 0:
        ct_nohgtqs_noshal_arr[i] = 'NaN'
    if ct_nohgtqs_nouvmix_arr[i] <0: 
        ct_nohgtqs_nouvmix_arr[i] = 'NaN'

#%%
plt.figure()
plt.hist(ct_nohgtqs_arr ,   bins =30 , color = 'red', label = 'HARMONIE noHGTQS', density = True, histtype = 'step', align = 'mid')
plt.hist(ct_nohgtqs_noshal_arr ,   bins = 30, color = 'blue', label = 'HARMONIE noHGTQS noSHAL', density = True, histtype = 'step', align = 'mid')
plt.hist(ct_nohgtqs_nouvmix_arr ,   bins = 30, color = 'green', label = 'HARMONIE noHGTQS noUVmix', density = True, histtype = 'step', align = 'mid')
plt.legend()
plt.title("Histogram of cloud top")
plt.xlabel('CT [m]')
plt.ylabel('Count')
plt.grid()
plt.show()


# =============================================================================
# 
# plt.figure()
# plt.hist(ct_nohgtqs_noshal_arr ,   bins = 20, color = 'blue', label = 'HARMONIE noHGTQS noSHAL')
# plt.legend()
# plt.title("Histogram of cloud top")
# plt.xlabel('CT [m]')
# plt.ylabel('Count')
# plt.grid()
# plt.show()
# 
# plt.figure()
# plt.hist(ct_nohgtqs_nouvmix_arr ,   bins = 20, color = 'green', label = 'HARMONIE noHGTQS noUVmix')
# plt.legend()
# plt.title("Histogram of cloud top")
# plt.xlabel('CT [m]')
# plt.ylabel('Count')
# plt.grid()
# plt.show()
# 
# =============================================================================
