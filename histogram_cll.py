# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 15:17:30 2023

@author: User
"""

import numpy
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import scipy as sc
from scipy.signal import argrelmin 
import datetime as dt
import inspect 
import cloudmetrics
import plotly.graph_objects as go
import seaborn as sn

# =============================================================================
# fERA43 = xr.open_mfdataset(r'C:\Users\User\Desktop\TU DELFT\Master_2022-2023\Thesis\Test_Data\cll_fERA_HARM43_cut\cll_*', combine = 'by_coords')
# noHGTQS_ECUME43 = xr.open_mfdataset(r'C:\Users\User\Desktop\TU DELFT\Master_2022-2023\Thesis\Test_Data\cll_noHGTQS_noSHAL43_cut\cll_*', combine = 'by_coords')
# noHGTQS43 = xr.open_mfdataset(r'C:\Users\User\Desktop\TU DELFT\Master_2022-2023\Thesis\Test_Data\cll_noHGTQS43_cut\cll_*', combine = 'by_coords')
# noSHAL43 = xr.open_mfdataset(r'C:\Users\User\Desktop\TU DELFT\Master_2022-2023\Thesis\Test_Data\cll_noSHAL43_cut\cll_*', combine = 'by_coords')
# noSHAL_noHGTQS = xr.open_mfdataset(r'C:\Users\User\Desktop\TU DELFT\Master_2022-2023\Thesis\Test_Data\cll_noHGTQS_noSHAL43_cut\cll_*', combine = 'by_coords')
# =============================================================================
noUVmix_noHGTQS = xr.open_mfdataset(r'C:\Users\User\Desktop\TU DELFT\Master_2022-2023\Thesis\Test_Data\cll_noUVmix_noHGTQS_cut\cll_*', combine = 'by_coords')

# =============================================================================
# fERA_hist, fERA_bins = np.histogram(fERA43.cll, bins = 10, range = [0,1] )
# noHGTQS_ECUME43_hist, noHGTQS_ECUME43_bins = np.histogram(noHGTQS_ECUME43.cll, bins = 10, range = [0,1] )
# noHGTQS43_hist, noHGTQS43_bins = np.histogram(noHGTQS43.cll, bins = 10, range = [0,1] )
# noSHAL43_hist, noSHAL43_bins = np.histogram(noSHAL43.cll, bins = 10, range = [0,1] )
# noSHAL_noHGTQS_hist, noSHAL_noHGTQS_bins = np.histogram(noSHAL_noHGTQS.cll, bins = 10, range = [0,1] )
# =============================================================================
noUVmix_noHGTQS_hist, noUVmix_noHGTQS_bins = np.histogram(noUVmix_noHGTQS.cll, bins = 10, range = [0,1] )



# =============================================================================
# 
# plt.figure()
# plt.title('Histogram HARMONIE 43 fERA Climate run ')
# plt.ylabel('Count [-]')
# plt.xlabel('Cloud Cover [-]')
# plt.grid()
# plt.bar(fERA_bins[:-1],fERA_hist,width=0.1, align = 'edge', edgecolor = 'blue')
# plt.ylim(0, 3e8)
# plt.show()
# 
# plt.figure()
# plt.title('Histogram HARMONIE 43 noHGTQS ECUME ')
# plt.ylabel('Count [-]')
# plt.xlabel('Cloud Cover [-]')
# plt.grid()
# plt.bar(noHGTQS_ECUME43_bins[:-1],noHGTQS_ECUME43_hist,width=0.1, align = 'edge',  edgecolor = 'blue')
# plt.ylim(0, 3e8)
# plt.show()
# 
# plt.figure()
# plt.title('Histogram HARMONIE 43 noHGTQS')
# plt.ylabel('Count [-]')
# plt.xlabel('Cloud Cover [-]')
# plt.grid()
# plt.bar(noHGTQS43_bins[:-1],noHGTQS43_hist,width=0.1, align = 'edge',  edgecolor = 'blue')
# plt.ylim(0, 3e8)
# plt.show()
# 
# plt.figure()
# plt.title('Histogram HARMONIE 43 noSHAL')
# plt.ylabel('Count [-]')
# plt.xlabel('Cloud Cover [-]')
# plt.grid()
# plt.bar(noSHAL43_bins[:-1],noSHAL43_hist,width=0.1, align = 'edge',  edgecolor = 'blue')
# plt.ylim(0, 3e8)
# plt.show()
# =============================================================================



# =============================================================================
# plt.figure()
# plt.title('Histogram HARMONIE 43 noSHAL noHGTQS')
# plt.ylabel('Count [-]')
# plt.xlabel('Cloud Cover [-]')
# plt.grid()
# plt.bar(noSHAL_noHGTQS_bins[:-1],noSHAL_noHGTQS_hist,width=0.1, align = 'edge',  edgecolor = 'blue')
# plt.ylim(0, 3e8)
# plt.show()
# =============================================================================

plt.figure()
plt.title('Histogram HARMONIE 43 noUVmix noHGTQS')
plt.ylabel('Count [-]')
plt.xlabel('Cloud Cover [-]')
plt.grid()
plt.bar(noUVmix_noHGTQS_bins[:-1],noUVmix_noHGTQS_hist,width=0.1, align = 'edge',  edgecolor = 'blue')
plt.ylim(0, 3e8)
plt.show()