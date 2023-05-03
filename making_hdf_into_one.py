# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 15:03:56 2023

@author: LENOVO
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
import datetime as dt
import h5py
import os

output = 'metric_observations.h5'
combined_data = []



directory = r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\Observations'
entries = os.listdir(directory)


for entry in entries:
    file_location_name = os.path.join(directory,entry)
    input_file = file_location_name
    data = pd.read_hdf(input_file)
    combined_data.append(data)
    
total_data = np.concatenate(combined_data, axis=0)
tot_data_pandas = pd.DataFrame(total_data, columns =data.keys())
tot_data_pandas['Time'] = pd.to_datetime(tot_data_pandas['Time'])
tot_data_pandas = tot_data_pandas.sort_values(tot_data_pandas.columns[20], ascending = True)
tot_data_pandas.set_index(tot_data_pandas['Time'], inplace=True)


tot_data_pandas.to_hdf('complete_observations.h5', key='df', mode='w')  

#kaminda no tin data ta paso tbt tin muchu high clouds anto nos a reject e data nan ey. Kaminda e 
#linanan ta kore recht bai topa e otro eynan falta data ku nos no a run ahinda. 

cloud_frac_grouphr = tot_data_pandas.groupby([tot_data_pandas.index.hour]).cloud_fraction.mean()
iorg_grouphr = tot_data_pandas.groupby([tot_data_pandas.index.hour]).iorg.mean()

plt.figure(figsize = (20,10))
plt.plot(tot_data_pandas.Time[1:], tot_data_pandas.cloud_fraction[1:])
plt.grid()
plt.title('Timeseries of cloud fraction')
plt.xlabel('Time')
plt.ylabel('CF [-]')



df_metrics_noHGTQS43 = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\HARMONIE\df_metrics_noHGTQS.h5')
df_metrics_noHGTQS_noSHAL = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\HARMONIE\df_metrics_noHGTQS_noSHAL.h5')
df_metrics_noHGTQS_noUVmix = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\HARMONIE\df_metrics_noHGTQS_noUVmix.h5')

noHGTQS_diurcf = df_metrics_noHGTQS43.cloud_fraction.groupby(df_metrics_noHGTQS43.index.hour).mean()
noHGTQS_noSHAL_diurcf = df_metrics_noHGTQS_noSHAL.cloud_fraction.groupby(df_metrics_noHGTQS_noSHAL.index.hour).mean()


plt.figure(figsize = (10,5))
plt.plot(cloud_frac_grouphr* 100, label = 'Observations', color = 'black', linewidth=2.5)
plt.plot(noHGTQS_diurcf * 100, label = 'HARMONIE noHGTQS')
plt.plot(noHGTQS_noSHAL_diurcf * 100, label = 'HARMONIE noHGTQS noSHAL')
plt.legend()
plt.grid()

hours = np.linspace(0,23,24)
plt.xticks(hours)

#only thing missing is checking time of this series. 

# =============================================================================
#     print(input_file)
#     with h5py.File(input_file, 'r') as f:
#         print(f)
#         data = f[entry][:]
#         #combined_data.append(data)
# =============================================================================
