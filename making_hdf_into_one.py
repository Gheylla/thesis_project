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



directory = r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\277_290'
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


tot_data_pandas.to_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\277_290\complete_observations_277290.h5', key='df', mode='w')  


