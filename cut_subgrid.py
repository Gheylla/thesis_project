# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 11:12:51 2023

@author: Gheylla Liberia 
Cut the datset for the subgrid required
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
import os 

latbegin = 10
latend = 20
lonbegin = -48
lonend = -61


path_of_file = r'C:\Users\User\Desktop\TU DELFT\Master_2022-2023\Thesis\Test_Data\cll_noUVmix_noHGTQS'

for filename in os.listdir(path_of_file):
        if '.nc' in filename:
                print(filename)
                f = os.path.join(path_of_file,filename)
                data = xr.open_dataset(f, engine = 'netcdf4', decode_times = True)
                save_dir = r'C:\Users\User\Desktop\TU DELFT\Master_2022-2023\Thesis\Test_Data\cll_noUVmix_noHGTQS_cut'
                save_path = os.path.join(save_dir,filename)
                
                distance_lat_begin = np.abs(data.lat - latbegin) 
                distance_lat_end = np.abs(data.lat - latend) 

                distance_lon_begin = np.abs(data.lon - lonbegin) 
                distance_lon_end = np.abs(data.lon - lonend)

                index_closest_begin_lon = np.argwhere(distance_lon_begin.values == np.min(distance_lon_begin.values))
                index_closest_end_lon = np.argwhere(distance_lon_end.values == np.min(distance_lon_end.values))
               
                data = data.isel(x = slice(index_closest_end_lon[0,1], index_closest_begin_lon[0,1] ), y = slice( index_closest_begin_lon[0,0], index_closest_end_lon[0,0]))
                print(data.lat.values)
                print(data.lon.values)
                print(len(data.x.values))
                print(len(data.y.values))
                
                data.to_netcdf(path = save_path)
