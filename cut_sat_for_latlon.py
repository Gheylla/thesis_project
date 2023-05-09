# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 16:35:55 2023

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
import seaborn as sn
#import pyhdf
#import GOES
from calculating_latlon import calculate_degrees
#from cut_subgrid_module import cut_sub_grid
import os


directory = r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\satellite'
entries = os.listdir(directory)


for entry in entries:
    python_indices  = [index for (index, item) in enumerate(entries) if item == entry]
    print(python_indices)
    file_location_name = os.path.join(directory,entry)

    '''Import data (brightness temperatures are L2 products'''
    
    data = xr.open_mfdataset(file_location_name)


    '''Regridding the area to lat and lon in degrees instead of fix radians'''
    #this is a function that is made in the calculate_degrees module
    lat, lon = calculate_degrees(data)
    
    '''Assign the dataset the lat and lon values'''
    
    ds = xr.Dataset({'C13_Brightness_Temp': (('y', 'x'), data.CMI.values)},
                      {'lat': (['y','x'],lat), 'lon': (['y','x'],lon), 'time': ( data.t.values)}, 
                      {'units': ('Kelvin'), 
                       'long_name': ('Brightness Temperatures for cloud and moisture')})
    
    #plt.figure(figsize = (10,10))
    ds = ds.isel(x = slice(3000,4000), y = slice(1000,4000)) #after converting to lat and lon there are nan values of lat and lon andthis makes it impossible to plot. so this is why we cut
    
    latbegin = 10
    latend = 20
    lonbegin = -48
    lonend = -61
    
    lat_begin_arr = []
    lat_end_arr = []
    for val in ds.lat.values[:,0]: #for latitude you need to change the rows because if you stay on one row but change cols the latitude stays the same. 
        distance_lat_begin = np.abs(val - latbegin)
        lat_begin_arr.append(distance_lat_begin)
        distance_lat_end = np.abs(val - latend)
        lat_end_arr.append(distance_lat_end)
     
    index_closest_begin_lat = np.argwhere(lat_begin_arr == np.min(lat_begin_arr))       
    index_closest_end_lat = np.argwhere(lat_end_arr == np.min(lat_end_arr))           
            
    ds = ds.isel(y = slice(index_closest_end_lat[0,0], (index_closest_begin_lat[0,0]) ))
    
    lon_begin_arr = []
    lon_end_arr = []
    for val in ds.lon.values[0,:]: #for lonigtude you need to change columns, if you stay on the same column the longitude remains the same.
        distance_lon_begin = np.abs(val - lonbegin)
        lon_begin_arr.append(distance_lon_begin)
        distance_lon_end = np.abs(val - lonend)
        lon_end_arr.append(distance_lon_end)
     
    index_closest_begin_lon = np.argwhere(lon_begin_arr == np.min(lon_begin_arr))       
    index_closest_end_lon = np.argwhere(lon_end_arr == np.min(lon_end_arr))   
    
    ds = ds.isel(x = slice(index_closest_end_lon[0,0], (index_closest_begin_lon[0,0]) ))
    
    #ds = cut_sub_grid(ds)
    ds = ds.expand_dims(dim = 'time')
    
    time = np.array(ds.time.values)
    
    x = len(ds.x)
    y = len(ds.y)
            
    if x == y: 
        print('Field is a square')
    else: 
        imin = np.argmin(ds.lat.shape)
        imin = ds.lat.shape[0]
        div = imin 
        if (div % 2) == 0 :
            ds = ds.isel(x = slice(0, (imin)), y = slice(0, (imin)))
        else:
            ds = ds.isel(x = slice(0, (imin-1)), y = slice(0, (imin-1)))
    
    save_directory = 'C:\\Users\\LENOVO\\Desktop\\TU_Delft\\thesis\\data\\cut_sat\\'
    filename = entry
    #filename = entry + '.nc'
    save_path = os.path.join(save_directory,filename)        
    ds.to_netcdf(save_path)
    


