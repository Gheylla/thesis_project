# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 14:45:25 2023

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
from cut_subgrid_module import cut_sub_grid
import os


directory = r'C:\Users\LENOVO\Desktop\sat'
entries = os.listdir(directory)

'''Compute metrics'''
metrics = [
           'cloud_fraction',
           'fractal_dimension',
           'open_sky',
           'cop',
           'iorg',
           'scai',
           'max_length_scale',
           'mean_eccentricity',
           'mean_length_scale',
           'mean_perimeter_length',
           'num_objects',
           'orientation',
           'spectral_length_moment',
           'spectral_anisotropy',
           'spectral_slope',
           'woi1',
           'woi2',
           'woi3',
           'mean',
           'var',
          ]

x = np.linspace(0,((len(entries)) - 1), ((len(entries))))
data_pd = pd.DataFrame(index = x, columns = metrics)
time_sat = []

for entry in entries:
    python_indices  = [index for (index, item) in enumerate(entries) if item == entry]
    print(python_indices)
    file_location_name = os.path.join(directory,entry)

    cmap = "Blues_r"
    '''Import data (brightness temperatures are L2 products'''
    
    ds = xr.open_mfdataset(file_location_name)
# =============================================================================
#     #plt.figure(figsize = (10,10))
#     #plt.contourf(data.x, data.y, data.CMI, cmap = cmap)
# 
#     '''Regridding the area to lat and lon in degrees instead of fix radians'''
#     #this is a function that is made in the calculate_degrees module
#     lat, lon = calculate_degrees(data)
#     
#     '''Assign the dataset the lat and lon values'''
#     
#     ds = xr.Dataset({'C13_Brightness_Temp': (('y', 'x'), data.CMI.values)},
#                       {'lat': (['y','x'],lat), 'lon': (['y','x'],lon), 'time': ( data.t.values)}, 
#                       {'units': ('Kelvin'), 
#                        'long_name': ('Brightness Temperatures for cloud and moisture')})
#     
#     #plt.figure(figsize = (10,10))
#     ds = ds.isel(x = slice(3000,4000), y = slice(1000,4000)) #after converting to lat and lon there are nan values of lat and lon andthis makes it impossible to plot. so this is why we cut
#     #ds.C13_Brightness_Temp.plot(x = 'lon', y = 'lat', cmap = cmap)
#     #plt.contourf(ds.lon.values, ds.lat.values, ds.C13_Brightness_Temp.values, cmap = cmap)
#     
#     latbegin = 10
#     latend = 20
#     lonbegin = -48
#     lonend = -61
#     
#     lat_begin_arr = []
#     lat_end_arr = []
#     for val in ds.lat.values[:,0]: #for latitude you need to change the rows because if you stay on one row but change cols the latitude stays the same. 
#         distance_lat_begin = np.abs(val - latbegin)
#         lat_begin_arr.append(distance_lat_begin)
#         distance_lat_end = np.abs(val - latend)
#         lat_end_arr.append(distance_lat_end)
#      
#     index_closest_begin_lat = np.argwhere(lat_begin_arr == np.min(lat_begin_arr))       
#     index_closest_end_lat = np.argwhere(lat_end_arr == np.min(lat_end_arr))           
#             
#     ds = ds.isel(y = slice(index_closest_end_lat[0,0], (index_closest_begin_lat[0,0]) ))
#     
#     lon_begin_arr = []
#     lon_end_arr = []
#     for val in ds.lon.values[0,:]: #for lonigtude you need to change columns, if you stay on the same column the longitude remains the same.
#         distance_lon_begin = np.abs(val - lonbegin)
#         lon_begin_arr.append(distance_lon_begin)
#         distance_lon_end = np.abs(val - lonend)
#         lon_end_arr.append(distance_lon_end)
#      
#     index_closest_begin_lon = np.argwhere(lon_begin_arr == np.min(lon_begin_arr))       
#     index_closest_end_lon = np.argwhere(lon_end_arr == np.min(lon_end_arr))   
#     
#     ds = ds.isel(x = slice(index_closest_end_lon[0,0], (index_closest_begin_lon[0,0]) ))
#     
#     #ds = cut_sub_grid(ds)
#     plt.figure(figsize = (10,10))
#     ds.C13_Brightness_Temp.plot(x = 'lon', y = 'lat', cmap = cmap)
#     
#     y = 13.1939
#     x = -59.5432
#     plt.plot(x, y, marker="o", markersize=20, markeredgecolor="red", markerfacecolor="red")
#     plt.text((x + 0.2), (y + 0.2), 'BCO')
#     
#     
#     ds = ds.expand_dims(dim = 'time')
#     
#     time = np.array(ds.time.values)
#     
#     x = len(ds.x)
#     y = len(ds.y)
#             
#     if x == y: 
#         print('Field is a square')
#     else: 
#         imin = np.argmin(ds.lat.shape)
#         imin = ds.lat.shape[0]
#         div = imin 
#         if (div % 2) == 0 :
#             ds = ds.isel(x = slice(0, (imin)), y = slice(0, (imin)))
#         else:
#             ds = ds.isel(x = slice(0, (imin-1)), y = slice(0, (imin-1)))
#     
#     save_directory = 'C:\\Users\\LENOVO\\Desktop\\TU_Delft\\thesis\\data\\cut_sat\\'
#     filename = entry
#     filename = entry + '.nc'
#     save_path = os.path.join(save_directory,filename)        
#     ds.to_netcdf(save_path)
#     
#     #plt.figure(figsize = (10,10))
#     #ds.C13_Brightness_Temp.plot(x = 'lon', y = 'lat', cmap = cmap)
#     
#     #y = 13.1939
#     #x = -59.5432
#     #plt.plot(x, y, marker="o", markersize=20, markeredgecolor="red", markerfacecolor="red")
#     #plt.text((x + 0.2), (y + 0.2), 'BCO')
# 
# 
# =============================================================================


    
    '''Compute metrics'''
    metrics = [
               'cloud_fraction',
               'fractal_dimension',
               'open_sky',
               'cop',
               'iorg',
               'scai',
               'max_length_scale',
               'mean_eccentricity',
               'mean_length_scale',
               'mean_perimeter_length',
               'num_objects',
               'orientation',
               'spectral_length_moment',
               'spectral_anisotropy',
               'spectral_slope',
               'woi1',
               'woi2',
               'woi3',
               'mean',
               'var',
              ]
    
    # These are the metrics you can choose from
    available_mask_metrics = dict(inspect.getmembers(cloudmetrics.mask, inspect.isfunction))
    available_object_metrics = dict(inspect.getmembers(cloudmetrics.objects, inspect.isfunction))
    available_scalar_metrics = dict(inspect.getmembers(cloudmetrics.scalar, inspect.isfunction))
    
    load_df = False

    
    if load_df:
        df_metrics = pd.read_hdf(r'C:\Users\User\Desktop\TU DELFT\Master_2022-2023\Thesis\Scripts\df_metrics.h5')
    else:
        df_metrics = pd.DataFrame(index =([ds.time.values]), columns=metrics)
        
        for i in range(ds.time.size):
            #print('it works')
            print('Processing ', 'scene ', i+1, '/', ds.time.size)
            while True:
                # Cut for time and subgrid
                cmap = "RdBu_r"
                scene = ds.C13_Brightness_Temp.isel(time = i)
                time_sat.append(scene.time.values)
                #scene = ds.isel(time = i)
                plt.figure()
                plt.imshow(scene, cmap=cmap)
                break
            T_cl_min = 277
            T_cl_max = 290
            
            q_cutoff = 25
            T_cutoff = 285
            
            T_q = np.percentile(scene, q_cutoff)
            
            cloud_mask = np.zeros(scene.shape,dtype=int)
            cloud_mask[(scene < T_cl_max) & (scene > T_cl_min)] = 1
            plt.imshow(cloud_mask, cmap=cmap)
            
            if np.percentile(scene, q_cutoff) < T_cutoff:
                print('Reject due to high clouds')
                continue
    
            computed_object_labels = False
            computed_spectra = False
            for j in range(len(metrics)):
                # Cloud object metrics
                if metrics[j] in available_object_metrics.keys():
                    print('Computing', metrics[j])
    
                    # Compute object labels if not done yet
                    if not computed_object_labels:
                        object_labels = cloudmetrics.objects.label(cloud_mask)
                        computed_object_labels = True
    
                    # Compute metric
                    fn_metric = available_object_metrics[metrics[j]]                
                    df_metrics.iloc[i, df_metrics.columns.get_loc(metrics[j])] = fn_metric(object_labels) 
    
    
                elif metrics[j] in available_mask_metrics.keys():
                    print('Computing', metrics[j])
                    fn_metric = available_mask_metrics[metrics[j]]
    
                    # Open sky exception - just take the mean open sky area (second function output)
                    if 'open_sky' in metrics[j]:
                        _, df_metrics.iloc[i, df_metrics.columns.get_loc(metrics[j])] = fn_metric(cloud_mask)
                    else:
                        df_metrics.iloc[i, df_metrics.columns.get_loc(metrics[j])] = fn_metric(cloud_mask)
    
                # Cloud scalar metrics
                elif metrics[j] in available_scalar_metrics.keys():
                    print('Computing', metrics[j])
                    fn_metric = available_scalar_metrics[metrics[j]]
    
    
                    # Spectral metrics exception
                    if 'spectral' in metrics[j]:
    
                        # Compute spectra if not done yet
                        if not computed_spectra:
                            wavenumbers, psd_1d_radial, psd_1d_azimuthal = cloudmetrics.scalar.compute_spectra(scene.values)
                            computed_spectra = True
    
                        # Compute metrics
                        if 'anisotropy' in metrics[j]:
                            df_metrics.iloc[i, df_metrics.columns.get_loc(metrics[j])] = fn_metric(psd_1d_azimuthal)
                        else:
                            df_metrics.iloc[i, df_metrics.columns.get_loc(metrics[j])] = fn_metric(wavenumbers, psd_1d_radial)
    
                    # All other scalar metrics computed normally
                    else:
                        df_metrics.iloc[i, df_metrics.columns.get_loc(metrics[j])] = fn_metric(scene.values) #issue is here
                # Store after each scene
                save_directory = 'C:\\Users\\LENOVO\\Desktop\\TU_Delft\\thesis\\metric_results\\Observations\\' #metric_results\\'
                number = str(python_indices) 
                data_pd.iloc[python_indices] = df_metrics
            print(df_metrics)


df_metrics
data_pd['Time'] = time_sat

save_directory = 'C:\\Users\\LENOVO\\Desktop\\TU_Delft\\thesis\\metric_results\\'
filename = 'metrics_obs_293_286_test.h5'
save_path = os.path.join(save_directory,filename)
data_pd.to_hdf(save_path, 'cloudmetrics', mode='w')




