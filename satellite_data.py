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
import pyhdf
#import GOES
from calculating_latlon import calculate_degrees
from cut_subgrid_module import cut_sub_grid
import os


'''Import data (brightness temperatures are L2 products'''
data = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\satellite\OR_ABI-L2-CMIPF-M6C13_G16_s20230580000206_e20230580009525_c20230580009596.nc')
data = data.isel(x = slice(2700,4000), y = slice(1000,2700)) #the plotting does not work if there are nan lat and lon values. So cut the x and y. 
#cll = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cll_noUVmix_noHGTQS\cll_noUVmix_noHGTQS_cut\cll_his_BES_HA43h22tg3_clim_noHGTQS_noUVmix_1hr_202001010000-202001102300.nc')


'''Regridding the area to lat and lon in degrees instead of fix radians'''
#this is a function that is made in the calculate_degrees module
lat, lon = calculate_degrees(data)

'''Assign the dataset the lat and lon values'''

ds = xr.Dataset({'C13_Brightness_Temp': (('y', 'x'), data.CMI.values)},
                  {'lat': (['y','x'],lat), 'lon': (['y','x'],lon), 'time': ( data.t.values)}, 
                  {'units': ('Kelvin'), 
                   'long_name': ('Brightness Temperatures for cloud and moisture')})



cmap = "Blues_r"
ds = cut_sub_grid(ds)
plt.figure()
ds.C13_Brightness_Temp.plot(x = 'lon', y = 'lat', cmap = cmap)
ds = ds.expand_dims(dim = 'time')

time = np.array(ds.time.values)

#ds = ds.to_array(dim = 'time')

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
            #scene = ds.isel(time = i)
            plt.figure()
            plt.imshow(scene, cmap=cmap)
            break
        T_cl_min = 280
        T_cl_max = 290
        
        q_cutoff = 25
        T_cutoff = 285
        
        T_q = np.percentile(scene, q_cutoff)
        
        cloud_mask = np.zeros(scene.shape,dtype=int)
        cloud_mask[(scene < T_cl_max) & (scene > T_cl_min)] = 1
        plt.imshow(cloud_mask, cmap=cmap)
        
        if np.percentile(scene, q_cutoff) < T_cutoff:
            print('Reject due to high clouds')
        #else:
            #print('Less than 25% high clouds')
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
                df_metrics.iloc[i, df_metrics.columns.get_loc(metrics[j])] = fn_metric(object_labels) # its supposed to save it in the dataframe but it doesnt


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
                    df_metrics.iloc[i, df_metrics.columns.get_loc(metrics[j])] = fn_metric(scene.values) 
                    
                    
            # Store after each scene
            #save_directory = 'C:\\Users\\LENOVO\\Desktop\\TU_Delft\\thesis\\' #metric_results\\'
            #number = str(i) 
            #filename = 'df_metrics_noHGTQS_noUVmix' + number + '.h5'
            #save_path = os.path.join(save_directory,filename)
            #df_metrics.to_hdf(save_path, 'cloudmetrics', mode='w')            
df_metrics



