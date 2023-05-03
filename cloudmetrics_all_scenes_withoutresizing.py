# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 14:44:38 2023

@author: Gheylla Liberia
testting out how to make cloud masks
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
from skimage.transform import resize
from PIL import Image 
import netCDF4 as nc

'''Step 1: Import data'''
print('Step 1: Import data')
'''For this specific case we are using low cloud cover fields from HARMONIE outputs'''

path_dataset = r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cll_noUVmix_noHGTQS\cll_noUVmix_noHGTQS_cut\cll_his*'
dataset = xr.open_mfdataset(path_dataset, combine = 'by_coords', engine = 'netcdf4')

#path_height_43 = r'C:\Users\User\Desktop\TU DELFT\Master_2022-2023\Additional_Thesis\Data\z_43_fullheights.csv'
#height_43 = pd.read_csv(path_height_43)

#path_height_40 = r'C:\Users\User\Desktop\TU DELFT\Master_2022-2023\Additional_Thesis\Data\z_40_fullheights.csv'
#height_40 = pd.read_csv(path_height_40)


'''The following are all of the metrics we would like to compute'''

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
print('Step 1 Completed ')

'''Step 3: Set up loop to determine metrics for all scenes'''
print('Step 3: Set-up loop for all scenes')

'''All days and hours will be looped  through in order to calculate the metrics, these will then be saved in a dataset. 
ach hour is set as a scene and this will be used to calculate the metrics'''


#set load_df to false in order to make the dataframe with all the metrics from scratch otherwise it will just load the h5 file
load_df = False

if load_df:
    df_metrics = pd.read_hdf(r'C:\Users\User\Desktop\TU DELFT\Master_2022-2023\Thesis\Scripts\df_metrics.h5')
else:
    df_metrics = pd.DataFrame(index=dataset.time.values, columns=metrics)
    
    for i in range(len(dataset.time.values)):
        print('Processing ', dataset.time[i].values, 'scene ', i+1, '/', len(dataset.time.values))
        while True:
            # Cut for time and subgrid
            cmap = "RdBu_r"
            scene = dataset.isel(time = i)
            
            x = len(scene.x)
            y = len(scene.y)
            
            if x == y: 
                print('Field is a square')
            else: 
                imin = np.argmin(scene.cll.shape)
                imin = scene.cll.shape[0]
                scene = scene.isel(x = slice(0, imin), y = slice(0, imin))
                
            
            test_sub = scene.cll
            break

        cll_threshold = 0.5
        cloud_mask = np.zeros(scene.cll.shape,dtype=int)
        cloud_mask[(scene.cll > cll_threshold)] = 1
        
        lat = scene.lat.values
        lon = scene.lon.values
        time = scene.time.values

        mask_da = xr.DataArray(cloud_mask, coords = (scene.y.values, scene.x.values,), dims = ['y','x'])
        mask_da.name = 'cll_mask'
        mask_da.attrs['units'] = '-'

        mask_da_ds = xr.Dataset({'cll_mask': (['y','x'], mask_da.data)}, 
                coords = {'x': scene.x.values, 'y': scene.y.values, 'time': scene.time.values})
        
        
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
                        wavenumbers, psd_1d_radial, psd_1d_azimuthal = cloudmetrics.scalar.compute_spectra(scene.cll.values)
                        computed_spectra = True

                    # Compute metrics
                    if 'anisotropy' in metrics[j]:
                        df_metrics.iloc[i, df_metrics.columns.get_loc(metrics[j])] = fn_metric(psd_1d_azimuthal)
                    else:
                        df_metrics.iloc[i, df_metrics.columns.get_loc(metrics[j])] = fn_metric(wavenumbers, psd_1d_radial)

                # All other scalar metrics computed normally
                else:
                    df_metrics.iloc[i, df_metrics.columns.get_loc(metrics[j])] = fn_metric(scene.cll.values)
            # Store after each scene
            #df_metrics.to_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\df_metrics_noHGTQS_noUVmix.h5', 'cloudmetrics', mode='w')
df_metrics


