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
#from calculating_latlon import calculate_degrees
#from cut_subgrid_module import cut_sub_grid
import os
from skimage.transform import resize
import time


directory = r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\cut_sat'
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
    print(python_indices) #in order to know where we are in the list
    file_location_name = os.path.join(directory,entry)

    cmap = "Blues_r"
    '''Import data (brightness temperatures are L2 products'''
    
    ds = xr.open_mfdataset(file_location_name)
    
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
                scene = ds.C13_Brightness_Temp.isel(time = i) #scene only has the brightness temperatures
                time_sat.append(scene.time.values)
                try:
                    # Resize
                    imin = np.argmin(scene.shape)
                    scene = resize(scene, (scene.shape[imin],scene.shape[imin]), anti_aliasing=True) #we are resizing but its already a square so nothing is supposed to change. 
                    assert(np.isnan(scene).any() == False)
                except:
                    save_directory = 'C:\\Users\\LENOVO\\Desktop\\TU_Delft\\thesis\\metric_results\\'
                    filename = 'metrics_obs_2000.h5'
                    save_path = os.path.join(save_directory,filename)
                    data_pd.to_hdf(save_path, 'cloudmetrics', mode='w')
                    print(df_metrics)
                    data_pd['Time'] = time_sat
                    #print('Unable to access scene, retrying...')
                    #time.sleep(5)
                    continue
                
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
            plt.figure()
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
                            wavenumbers, psd_1d_radial, psd_1d_azimuthal = cloudmetrics.scalar.compute_spectra(scene)
                            computed_spectra = True
    
                        # Compute metrics
                        if 'anisotropy' in metrics[j]:
                            df_metrics.iloc[i, df_metrics.columns.get_loc(metrics[j])] = fn_metric(psd_1d_azimuthal)
                        else:
                            df_metrics.iloc[i, df_metrics.columns.get_loc(metrics[j])] = fn_metric(wavenumbers, psd_1d_radial)
    
                    # All other scalar metrics computed normally
                    else:
                        df_metrics.iloc[i, df_metrics.columns.get_loc(metrics[j])] = fn_metric(scene) #issue is here
                # Store after each scene
                save_directory = 'C:\\Users\\LENOVO\\Desktop\\TU_Delft\\thesis\\metric_results\\Observations\\' #metric_results\\'
                number = str(python_indices) 
                data_pd.iloc[python_indices] = df_metrics
                
                


    save_directory = 'C:\\Users\\LENOVO\\Desktop\\TU_Delft\\thesis\\metric_results\\'
    filename = 'metrics_obs_2000.h5'
    save_path = os.path.join(save_directory,filename)
    data_pd.to_hdf(save_path, 'cloudmetrics', mode='w')
    print(df_metrics)


data_pd['Time'] = time_sat




