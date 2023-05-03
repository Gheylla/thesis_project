# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 17:06:18 2023

@author: LENOVO
"""

import eurec4a
from cut_subgrid_module import cut_sub_grid
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import cloudmetrics
import inspect
from skimage.transform import resize



from intake import open_catalog
cat = open_catalog("https://raw.githubusercontent.com/eurec4a/eurec4a-intake/master/catalog.yml")

cat = eurec4a.get_intake_catalog()
data = xr.concat([cat.satellites.GOES16['latlongrid'](date=d).to_dask().chunk() for d in pd.date_range('2020-01-12','2020-01-20')], dim='time')
#data = xr.concat([cat.satellites.GOES16['latlongrid'](date=d).to_dask().chunk() for d in pd.date_range('2020-01-22','2020-01-30')], dim='time')


#data = cat.satellites.GOES16['latlongrid'].to_dask()


vmin = 270
vmax = 300
cmap = "RdBu_r"


#data = data.C13


#from harmonie we used these lat from 10.049367 to 20.32508 and lon from -61.785774 to -49.850334

lat_begin = 10
lat_end = 20

min_lat_beg = []
min_lat_end = []

for i in range(len(data.lat.values)):
    diff_lat_beg = np.abs(data.lat.values[i] - lat_begin)
    min_lat_beg.append(diff_lat_beg)
    diff_lat_end = np.abs(data.lat.values[i] - lat_end)
    min_lat_end.append(diff_lat_end)
    
    
min_lat_beg_val = np.min(min_lat_beg)
min_lat_end_val = np.min(min_lat_end) 

loc_lat_beg = np.where(min_lat_beg == min_lat_beg_val)[0][0]
loc_lat_end = np.where(min_lat_end == min_lat_end_val)[0][0]      
   


data = data.isel(lat = slice(loc_lat_end, loc_lat_beg))

lon_begin = -48
lon_end = -61

min_lon_beg = []
min_lon_end = []

for i in range(len(data.lon.values)):
    diff_lon_beg = np.abs(data.lon.values[i] - lon_begin)
    min_lon_beg.append(diff_lon_beg)
    diff_lon_end = np.abs(data.lon.values[i] - lon_end)
    min_lon_end.append(diff_lon_end)
    
    
min_lon_beg_val = np.min(min_lon_beg)
min_lon_end_val = np.min(min_lon_end) 

loc_lon_beg = np.where(min_lon_beg == min_lon_beg_val)[0][0]
loc_lon_end = np.where(min_lon_end == min_lon_end_val)[0][0]      
   
data = data.isel(lon = slice(loc_lon_end, loc_lon_beg))


x = len(data.lat)
y = len(data.lon)
        
if x == y: 
    print('Field is a square')
else: 
    imin = min(x,y)
    if (imin %2) == 0 :
        data = data.isel(lat= slice(0, (imin)), lon = slice(0, (imin)))
    else:
        data = data.isel(lat = slice(0, (imin-1)), lon = slice(0, (imin-1)))
        
#this results in a square with even numbers
  
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
    df_metrics = pd.DataFrame(index=data.time.values, columns=metrics)
    
    for i in range(len(data.time.values)):
        print('Processing ', data.time[i].values, 'scene ', i+1, '/', len(data.time.values))
        while True:
            # Cut for time and subgrid
            cmap = "RdBu_r"
            scene = data.C13.isel(time = i)
            try:
                # Resize
                imin = np.argmin(scene.shape)
                scene = resize(scene, (scene.shape[imin],scene.shape[imin]), anti_aliasing=True)
                #assert(np.isnan(scene).any() == False)
            except:
                print('Unable to access scene, retrying...')
                continue
            break
            plt.imshow(scene, cmap=cmap, vmin=vmin, vmax=vmax)
            plt.show()
            break

        T_cl_min = 280
        T_cl_max = 290
        
        q_cutoff = 25
        T_cutoff = 285
        
        #T_q = np.percentile(scene, q_cutoff)
        
        cmap = "RdBu_r"
        cloud_mask = np.zeros(scene.shape,dtype=int)
        cloud_mask[(scene < T_cl_max) & (scene > T_cl_min)] = 1
        plt.imshow(cloud_mask, cmap=cmap)
        
        scene = resize(scene, (scene.shape[0],scene.shape[1]), anti_aliasing=True)
       
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
                        wavenumbers, psd_1d_radial, psd_1d_azimuthal = cloudmetrics.scalar.compute_spectra(scene)
                        computed_spectra = True

                    # Compute metrics
                    if 'anisotropy' in metrics[j]:
                        df_metrics.iloc[i, df_metrics.columns.get_loc(metrics[j])] = fn_metric(psd_1d_azimuthal)
                    else:
                        df_metrics.iloc[i, df_metrics.columns.get_loc(metrics[j])] = fn_metric(wavenumbers, psd_1d_radial)

                # All other scalar metrics computed normally
                else:
                    df_metrics.iloc[i, df_metrics.columns.get_loc(metrics[j])] = fn_metric(scene)
            # Store after each scene
            #df_metrics.to_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\df_metrics_noHGTQS_noUVmix.h5', 'cloudmetrics', mode='w')
df_metrics





