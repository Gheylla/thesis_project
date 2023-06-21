# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 09:24:09 2023

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
import os 
import datetime as dt
from cut_subgrid_module import cut_sub_grid

#%% 
'''Import all data'''
evs_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\evs_data\evspsbl_noHGTQS_cut.nc', combine ='by_coords')
evs_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\evs_data\evspsbl_noHGTQS_noSHAL_cut.nc', combine ='by_coords')
evs_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\evs_data\evspsbl_noHGTQS_noUVmix_cut.nc', combine ='by_coords')

hfls_noHGTQS = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\hfls_data\hfls_eva_noHGTQS_cut.nc', combine ='by_coords')
hfls_noHGTQS_noSHAL = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\hfls_data\hfls_eva_noHGTQS_noSHAL_cut.nc', combine ='by_coords')
hfls_noHGTQS_noUVmix = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\hfls_data\hfls_eva_noHGTQS_noUVmix_cut.nc', combine ='by_coords')
#%%
'''Get area'''
x = evs_noHGTQS.x.shape[0]
y = evs_noHGTQS.y.shape[0]

res = 2.5 #km
sizex = x * res
sizey = y * res

area = sizex * sizey #in km2

#%%
'''Get the mean of the domain'''
evs_noHGTQS = evs_noHGTQS.mean(dim = ('x', 'y'))
evs_noHGTQS_noSHAL = evs_noHGTQS_noSHAL.mean(dim = ('x', 'y'))
evs_noHGTQS_noUVmix = evs_noHGTQS_noUVmix.mean(dim = ('x', 'y'))

hfls_noHGTQS = hfls_noHGTQS.mean(dim = ('x', 'y'))
hfls_noHGTQS_noSHAL = hfls_noHGTQS_noSHAL.mean(dim = ('x', 'y'))
hfls_noHGTQS_noUVmix = hfls_noHGTQS_noUVmix.mean(dim = ('x', 'y'))

#%%
'''Determine time passed'''
#feb 1st is at index 744
#everytime step is another hour added
#The units are joules per square metre (J m-2). To convert to watts per square metre (W m-2), the accumulated values should be divided by the accumulation period expressed in seconds. ECMWF

jan = np.linspace(1,743,743)
feb = np.linspace(744,1440,697)

accum_time_sec_jan = []
hour_jan = 1

for m in range(len(jan)):
    time_sec_j = hour_jan * 60 * 60
    hour_jan += 1
    accum_time_sec_jan.append(time_sec_j)
    
    
    
accum_time_sec_feb = []
hour_feb = 0

for k in range(len(feb)):
    time_sec_f = hour_feb * 60 * 60
    hour_feb += 1
    accum_time_sec_feb.append(time_sec_f)
    

#%%
'''Get W/m2 for latent heat'''
latent_jan_nohgtqs = hfls_noHGTQS.hfls_eva.values[0:743] / accum_time_sec_jan
latent_feb_nohgtqs = hfls_noHGTQS.hfls_eva.values[743:1440] / accum_time_sec_feb


latent_jan_nohgtqs_noshal = hfls_noHGTQS_noSHAL.hfls_eva.values[0:743] / accum_time_sec_jan
latent_feb_nohgtqs_noshal = hfls_noHGTQS_noSHAL.hfls_eva.values[743:1440] / accum_time_sec_feb


latent_jan_nohgtqs_nouvmix = hfls_noHGTQS_noUVmix.hfls_eva.values[0:743] / accum_time_sec_jan
latent_feb_nohgtqs_nouvmix = hfls_noHGTQS_noUVmix.hfls_eva.values[743:1440] / accum_time_sec_feb



#%%
total_lat_nohgtqs = np.concatenate((latent_jan_nohgtqs, latent_feb_nohgtqs))
total_lat_nohgtqs_noshal = np.concatenate((latent_jan_nohgtqs_noshal, latent_feb_nohgtqs_noshal))
total_lat_nohgtqs_nouvmix = np.concatenate((latent_jan_nohgtqs_nouvmix, latent_feb_nohgtqs_nouvmix))


#%%
#total_lat_nohgtqs[0] = 'NaN'
total_lat_nohgtqs_noshal[0] = 'NaN'
total_lat_nohgtqs_nouvmix[0] ='NaN'

total_lat_nohgtqs[743] = 'NaN'
total_lat_nohgtqs_noshal[743] = 'NaN'
total_lat_nohgtqs_nouvmix[743] ='NaN'

#%%
xr_latent_noHGTQS = xr.Dataset({'latent_heat': (('time'), total_lat_nohgtqs )},
                  {'time': ( hfls_noHGTQS.time.values)}, 
                  {'units': ('W/m2'), 
                   'long_name': ('Latent heat flux HARMONIE noHGTQS W/m2')})

xr_latent_noHGTQS_noSHAL = xr.Dataset({'latent_heat': (('time'), total_lat_nohgtqs_noshal )},
                  {'time': ( hfls_noHGTQS_noSHAL.time.values)}, 
                  {'units': ('W/m2'), 
                   'long_name': ('Latent heat flux HARMONIE noHGTQS noSHAL W/m2')})

xr_latent_noHGTQS_noUVmix = xr.Dataset({'latent_heat': (('time'), total_lat_nohgtqs_nouvmix )},
                  {'time': ( hfls_noHGTQS_noUVmix.time.values)}, 
                  {'units': ('W/m2'), 
                   'long_name': ('Latent heat flux HARMONIE noHGTQS noUVmix W/m2')})

#%%
#da.dropna(dim="Y", how="any")

xr_latent_noHGTQS = xr_latent_noHGTQS.dropna(dim = 'time' )
xr_latent_noHGTQS_noSHAL = xr_latent_noHGTQS_noSHAL.dropna(dim = 'time' )
xr_latent_noHGTQS_noUVmix = xr_latent_noHGTQS_noUVmix.dropna(dim = 'time' )


#%%

lat_nohgtqs_diur = xr_latent_noHGTQS.latent_heat.groupby(xr_latent_noHGTQS.time.dt.hour).mean()
lat_nohgtqs_noshal_diur = xr_latent_noHGTQS_noSHAL.latent_heat.groupby(xr_latent_noHGTQS_noSHAL.time.dt.hour).mean()
lat_nohgtqs_nouvmix_diur = xr_latent_noHGTQS_noUVmix.latent_heat.groupby(xr_latent_noHGTQS_noUVmix.time.dt.hour).mean()


#%%
hours = np.linspace(00,23, 24)
lt = [20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

fig, ax = plt.subplots(layout='constrained', figsize = (15,5))
#ax.plot(lat_nohgtqs_diur.hour.values, lat_nohgtqs_diur, color = 'red', label = 'HARMONIE noHGTQS')
ax.plot(lat_nohgtqs_diur.hour.values, lat_nohgtqs_noshal_diur, color = 'blue', label = 'HARMONIE noHGTQS noSHAL')
#ax.plot(lat_nohgtqs_diur.hour.values, lat_nohgtqs_nouvmix_diur, color = 'green', label = 'HARMONIE noHGTQS noUVmix')
 
ax.grid()
ax.set_xticks(hours)
ax.set_xlabel('UTC Time [hr]')
secax = ax.secondary_xaxis('top')
secax.set_xticks(hours, lt)
secax.set_xlabel('Local Time [hr]')

plt.title('Composite diurnal cycle of latent heat flux')
plt.xlabel('Time [UTC]')
plt.ylabel('Accumulated latent heat flux [W/m2]')
plt.legend()
#plt.ylim([205, 230])




#%%
hours = np.linspace(00,23, 24)
lt = [20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

fig, ax = plt.subplots(layout='constrained', figsize = (15,5))
ax.plot(hfls_noHGTQS.time.values, total_lat_nohgtqs, color = 'red')
ax.plot(hfls_noHGTQS_noSHAL.time.values, total_lat_nohgtqs_noshal, color = 'blue')
ax.plot(hfls_noHGTQS_noUVmix.time.values, total_lat_nohgtqs_nouvmix, color = 'green')

ax.grid()
ax.set_xticks(hours)
ax.set_xlabel('UTC Time [hr]')
secax = ax.secondary_xaxis('top')
secax.set_xticks(hours, lt)
secax.set_xlabel('Local Time [hr]')

plt.title('Accumulated latent heat flux')
plt.xlabel('Time [UTC]')
plt.ylabel('Accumulated latent heat flux [W/m2]')



#%%
'''Remove accumulated latent heat and make it hourly'''

lat_nohgtqs_hrly = []

for i in range(len(latent_jan_nohgtqs)-1):
    if i == 0:
        lat_hr_nohgtqs = latent_jan_nohgtqs[0]
        lat_nohgtqs_hrly.append(lat_hr_nohgtqs)
    else:
        lat_hr_nohgtqs = latent_jan_nohgtqs[i + 1] - latent_jan_nohgtqs[i]
        lat_nohgtqs_hrly.append(lat_hr_nohgtqs)
        
plt.figure()
plt.plot(lat_nohgtqs_hrly)

#this is not correct. 

#%%

#evs_noHGTQS
hrly_eva_nohgtqs = []
for f in range(len(evs_noHGTQS.time.values)):
    if f == 0 :
        hrly_eva_nohgtqs.append(evs_noHGTQS.evspsbl.values[0])
    if f == 743:
        hrly_eva_nohgtqs.append(evs_noHGTQS.evspsbl.values[743])
    else:
        hrly_eva_nohgtqs.append(evs_noHGTQS.evspsbl.values[f + 1] - evs_noHGTQS.evspsbl.values[f])
        
        


























# =============================================================================
# 
# #%%
# '''Plot the accumulated evaporation and latent heat'''
# 
# plt.figure(figsize = (15,10))
# plt.plot(evs_noHGTQS.time.values, evs_noHGTQS.values, label = 'HARMONIE noHGTQS', color = 'red')
# plt.plot(evs_noHGTQS_noSHAL.time.values, evs_noHGTQS_noSHAL.values, label = 'HARMONIE noHGTQS noSHAL', color = 'blue')
# plt.plot(evs_noHGTQS_noUVmix.time.values, evs_noHGTQS_noUVmix.values, label = 'HARMONIE noHGTQS noUVmix', color = 'green')
# plt.grid()
# plt.title('Accumulated Surface Water Evaporation Flux')
# plt.xlabel('Time [hr]')
# plt.ylabel('Surface Water Evaporation Flux [kg m-2]')
# plt.legend()
# 
# plt.figure(figsize = (15,10))
# plt.plot(hfls_noHGTQS.time.values, hfls_noHGTQS.hfls_eva.values, label = 'HARMONIE noHGTQS', color = 'red' )
# plt.plot(hfls_noHGTQS_noSHAL.time.values, hfls_noHGTQS_noSHAL.hfls_eva.values, label = 'HARMONIE noHGTQS noSHAL', color = 'blue')
# plt.plot(hfls_noHGTQS_noUVmix.time.values, hfls_noHGTQS_noUVmix.hfls_eva.values, label = 'HARMONIE noHGTQS noUVmix', color = 'green')
# plt.grid()
# plt.title('Accumulated Surface Upward Latent Heat Flux Due To Evavaporation')
# plt.xlabel('Time [hr]')
# plt.ylabel('Surface Upward Latent Heat Flux [J m-2]')
# plt.legend()
# 
# 
# #%% 
# '''Get the not accumulated '''
# 
# 
# min_evs_noHGTQS = np.min(evs_noHGTQS.values)
# min_evs_noHGTQS_noSHAL = np.min(evs_noHGTQS_noSHAL.values)
# min_evs_noHGTQS_noUVmix = np.min(evs_noHGTQS_noUVmix.values)
# 
# loc_min_noHGTQS = np.where(evs_noHGTQS.values == min_evs_noHGTQS )
# loc_min_noHGTQS_noSHAL = np.where(evs_noHGTQS_noSHAL.values == min_evs_noHGTQS_noSHAL )
# loc_min_noHGTQS_noSHAL = np.where(evs_noHGTQS_noUVmix.values == min_evs_noHGTQS_noUVmix )
# 
# 
# hr_evs_noHGTQS = np.zeros(len(evs_noHGTQS.time.values))
# hr_evs_noHGTQS_noSHAL = np.zeros(len(evs_noHGTQS_noSHAL.time.values))
# hr_evs_noHGTQS_noUVmix = np.zeros(len(evs_noHGTQS_noUVmix.time.values))
# 
# for i in range((len(evs_noHGTQS.time.values)) - 1):
#     print(i)
#     if i == 0: 
#         hr_evs_noHGTQS[i] = evs_noHGTQS.values[i]
#         hr_evs_noHGTQS_noSHAL[i] = evs_noHGTQS_noSHAL.values[i]
#         hr_evs_noHGTQS_noUVmix[i] = evs_noHGTQS_noUVmix.values[i]
#     else:
#         evs_noHGTQS_before = evs_noHGTQS.values[i - 1]
#         evs_noHGTQS_noSHAL_before = evs_noHGTQS_noSHAL.values[i - 1]
#         evs_noHGTQS_noUVmix_before = evs_noHGTQS_noUVmix.values[i - 1]
#         evs_noHGTQS_now = evs_noHGTQS.values[i]
#         evs_noHGTQS_noSHAL_now = evs_noHGTQS_noSHAL.values[i]
#         evs_noHGTQS_noUVmix_now = evs_noHGTQS_noUVmix.values[i]
#         hr_evs_noHGTQS[i] =  evs_noHGTQS_now -  evs_noHGTQS_before
#         hr_evs_noHGTQS_noSHAL[i] =  evs_noHGTQS_noSHAL_now -  evs_noHGTQS_noSHAL_before
#         hr_evs_noHGTQS_noUVmix[i] =  evs_noHGTQS_noUVmix_now -  evs_noHGTQS_noUVmix_before
#         if i == loc_min_noHGTQS[0][0]:
#             hr_evs_noHGTQS[i] = evs_noHGTQS.values[i]
#             hr_evs_noHGTQS_noSHAL[i] = evs_noHGTQS_noSHAL.values[i]
#             hr_evs_noHGTQS_noUVmix[i] = evs_noHGTQS_noUVmix.values[i]
#             
#             
# #%%
# hr_hfls_noHGTQS = np.zeros(len(hfls_noHGTQS.time.values))
# hr_hfls_noHGTQS_noSHAL = np.zeros(len(hfls_noHGTQS_noSHAL.time.values))
# hr_hfls_noHGTQS_noUVmix = np.zeros(len(hfls_noHGTQS_noUVmix.time.values))
# 
# for i in range((len(hfls_noHGTQS.time.values)) - 1):
#     print(i)
#     if i == 0: 
#         hr_hfls_noHGTQS[i] = hfls_noHGTQS.hfls_eva.values[i]
#         hr_hfls_noHGTQS_noSHAL[i] = hfls_noHGTQS_noSHAL.hfls_eva.values[i]
#         hr_hfls_noHGTQS_noUVmix[i] = hfls_noHGTQS_noUVmix.hfls_eva.values[i]
#     else:
#         hfls_noHGTQS_before = hfls_noHGTQS.hfls_eva.values[i - 1]
#         hfls_noHGTQS_noSHAL_before = hfls_noHGTQS_noSHAL.hfls_eva.values[i - 1]
#         hfls_noHGTQS_noUVmix_before = hfls_noHGTQS_noUVmix.hfls_eva.values[i - 1]
#         hfls_noHGTQS_now = hfls_noHGTQS.hfls_eva.values[i]
#         hfls_noHGTQS_noSHAL_now = hfls_noHGTQS_noSHAL.hfls_eva.values[i]
#         hfls_noHGTQS_noUVmix_now = hfls_noHGTQS_noUVmix.hfls_eva.values[i]
#         hr_hfls_noHGTQS[i] =  hfls_noHGTQS_now -  hfls_noHGTQS_before
#         hr_hfls_noHGTQS_noSHAL[i] = hfls_noHGTQS_noSHAL_now - hfls_noHGTQS_noSHAL_before
#         hr_hfls_noHGTQS_noUVmix[i] = hfls_noHGTQS_noUVmix_now - hfls_noHGTQS_noUVmix_before
#         if i == loc_min_noHGTQS[0][0]:
#             hr_hfls_noHGTQS[i] = hfls_noHGTQS.hfls_eva.values[i]
#             hr_hfls_noHGTQS_noSHAL[i] = hfls_noHGTQS_noSHAL.hfls_eva.values[i]
#             hr_hfls_noHGTQS_noUVmix[i] = hfls_noHGTQS_noUVmix.hfls_eva.values[i]
# 
# 
# #%%
# '''Compute evaporation for the area'''
# 
# Le = 2.25 * 10 **6 #in J/kg
# area_m2 = area * 1000000
# 
# hr_evs_noHGTQS_complarea = hr_evs_noHGTQS * 1000000 #make it per km2
# hr_evs_noHGTQS_complarea = hr_evs_noHGTQS_complarea * area #remove the per km2 you are left with kg/hr
# hr_evs_noHGTQS_complarea = hr_evs_noHGTQS_complarea * Le #you are left with J/hr 
# hr_evs_noHGTQS_complarea = hr_evs_noHGTQS_complarea / area_m2 #this is j/(hr m2)
# 
# 
# 
# 
# #%%
# '''Plot not accumulated'''
# plt.figure(figsize = (15,10))
# plt.plot(evs_noHGTQS.time.values, hr_evs_noHGTQS, label = 'HARMONIE noHGTQS ' )          
# plt.plot(evs_noHGTQS.time.values, hr_evs_noHGTQS_noSHAL, label = 'HARMONIE noHGTQS noSHAL' )   
# plt.plot(evs_noHGTQS.time.values, hr_evs_noHGTQS_noUVmix, label = 'HARMONIE noHGTQS noUVmix')     
# plt.grid()
# plt.title('Surface Water Evaporation Flux')
# plt.xlabel('Time [hr]')
# plt.ylabel('Surface Water Evaporation Flux [kg m-2]')
# plt.legend()
# 
# plt.figure(figsize = (15,10))
# plt.plot(hfls_noHGTQS.time.values, hr_hfls_noHGTQS, label = 'HARMONIE noHGTQS')          
# plt.plot(hfls_noHGTQS_noSHAL.time.values, hr_hfls_noHGTQS_noSHAL, label = 'HARMONIE noHGTQS noSHAL')   
# plt.plot(hfls_noHGTQS_noUVmix.time.values, hr_hfls_noHGTQS_noUVmix, label = 'HARMONIE noHGTQS noUVmix')     
# plt.grid()
# plt.title('Surface Upward Latent Heat Flux Due To Evavaporation')
# plt.xlabel('Time [hr]')
# plt.ylabel('Surface Upward Latent Heat Flux [J m-2]')
# plt.legend()
# 
# =============================================================================
