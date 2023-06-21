# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 17:58:40 2023

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
import plotly.graph_objects as go
import os 
import datetime as dt

#%%
'''Import data'''
conv_dry_nohgtqs = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\thflx_data\thflx_conv_dry_noHGTQS.nc')
conv_dry_nohgtqs_noshal = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\thflx_data\thflx_conv_dry_noHGTQS_noSHAL.nc')
conv_moist_nohgtqs = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\thflx_data\thflx_conv_moist_noHGTQS.nc')
conv_moist_nohgtqs_noshal = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\thflx_data\thflx_conv_moist_noHGTQS_noSHAL.nc')
turb_nohgtqs = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\thflx_data\thflx_turb_noHGTQS.nc')
turb_nohgtqs_noshal = xr.open_mfdataset(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\thflx_data\thflx_turb_noHGTQS_noSHAL.nc')

heights = pd.read_csv(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\data\z_43_fullheights.csv')


#%%
'''Get the mean of the area'''
conv_dry_nohgtqs_mean = conv_dry_nohgtqs.mean(dim = ('x', 'y'))
conv_dry_nohgtqs_noshal_mean = conv_dry_nohgtqs_noshal.mean(dim = ('x', 'y'))
conv_moist_nohgtqs_mean = conv_moist_nohgtqs.mean(dim = ('x', 'y'))
conv_moist_nohgtqs_noshal_mean = conv_moist_nohgtqs_noshal.mean(dim = ('x', 'y'))
turb_nohgtqs_mean = turb_nohgtqs.mean(dim = ('x', 'y'))
turb_nohgtqs_noshal_mean = turb_nohgtqs_noshal.mean(dim = ('x', 'y'))

#%%
'''Load all values to make computation faster'''
values_conv_dry_nohgtqs = conv_dry_nohgtqs_mean.Thlflx_conv_dry.values
values_conv_dry_nohgtqs_noshal = conv_dry_nohgtqs_noshal_mean.Thlflx_conv_dry.values
values_conv_moist_nohgtqs = conv_moist_nohgtqs_mean.Thlflx_conv_mois.values
values_conv_moist_nohgtqs_noshal =conv_moist_nohgtqs_noshal_mean.Thlflx_conv_mois.values
values_turb_nohgtqs = turb_nohgtqs_mean.Thlflx_turb.values
values_turb_nohgtqs_noshal = turb_nohgtqs_noshal_mean.Thlflx_turb.values


#%%
'''Remove the accumulation'''
shape_convdry_nohgtqs = conv_dry_nohgtqs_mean.Thlflx_conv_dry.shape
shape_convdry_nohgtqsnoshal = conv_dry_nohgtqs_noshal_mean.Thlflx_conv_dry.shape
shape_convmoist_nohgtqs = conv_moist_nohgtqs_mean.Thlflx_conv_mois.shape
shape_convmoist_nohgtqs_noshal = conv_moist_nohgtqs_noshal_mean.Thlflx_conv_mois.shape
shape_turb_nohgtqs = turb_nohgtqs_mean.Thlflx_turb.shape
shape_turb_nohgtqs_noshal = turb_nohgtqs_noshal_mean.Thlflx_turb.shape


nonaccum_convdry_nohgtqs = np.zeros((shape_convdry_nohgtqs[0], shape_convdry_nohgtqs[1]))
nonaccum_convdry_nohgtqs_noshal = np.zeros((shape_convdry_nohgtqsnoshal[0], shape_convdry_nohgtqsnoshal[1]))
nonaccum_convmoist_nohgtqs = np.zeros((shape_convmoist_nohgtqs[0], shape_convmoist_nohgtqs[1]))
nonaccum_convmoist_nohgtqs_noshal = np.zeros((shape_convmoist_nohgtqs_noshal[0], shape_convmoist_nohgtqs_noshal[1]))
nonaccum_turb_nohgtqs = np.zeros((shape_turb_nohgtqs[0], shape_turb_nohgtqs[1]))
nonaccum_turb_nohgtqs_noshal = np.zeros((shape_turb_nohgtqs_noshal[0], shape_turb_nohgtqs_noshal[1]))

#%% done these work perfectly to get the non accumulated data
'''Make non accumulated data'''

for i in range(len(conv_dry_nohgtqs_mean.lev.values)):
    for ii in range(len(conv_dry_nohgtqs_mean.time.values)):
        if ii == 0 or ii == 744: #the accumulation starts again at i = 744
            nonaccum_convdry_nohgtqs[ii,i] = values_conv_dry_nohgtqs[ii,i]
        else: 
            nonaccum_convdry_nohgtqs[ii,i] = values_conv_dry_nohgtqs[ii,i] - values_conv_dry_nohgtqs[ii -1 ,i]
            
            

xr_conv_dry_nohgtqs = xr.Dataset({'Thlflx_conv_dry': (('time', 'lev'), nonaccum_convdry_nohgtqs )},
                  {'time': ( conv_dry_nohgtqs_mean.time.values)}, 
                  {'units': ('-'), 
                   'long_name': ('Theta l flux non accum')})


            
#%% done
for k in range(len(conv_dry_nohgtqs_noshal_mean.lev.values)):
    for kk in range(len(conv_dry_nohgtqs_noshal_mean.time.values)):
        if kk== 0 or kk == 744:
            nonaccum_convdry_nohgtqs_noshal[kk,k] = values_conv_dry_nohgtqs_noshal[kk,k]
        else: 
            nonaccum_convdry_nohgtqs_noshal[kk,k] = values_conv_dry_nohgtqs_noshal[kk,k] - values_conv_dry_nohgtqs_noshal[kk -1 ,k]
            
            
xr_conv_dry_nohgtqs_noshal = xr.Dataset({'Thlflx_conv_dry': (('time', 'lev'), nonaccum_convdry_nohgtqs_noshal )},
                  {'time': ( conv_dry_nohgtqs_noshal_mean.time.values)}, 
                  {'units': ('-'), 
                   'long_name': ('Theta l flux non accum')})
            
#%% done
for l in range(len(conv_moist_nohgtqs_mean.lev.values)):
    for ll in range(len(conv_moist_nohgtqs_mean.time.values)):
        if ll == 0 or ll == 744:
            nonaccum_convmoist_nohgtqs[ll,l] = values_conv_moist_nohgtqs[ll,l]
        else: 
            nonaccum_convmoist_nohgtqs[ll,l] = values_conv_moist_nohgtqs[ll,l] - values_conv_moist_nohgtqs[ll -1 ,l]
            
            
xr_conv_moist_nohgtqs = xr.Dataset({'Thlflx_conv_moist': (('time', 'lev'), nonaccum_convmoist_nohgtqs )},
                  {'time': ( conv_moist_nohgtqs_mean.time.values)}, 
                  {'units': ('-'), 
                   'long_name': ('Theta l flux non accum')})
            
#%% done
for f in range(len(conv_moist_nohgtqs_noshal_mean.lev.values)):
    for ff in range(len(conv_moist_nohgtqs_noshal_mean.time.values)):
        if ff == 0 or ff == 744:
            nonaccum_convmoist_nohgtqs_noshal[ff,f] = values_conv_moist_nohgtqs_noshal[ff,f]
        else: 
            nonaccum_convmoist_nohgtqs_noshal[ff,f] = values_conv_moist_nohgtqs_noshal[ff,f] - values_conv_moist_nohgtqs_noshal[ff -1 ,f]
            
xr_conv_moist_nohgtqs_noshal = xr.Dataset({'Thlflx_conv_moist': (('time', 'lev'), nonaccum_convmoist_nohgtqs_noshal )},
                  {'time': ( conv_moist_nohgtqs_noshal_mean.time.values)}, 
                  {'units': ('-'), 
                   'long_name': ('Theta l flux non accum')})
            
            
#%% done             
for m in range(len(turb_nohgtqs_mean.lev.values)):
    for mm in range(len(turb_nohgtqs_mean.time.values)):
        if mm == 0 or mm == 744:
            nonaccum_turb_nohgtqs[mm,m] = values_turb_nohgtqs[mm,m]
        else: 
            nonaccum_turb_nohgtqs[mm,m] = values_turb_nohgtqs[mm,m] - values_turb_nohgtqs[mm -1 ,m]
            
xr_turb_nohgtqs = xr.Dataset({'Thlflx_turb': (('time', 'lev'), nonaccum_turb_nohgtqs )},
                  {'time': ( turb_nohgtqs_mean.time.values)}, 
                  {'units': ('-'), 
                   'long_name': ('Theta l flux non accum')})
            
    
    
#%%
for n in range(len(turb_nohgtqs_noshal_mean.lev.values)):
    for nn in range(len(turb_nohgtqs_noshal_mean.time.values)):
        if nn == 0 or nn == 744:
            nonaccum_turb_nohgtqs_noshal[nn,n] = values_turb_nohgtqs_noshal[nn,n]
        else: 
            nonaccum_turb_nohgtqs_noshal[nn,n] = values_turb_nohgtqs_noshal[nn,n] - values_turb_nohgtqs_noshal[nn -1 ,n]
            
            
xr_turb_nohgtqs_noshal = xr.Dataset({'Thlflx_turb': (('time', 'lev'), nonaccum_turb_nohgtqs_noshal )},
                  {'time': ( turb_nohgtqs_noshal_mean.time.values)}, 
                  {'units': ('-'), 
                   'long_name': ('Theta l flux non accum')})
            
#%%    
'''Compute the composite diurnal cycle of all'''
xr_conv_dry_nohgtqs_diur = xr_conv_dry_nohgtqs.groupby(xr_conv_dry_nohgtqs.time.dt.hour).mean()
xr_conv_dry_nohgtqs_noshal_diur = xr_conv_dry_nohgtqs_noshal.groupby(xr_conv_dry_nohgtqs_noshal.time.dt.hour).mean()
xr_conv_moist_nohgtqs_diur = xr_conv_moist_nohgtqs.groupby(xr_conv_moist_nohgtqs.time.dt.hour).mean()
xr_conv_moist_nohgtqs_noshal_diur = xr_conv_moist_nohgtqs_noshal.groupby(xr_conv_moist_nohgtqs_noshal.time.dt.hour).mean()
xr_turb_nohgtqs_diur = xr_turb_nohgtqs.groupby(xr_turb_nohgtqs.time.dt.hour).mean()
xr_turb_nohgtqs_noshal_diur = xr_turb_nohgtqs_noshal.groupby(xr_turb_nohgtqs_noshal.time.dt.hour).mean()

#%%
'''Transpose the diurnal data'''
#we need to transpose the data in order to get the hour values on the x-axis
xr_conv_dry_nohgtqs_diur_trans = xr_conv_dry_nohgtqs_diur.transpose()
xr_conv_dry_nohgtqs_noshal_diur_trans = xr_conv_dry_nohgtqs_noshal_diur.transpose()
xr_conv_moist_nohgtqs_diur_trans = xr_conv_moist_nohgtqs_diur.transpose()
xr_conv_moist_nohgtqs_noshal_diur_trans = xr_conv_moist_nohgtqs_noshal_diur.transpose()
xr_turb_nohgtqs_diur_trans = xr_turb_nohgtqs_diur.transpose()
xr_turb_nohgtqs_noshal_diur_trans = xr_turb_nohgtqs_noshal_diur.transpose()

#%%
'''Plot the diurnal cycle'''
#we need to flip the data because the highest level is the one closest to the ground
hours = np.linspace(00,23, 24)
lt = [20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

fig, ax = plt.subplots(layout='constrained', figsize = (15,5))
im = ax.contourf(xr_conv_dry_nohgtqs_diur_trans.hour.values, heights['Full heights'], np.flip(xr_conv_dry_nohgtqs_diur_trans.Thlflx_conv_dry.values, axis = 0))
plt.ylim([0,3000])
plt.title('Composite diurnal cycle of $\Theta_l$ flux dry convection for HARMONIE noHGTQS')
ax.grid()
ax.set_xticks(hours)
ax.set_xlabel('UTC Time [hr]')
secax = ax.secondary_xaxis('top')
secax.set_xticks(hours, lt)
secax.set_xlabel('Local Time [hr]')
plt.ylabel('Height [m]')
fig.colorbar(im, label = '$\Theta_l$ flux [K m-1 s-1]')
plt.show()



fig1, ax1 = plt.subplots(layout='constrained', figsize = (15,5))
im1 = ax1.contourf(xr_conv_dry_nohgtqs_noshal_diur_trans.hour.values, heights['Full heights'], np.flip(xr_conv_dry_nohgtqs_noshal_diur_trans.Thlflx_conv_dry.values, axis = 0))
plt.ylim([0,3000])
plt.title('Composite diurnal cycle of $\Theta_l$ flux dry convection for HARMONIE noHGTQS noSHAL')
ax1.grid()
ax1.set_xticks(hours)
ax1.set_xlabel('UTC Time [hr]')
secax1 = ax1.secondary_xaxis('top')
secax1.set_xticks(hours, lt)
secax1.set_xlabel('Local Time [hr]')
plt.ylabel('Height [m]')
fig1.colorbar(im1, label = '$\Theta_l$ flux [K m-1 s-1]')
plt.show()

fig2, ax2 = plt.subplots(layout='constrained', figsize = (15,5))
im2 = ax2.contourf(xr_conv_moist_nohgtqs_diur_trans.hour.values, heights['Full heights'], np.flip(xr_conv_moist_nohgtqs_diur_trans.Thlflx_conv_moist.values, axis = 0))
plt.ylim([0,3000])
plt.title('Composite diurnal cycle of $\Theta_l$ flux moist convection for HARMONIE noHGTQS')
ax2.grid()
ax2.set_xticks(hours)
ax2.set_xlabel('UTC Time [hr]')
secax2 = ax2.secondary_xaxis('top')
secax2.set_xticks(hours, lt)
secax2.set_xlabel('Local Time [hr]')
plt.ylabel('Height [m]')
fig2.colorbar(im2, label = '$\Theta_l$ flux [K m-1 s-1]')
plt.show()


fig3, ax3 = plt.subplots(layout='constrained', figsize = (15,5))
im3 = ax3.contourf(xr_conv_moist_nohgtqs_noshal_diur_trans.hour.values, heights['Full heights'], np.flip(xr_conv_moist_nohgtqs_noshal_diur_trans.Thlflx_conv_moist.values, axis = 0))
plt.ylim([0,3000])
plt.title('Composite diurnal cycle of $\Theta_l$ flux moist convection for HARMONIE noHGTQS noSHAL')
ax3.grid()
ax3.set_xticks(hours)
ax3.set_xlabel('UTC Time [hr]')
secax3 = ax3.secondary_xaxis('top')
secax3.set_xticks(hours, lt)
secax3.set_xlabel('Local Time [hr]')
plt.ylabel('Height [m]')
fig3.colorbar(im3, label = '$\Theta_l$ flux [K m-1 s-1]')
plt.show()


#%%
levels = np.linspace(-60,60,31)


fig4, ax4 = plt.subplots(layout='constrained', figsize = (15,5))
im4 = ax4.contourf(xr_turb_nohgtqs_diur_trans.hour.values, heights['Full heights'], np.flip(xr_turb_nohgtqs_diur_trans.Thlflx_turb.values, axis = 0), levels = levels, vmin = -60, vmax = 60)
plt.ylim([0,3000])
plt.title('Composite diurnal cycle of $\Theta_l$ turbulent for HARMONIE noHGTQS')
ax4.grid()
ax4.set_xticks(hours)
ax4.set_xlabel('UTC Time [hr]')
secax4 = ax4.secondary_xaxis('top')
secax4.set_xticks(hours, lt)
secax4.set_xlabel('Local Time [hr]')
plt.ylabel('Height [m]')
fig4.colorbar(im4, label = '$\Theta_l$ flux [K m-1 s-1]')
plt.show()



fig5, ax5 = plt.subplots(layout='constrained', figsize = (15,5))
im5 = ax5.contourf(xr_turb_nohgtqs_noshal_diur_trans.hour.values, heights['Full heights'], np.flip(xr_turb_nohgtqs_noshal_diur_trans.Thlflx_turb.values, axis = 0), levels = levels, vmin = -60, vmax = 60)
plt.ylim([0,3000])
plt.title('Composite diurnal cycle of $\Theta_l$ turbulent flux for HARMONIE noHGTQS noSHAL')
ax5.grid()
ax5.set_xticks(hours)
ax5.set_xlabel('UTC Time [hr]')
secax5 = ax5.secondary_xaxis('top')
secax5.set_xticks(hours, lt)
secax5.set_xlabel('Local Time [hr]')
plt.ylabel('Height [m]')
fig5.colorbar(im5, label = '$\Theta_l$ flux [K m-1 s-1]')
plt.show()


#%%
'''Transpose the data to plot complete'''
            
nonaccum_convdry_nohgtqs_trans = nonaccum_convdry_nohgtqs.transpose()
nonaccum_convdry_nohgtqs_noshal_trans = nonaccum_convdry_nohgtqs_noshal.transpose()
nonaccum_convmoist_nohgtqs_trans = nonaccum_convmoist_nohgtqs.transpose()
nonaccum_convmoist_nohgtqs_noshal_trans = nonaccum_convmoist_nohgtqs_noshal.transpose()
nonaccum_turb_nohgtqs_trans = nonaccum_turb_nohgtqs.transpose()
nonaccum_turb_nohgtqs_noshal_trans = nonaccum_turb_nohgtqs_noshal.transpose()


#%%
'''Plot the timeseries of the non accumulated data'''
plt.figure(figsize = (20,10))
plt.contourf(conv_dry_nohgtqs_mean.time.values, heights['Full heights'], np.flip(nonaccum_convdry_nohgtqs_trans, axis = 0))
plt.ylim([0,3000])
plt.title('Timeseries of $\Theta_l$ flux dry convection for HARMONIE noHGTQS')
plt.ylabel('Height [m]')
plt.xlabel('Time [UTC]')
plt.colorbar(label = '$\Theta_l$ flux [K m-1 s-1]')

# =============================================================================
# plt.figure(figsize = (15,10))
# plt.contourf(conv_dry_nohgtqs_noshal_mean.time.values, heights['Full heights'], np.flip(nonaccum_convdry_nohgtqs_noshal_trans, axis = 0))
# plt.ylim([0,3000])
# plt.colorbar()
# =============================================================================

plt.figure(figsize = (20,10))
plt.contourf(conv_dry_nohgtqs_mean.time.values, heights['Full heights'], np.flip(nonaccum_convmoist_nohgtqs_trans, axis = 0))
plt.ylim([0,3000])
plt.colorbar(label = '$\Theta_l$ flux [K m-1 s-1]')
plt.ylabel('Height [m]')
plt.xlabel('Time [UTC]')
plt.title('Timeseries of $\Theta_l$ flux moist convection for HARMONIE noHGTQS')

# =============================================================================
# plt.figure(figsize = (15,10))
# plt.contourf(conv_dry_nohgtqs_mean.time.values, heights['Full heights'], np.flip(nonaccum_convmoist_nohgtqs_noshal_trans, axis = 0))
# plt.ylim([0,3000])
# plt.colorbar()
# =============================================================================
levels1 = np.linspace(-150, 150, 11)
plt.figure(figsize = (20,10))
plt.contourf(conv_dry_nohgtqs_mean.time.values, heights['Full heights'], np.flip(nonaccum_turb_nohgtqs_trans, axis = 0), vmin = -150, vmax = 150, levels = levels1)
plt.ylim([0,3000])
plt.colorbar(label = '$\Theta_l$ flux [K m-1 s-1]')
plt.ylabel('Height [m]')
plt.xlabel('Time [UTC]')
plt.title('Timeseries of $\Theta_l$ turbulent flux for HARMONIE noHGTQS')

plt.figure(figsize = (20,10))
plt.contourf(conv_dry_nohgtqs_mean.time.values, heights['Full heights'], np.flip(nonaccum_turb_nohgtqs_noshal_trans, axis = 0), vmin = -150, vmax = 150, levels = levels1)
plt.ylim([0,3000])
plt.colorbar(label = '$\Theta_l$ flux [K m-1 s-1]')
plt.ylabel('Height [m]')
plt.xlabel('Time [UTC]')
plt.title('Timeseries of $\Theta_l$ turbulent flux for HARMONIE noHGTQS noSHAL')



#%%
'''Tranpose the data HERE STARTS ACCUMULATED DATA'''
conv_dry_nohgtqs_mean_trans = conv_dry_nohgtqs_mean.transpose()
conv_dry_nohgtqs_noshal_trans = conv_dry_nohgtqs_noshal_mean.transpose()
conv_moist_nohgtqs_trans = conv_moist_nohgtqs_mean.transpose()
conv_moist_nohgtqs_noshal_trans = conv_moist_nohgtqs_noshal_mean.transpose()
turb_nohgtqs_trans = turb_nohgtqs_mean.transpose()
turb_nohgtqs_noshal_trans = turb_nohgtqs_noshal_mean.transpose()


#%%
'''Plot the data'''
plt.figure(figsize = (15,10))
plt.contourf(conv_dry_nohgtqs_mean_trans.time.values, heights['Full heights'], np.flip(conv_dry_nohgtqs_mean_trans.Thlflx_conv_dry.values, axis = 0))
plt.colorbar(label = 'theta_l Flux [K/ms]')
plt.title('Accumulated theta_l Flux Dry Convection for HARMONIE noHGTQS')
plt.ylabel('Height [m]')
plt.xlabel('Time')
plt.ylim([0,3000])
plt.show()

plt.figure(figsize = (15,10))
plt.contourf(conv_dry_nohgtqs_noshal_trans.time.values, heights['Full heights'], np.flip(conv_dry_nohgtqs_noshal_trans.Thlflx_conv_dry.values, axis = 0))
plt.colorbar(label = 'theta_l Flux [K/ms]')
plt.title('Accumulated theta_l Flux Dry Convection for HARMONIE noHGTQS noSHAL')
plt.ylabel('Height [m]')
plt.xlabel('Time')
plt.ylim([0,3000])
plt.show()


plt.figure(figsize = (15,10))
plt.contourf(conv_moist_nohgtqs_trans.time.values, heights['Full heights'], np.flip(conv_moist_nohgtqs_trans.Thlflx_conv_mois.values, axis = 0))
plt.colorbar(label = 'theta_l Flux [K/ms]')
plt.title('Accumulated theta_l Flux Moist Convection for HARMONIE noHGTQS')
plt.ylabel('Height [m]')
plt.xlabel('Time')
plt.ylim([0,3000])
plt.show()


plt.figure(figsize = (15,10))
plt.contourf(conv_moist_nohgtqs_noshal_trans.time.values, heights['Full heights'], np.flip(conv_moist_nohgtqs_noshal_trans.Thlflx_conv_mois.values, axis = 0))
plt.colorbar(label = 'theta_l Flux [K/ms]')
plt.title('Accumulated theta_l Flux moist Convection for HARMONIE noHGTQS noSHAL')
plt.ylabel('Height [m]')
plt.xlabel('Time')
plt.ylim([0,3000])
plt.show()


plt.figure(figsize = (15,10))
plt.contourf(turb_nohgtqs_trans.time.values, heights['Full heights'], np.flip(turb_nohgtqs_trans.Thlflx_turb.values, axis = 0))
plt.colorbar(label = 'theta_l Flux [K/ms]')
plt.title('Accumulated theta_l Flux turbulent for HARMONIE noHGTQS')
plt.ylabel('Height [m]')
plt.xlabel('Time')
plt.ylim([0,3000])
plt.show()


plt.figure(figsize = (15,10))
plt.contourf(turb_nohgtqs_noshal_trans.time.values, heights['Full heights'], np.flip(turb_nohgtqs_noshal_trans.Thlflx_turb.values, axis = 0), vmin = -200, vmax = 150)
plt.colorbar(label = 'theta_l Flux [K/ms]')
plt.title('Accumulated theta_l Flux turbulent HARMONIE noHGTQS noSHAL')
plt.ylabel('Height [m]')
plt.xlabel('Time')
plt.ylim([0,3000])
plt.show()
        
