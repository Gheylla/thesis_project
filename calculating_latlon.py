# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 20:34:13 2023

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
#from pyhdf.SD import SD, SDC
#import GOES

def calculate_degrees(data):
    proj_info = data.variables['goes_imager_projection'].attrs
    lon_origin = proj_info['longitude_of_projection_origin']
    H = proj_info['perspective_point_height'] + proj_info['semi_major_axis']
    r_eq = proj_info['semi_major_axis']
    r_pol = proj_info['semi_minor_axis']

    lat_rad_1d = data.variables['x'][:]
    lon_rad_1d = data.variables['y'][:]

    lat_rad,lon_rad = np.meshgrid(lat_rad_1d,lon_rad_1d)

    lambda_0 = (lon_origin*np.pi)/180.0

    a_var = np.power(np.sin(lat_rad),2.0) + (np.power(np.cos(lat_rad),2.0)*(np.power(np.cos(lon_rad),2.0)+(((r_eq*r_eq)/(r_pol*r_pol))*np.power(np.sin(lon_rad),2.0))))
    b_var = -2.0*H*np.cos(lat_rad)*np.cos(lon_rad)
    c_var = (H**2.0)-(r_eq**2.0)

    r_s = (-1.0*b_var - np.sqrt((b_var**2)-(4.0*a_var*c_var)))/(2.0*a_var)

    s_x = r_s*np.cos(lat_rad)*np.cos(lon_rad)
    s_y = - r_s*np.sin(lat_rad)
    s_z = r_s*np.cos(lat_rad)*np.sin(lon_rad)

    lat = (180.0/np.pi)*(np.arctan(((r_eq*r_eq)/(r_pol*r_pol))*((s_z/np.sqrt(((H-s_x)*(H-s_x))+(s_y*s_y))))))
    lon = (lambda_0 - np.arctan(s_y/(H-s_x)))*(180.0/np.pi)

    # Ignore numpy errors for sqrt of negative number; occurs for GOES-16 ABI CONUS sector data
    np.seterr(all='ignore')

    
    return lat, lon