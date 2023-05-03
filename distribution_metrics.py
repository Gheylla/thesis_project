# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 10:24:18 2023

@author: Gheylla Liberia
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


'''Import data'''
data_noHGTQS = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\df_metrics_noHGTQS.h5')
data_noSHAL_noHGTQS = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\df_metrics_noHGTQS_noSHAL.h5')
data_noUVmix_noHGTQS = pd.read_hdf(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\metric_results\df_metrics_noHGTQS_noUVmix.h5')

metrics = data_noHGTQS.columns

# =============================================================================
# 
# for i in range(len(metrics)): 
#     plt.figure()
#     sn.distplot(data_noHGTQS[metrics[i]], label = 'noHGTQS', hist = False)
#     sn.distplot(data_noSHAL_noHGTQS[metrics[i]], label = 'noSHAL_noHGTQS', hist = False)
#     sn.distplot(data_noUVmix_noHGTQS[metrics[i]], label = 'noUVmix_noHGTQS', hist = False)
#     plt.title(f' Distribution plot of {metrics[i]}')
#     plt.grid()
#     plt.legend()
# =============================================================================
