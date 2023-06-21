# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 14:58:42 2023

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
from cut_subgrid_module import cut_sub_grid
import cv2
import glob
import re

img_array = []
numbers = re.compile(r'(\d+)')

                     
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

for filename in sorted(glob.glob(r'C:\Users\LENOVO\Desktop\TU_Delft\thesis\result_plots\mass_flux_nouvmix\*.png') , key=numericalSort):
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width,height)
    img_array.append(img)

out = cv2.VideoWriter('massfulx_noUVmix_video2.mp4', cv2.VideoWriter_fourcc(*'DIVX'), 15, size)

for i in range(len(img_array)):
    out.write(img_array[i])
out.release()

