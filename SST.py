# -*- coding: utf-8 -*-
"""
Created on Wed May  3 10:08:57 2023

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
data = cat.RonBrown.MAERI['SSTskin'].to_dask()

SST = data.sea_surface_temperature.values
mean_SST = np.mean(SST)
