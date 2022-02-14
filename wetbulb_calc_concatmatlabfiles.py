#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''Data Preparation  - Put togethet matlab wetbulb output files
Date: 11/02/2022
Author: Sarah C

            
'''

#%% Import packages

import matplotlib.pyplot as plt 

import iris
import iris.plot as iplt
import iris.quickplot as qplt
from iris.experimental.equalise_cubes import equalise_attributes
import iris.coord_categorisation
from cf_units import Unit

import cartopy.crs as ccrs

import glob

import copy

import numpy as np

import metpy
import metpy.calc as mpcalc
from metpy.units import units

import itertools

import sys
sys.path.append('/nfs/see-fs-02_users/earsch/Documents/Leeds/Repos/floods_heatwaves/wet_bulb_func')
import WetBulb as wb

proj = ccrs.PlateCarree(central_longitude = 38)

#%% Import data
#update model name and scenario name in filenames, and saved path at end of script

#temp_orig = iris.load_cube('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_p25_rcp85.nc')
temp_orig = iris.load_cube('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_cp4_historical_p25grid.nc')


file_path = ('/nfs/a321/earsch/floods_heatwaves/processed/wetbulb_temp/')
filenames = glob.glob(file_path + '*cp4_hist*part?.nc')
files = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    files.append(x)
    
#%%

concat_files = files.concatenate_cube()
concat_files.units = 'celsius'
concat_files.coord('latitude').units = Unit('degrees')
concat_files.coord('longitude').units = Unit('degrees')


#%% regrid and extarct to west africa

min_lat = 3.5
max_lat = 20.0
min_lon = -20.0
max_lon = 16.0 


cons = iris.Constraint(latitude = lambda cell: min_lat < cell < max_lat,
                       longitude = lambda cell: min_lon < cell < max_lon)

temp_orig = temp_orig.extract(cons)

try:
    temp_orig.coord('longitude').guess_bounds()
    temp_orig.coord('latitude').guess_bounds()
except:
    print('has bounds')
    
concat_files.coord('longitude').guess_bounds()
concat_files.coord('latitude').guess_bounds()

cs = temp_orig.coord_system(iris.coord_systems.CoordSystem)
concat_files.coord('longitude').coord_system = cs
concat_files.coord('latitude').coord_system = cs

concat_regrid = concat_files.regrid(temp_orig[0], iris.analysis.AreaWeighted())


#%%
#put data fromc conat files into P25 temp, so get coord info

temp_orig.data[0:3600, :,:] = concat_regrid.data


#%%

iris.save(temp_orig, '/nfs/a321/earsch/floods_heatwaves/processed/wetbulb_temp/wb_cp4_daily_historical_wa.nc')

