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

import os

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

model = 'cp4'
scen = 'histo'
m_scen = model + '_' + scen

if scen == 'histo':
    t_scen = 'historical'
else:
    t_scen = scen

if model == 'p25':
    temp_orig = iris.load_cube('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_p25_' + t_scen + '.nc')
elif  model == 'cp4':
    temp_orig = iris.load_cube('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_cp4_' + t_scen + '_p25grid.nc')


file_path = ('/nfs/a321/earsch/floods_heatwaves/processed/wetbulb_temp/')
filenames = glob.glob(file_path + '*{}*part?.nc'.format(m_scen))
files = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    files.append(x)    

print('Loading:')
for f in filenames:
    print(f)

    
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
temp_new = temp_orig[0:3600]


#%%

save_name = '/nfs/a321/earsch/floods_heatwaves/processed/wetbulb_temp/wb_' + model + '_daily_' + t_scen + '_wa.nc'

print('Saving as ', save_name)

iris.save(temp_new, save_name)

