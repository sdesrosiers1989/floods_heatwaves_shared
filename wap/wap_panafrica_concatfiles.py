#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''Data Preparation  - Put togethet pan-africa wap files
Date: 11/02/2022
Author: Sarah C

            
'''

#%% Import packages

import matplotlib.pyplot as plt 

import iris
import iris.plot as iplt
import iris.quickplot as qplt
from iris.util import equalise_attributes
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

def add_att_from_filename(cube, field, filename):
    file = filename[62:]
    mod = filename[62:65]
    #split filename into sections by '_', find second partition and split again
    sim = file.partition('_')[2].partition('_')[2].partition('_')[0]
    
    cube.add_aux_coord(iris.coords.AuxCoord('um', long_name = 'gcm'))
    cube.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
    cube.add_aux_coord(iris.coords.AuxCoord(sim, long_name = 'sim'))

proj = ccrs.PlateCarree(central_longitude = 38)

#%% Import data
#update model name and scenario name in filenames, and saved path at end of script

model = 'cp4'
scen = 'rcp85'
scen = 'histo'
wind_alpha = 'w44_a0.9'
m_scen = model + '_' + scen + '_' + wind_alpha

file_path = ('/nfs/a321/earsch/floods_heatwaves/processed/wap/pan_africa/wap/partial_files/')
filenames = glob.glob(file_path + '*{}*part?.nc'.format(m_scen))
files = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    files.append(x)    

print('Loading:')
for f in filenames:
    print(f)


    
#%%

equalise_attributes(files)
concat_files = files.concatenate_cube()


#%%

save_name = '/nfs/a321/earsch/floods_heatwaves/processed/wap/pan_africa/wap/' + m_scen + '.nc'

print('Saving as ', save_name)

iris.save(concat_files, save_name)

