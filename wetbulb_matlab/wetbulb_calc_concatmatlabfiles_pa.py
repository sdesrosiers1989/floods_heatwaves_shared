#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''Data Preparation  - Put togethet matlab wetbulb output files
Date: 11/02/2022
Author: Sarah C

            
'''

#%% Import packages

import iris
import iris.plot as iplt
import iris.quickplot as qplt
from iris.experimental.equalise_cubes import equalise_attributes
import iris.coord_categorisation
from cf_units import Unit

import cartopy.crs as ccrs

import glob



#%% Import data
#update model name and scenario name in filenames, and saved path at end of script

model = 'cp4'
#scen = 'histo'
scen = 'rcp85'
m_scen = model + '_' + scen

if scen == 'histo':
    t_scen = 'historical'
else:
    t_scen = scen

if model == 'p25':
    temp_orig = iris.load_cube('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_p25_' + t_scen + '.nc')
elif  model == 'cp4':
    temp_orig = iris.load_cube('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_cp4_' + t_scen + '_p25grid.nc')


file_path = ('/nfs/a321/earsch/floods_heatwaves/processed/wetbulb_temp/pan_africa/ind_years/')
filenames = glob.glob(file_path + '*{}*part?_pa.nc'.format(m_scen))
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


#%% regrid 

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

save_name = '/nfs/a321/earsch/floods_heatwaves/processed/wetbulb_temp/pan_africa/wb_' + model + '_daily_' + t_scen + '_pa.nc'

print('Saving as ', save_name)

iris.save(temp_new, save_name)

