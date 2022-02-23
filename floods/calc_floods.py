#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Get floods
-number of flood days per year
-duration of heatwaves
-location of floods

Defined as:
WAP >= 95th WAP


Created on Wed Mar 17 09:28:31 2021

@author: earsch
"""


# set wd and import packages
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import iris.coord_categorisation

from cf_units import Unit

import numpy as np
import numpy.ma as ma

import cartopy.crs as ccrs

import copy

import glob


#Import my functions
import sys
sys.path.append('/nfs/see-fs-02_users/earsch/Documents/Leeds/Repos/floods_heatwaves/floods_heatwaves_shared/floods')
import flood_functions as f


proj = ccrs.PlateCarree(central_longitude = 38)



#%% Load data

# set up save details based on laoded data
save_path = '/nfs/a321/earsch/floods_heatwaves/processed/wap/pan_africa/floods/restrict/'  
save_path_thres = '/nfs/a321/earsch/floods_heatwaves/processed/wap/pan_africa/thres/'  
mod = 'p25'
scen = 'rcp85'
area = 'west_africa'

fname = '/nfs/a321/earsch/floods_heatwaves/processed/wap/' + area + '/wap/' + mod + '_' + scen + '*.nc'
wap = iris.load_cube(fname)




#add aux coords
wap.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
wap.add_aux_coord(iris.coords.AuxCoord(scen, long_name = 'sim'))



#import landsea mask
cs = wap.coord_system(iris.coord_systems.CoordSystem)
ls = iris.load_cube('/nfs/a277/IMPALA/data/4km/ANCILS/landseamask_ancil_4km_regrid.nc')
ls.coord('longitude').points = ls.coord('longitude').points - 360
ls.coord('longitude').guess_bounds()
ls.coord('latitude').guess_bounds()
ls.coord(axis='x').coord_system = cs
ls.coord(axis='y').coord_system = cs
ls_regrid = ls.regrid(wap, iris.analysis.AreaWeighted())
ls_regrid =ls_regrid[0,0]

#%% Extract area
# change save path if run

if area == 'west_africa':

    save_path = '/nfs/a321/earsch/floods_heatwaves/processed/wap/west_africa/floods/' 
    save_path_thres = '/nfs/a321/earsch/floods_heatwaves/processed/wap/west_africa/thres/'
    
    min_lat = 3.5
    max_lat = 20.0
    min_lon = -20.0
    max_lon = 16.0 
    
    
    cons = iris.Constraint(latitude = lambda cell: min_lat < cell < max_lat,
                           longitude = lambda cell: min_lon < cell < max_lon)
    
    wap = wap.extract(cons)
    
    ls_regrid = ls_regrid.extract(cons)
#%% drop first half of 1997 and second hafl of 2006
# 1st June 1997 = 150
# 1st June 2006 = 3390

#will give 9 years exactly, starting 1st juen 1997, ending 1st june (exclusive)
wap = wap[150:3390,:,:]



#%% #calculate threshold

#thres    
per_95 = f.calc_thres(wap, 95, save_path_thres, mod, scen)

#%% calcualte floods

restrict = 10.0

#floods
output = f.get_floods(wap, per_95, ls_regrid)


#save data
var_names = ['floods', 'duration', 'ndays', 'nevents', 'start', 'end']


for i in np.arange(len(output)):

    save_name = save_path +  var_names[i] + '_' + mod + '_' + scen + '_r' + str(restrict) +'.nc'

    iris.save(output[i], save_name)
