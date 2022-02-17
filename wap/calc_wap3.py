#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate weighted avg precipitation
and flood days
Based on Lu 2009 and CHen 2021 paper


Defined as:
    WAP  = (1 - a)SUM1 to N, a^n * Pn
    where a is a weighting factor - start at a = 0.9
    N = 44 (weighted avg over 44 days)

1. Import daily rainfall data
2. Calculate WAP
3. Calculate flood days (WAP = 95% percentile)


Created on Jan 31 2022

@author: earsch
"""


# set wd and import packages
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import iris.coord_categorisation
from iris.util import rolling_window
from iris.analysis import Aggregator

from cf_units import Unit

import numpy as np
import numpy.ma as ma

import cartopy.crs as ccrs

import copy

import glob

import pandas as pd

import matplotlib.pyplot as plt

#Import my functions
import sys
sys.path.append('/nfs/see-fs-02_users/earsch/Documents/Leeds/Repos/Tanga/Plot_functions')
import tanzania1 as tp
sys.path.append('/nfs/see-fs-02_users/earsch/Documents/Leeds/Repos/Tanga/Onset_functions')
from onset_functions import masking
from onset_functions import find_season

proj = ccrs.PlateCarree(central_longitude = 38)



#%% Load data
   
# set up save details based on laoded data
save_path = '/nfs/a321/earsch/floods_heatwaves/processed/wap/pan_africa/wap/'     
mod = 'p25'
scen = 'rcp85'
area = 'pan_africa'
#area = 'wa'

print('Loading...')
print(mod , scen, sep = '')

#P25 

file_name = '/nfs/a321/earsch/floods_heatwaves/input_data/pr/pr_' + mod + '_daily_' + scen + '.nc'

pr = iris.load_cube(file_name)

#pr = iris.load_cube('/nfs/a321/earsch/floods_heatwaves/input_data/pr/pr_p25_daily_histo.nc')
#pr = iris.load_cube('/nfs/a321/earsch/floods_heatwaves/input_data/pr/pr_p25_daily_rcp85.nc')

#CP4
#pr = iris.load_cube('/nfs/a321/earsch/floods_heatwaves/input_data/pr/pr_cp4_p25grid_histo.nc')
#pr = iris.load_cube('/nfs/a321/earsch/floods_heatwaves/input_data/pr/pr_cp4_p25grid_rcp85.nc')

#convert to mm/day
pr.convert_units('kg m-2 day-1')

#import landsea mask
cs = pr.coord_system(iris.coord_systems.CoordSystem)


ls = iris.load_cube('/nfs/a277/IMPALA/data/4km/ANCILS/landseamask_ancil_4km_regrid.nc')
ls.coord('longitude').points = ls.coord('longitude').points - 360
ls.coord('longitude').guess_bounds()
ls.coord('latitude').guess_bounds()
try:
    pr.coord('longitude').guess_bounds()
    pr.coord('latitude').guess_bounds()
except:
    print('has bounds')
ls.coord(axis='x').coord_system = cs
ls.coord(axis='y').coord_system = cs
ls_regrid = ls.regrid(pr, iris.analysis.AreaWeighted())
ls_regrid =ls_regrid[0,0]


#%% Extract area

if area == 'wa':
    print('Extracting to West Africa')
    #chagne save path if extra to region
    save_path = '/nfs/a321/earsch/floods_heatwaves/processed/wap/west_africa/wap/'    
    
    
    min_lat = 3.5
    max_lat = 20.0
    min_lon = -20.0
    max_lon = 16.0 
    
    
    cons = iris.Constraint(latitude = lambda cell: min_lat < cell < max_lat,
                           longitude = lambda cell: min_lon < cell < max_lon)
    
    pr = pr.extract(cons)
    
    ls_regrid = ls_regrid.extract(cons)
else:
    print('Preparing for pan-africa')


#%% calc WAP



def wap(pr_in, window, alpha, weight_array):
    ''' Calculate wap based on a time series of precip data, in numpy array form
    
    wap(n) = (1 - a) * [SUM from day n to N] a^n * Pn
    N = window length, refers to previous days
    Days further in past wieghted lower
    '''

    #create aray to save data
    wap = np.zeros(pr_in.shape)
    
    #number of time points in input data
    t = pr_in.shape[0]
   
    #define index to use based on window
    window_index = window -1 #because python starts at 0
    
    #calculate wap - 
    # calculation needs to start at index correspond to window
    # i.e., if window is 44, can't calcualte wap at  33rd day in dataset (not enogun prior data)    
    
    for i in np.arange(window_index, t):
        #for each day, select prior N days, where N = window
        bin_dat = pr_in[i-window_index: i + 1] # data of same length as window
        weighted_bin = bin_dat * weight_array # multiply data by weights
        weighted_sum = np.sum(weighted_bin) # sum to get weighted sum
        
        wap_i = weighted_sum * (1 - alpha)
        
        wap[i] = wap_i # assign wap of day i
        
    return wap
    
    
def apply_wap(pr_cube, window, alpha, ls):
    ''' Apply wap function to an input pr cube, only for alnd areas
    '''
    
    dims = pr_cube.shape
    
    #create cube for saving data
    wap_cube = copy.deepcopy(pr_cube)
    wap_cube.data = np.zeros((dims[0], dims[1], dims[2]))

    #create array of weights
    # a^n, for n in 1 to window (previous days weighted lower)
    n_array = np.arange(1, window + 1)
    weight_array = alpha ** n_array
    #reverse weights array, so that n = 1 referes to first day before present day, rather than
    # meaning first day within window
    weight_array = weight_array[::-1]
    

    #cycle through each lot lon land point, caclulate wap and save
    for lat in np.arange(dims[1]):
        for lon in np.arange(dims[2]):
            
            if ls[lat,lon].data > 0.5:
            
                pr_in = pr_cube[:, lat, lon].data
                #calculate wap as lat, lon location
                wap_loc = wap(pr_in, window, alpha, weight_array)
                
                wap_cube.data[:, lat, lon] = wap_loc
                            

    return wap_cube




#%% calc wap
    
window = 44
alpha = 0.9
    
print('calculating wap with:')
print('Window', window, sep = ' ')
print('Alpha', alpha, sep = '')


if area == 'wa':
    wap_cube = apply_wap(pr, window, alpha, ls_regrid)

    #save data
    
    print('Saving data')
    save_name = save_path + mod + '_' + scen + '_w' + str(window) + '_a' + str(alpha) + '.nc'
    
    iris.save(wap_cube, save_name)
else:
    #for pan africa split into parts (need to do all timesteps at once for wap due to window)
    dims = pr.shape
    
    start_idx = np.arange(0, dims[1], 50)
    end_idx = start_idx + 50
    end_idx[-1] = dims[1]
    
    n_parts = len(start_idx)
    for k in np.arange(0, n_parts):
        print('Starting wap part ', k)
        new_pr = pr[:, start_idx[k]:end_idx[k], :]
        wap_cube = apply_wap(new_pr, window, alpha, ls_regrid)
        
        print('Saving wap part ', k)
        save_name = save_path + 'partial_files/' + mod + '_' + scen + '_w' + str(window) + '_a' + str(alpha) + '_part' + str(k) + '.nc'
        
        iris.save(wap_cube, save_name)
    

    


