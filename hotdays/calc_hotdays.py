#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Get number of hotdays
    

Defined as:
    Days with Tav >= 95th percentile Tmean, defined over all days in dataset
        -check also 97.5th percentile


Created on Wed Mar 17 09:28:31 2021

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



#%% User inputs - model, temp_type and scenario

# set up save details based on laoded data
temp_type = 'td'
#temp_type = 'tw' 
mod = 'p25'
scen = 'rcp85'
scen = 'historical'
area = 'pan_africa'
hd_thres = '95' # threshold which defines hotdays

#%% Load data basedon user inputs

save_path = '/nfs/a321/earsch/floods_heatwaves/processed/hotdays/pan_africa/' + temp_type + '/' 
file_path_thres = '/nfs/a321/earsch/floods_heatwaves/processed/heatwaves/pan_africa/' + temp_type + '/thres/'  
#95th percentile already calcualted as part of hetawve script, don't need to redo

## td temps
#P25 
if temp_type == 'td':
    if mod == 'p25':
        fname = '/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_p25_' + scen + '.nc'
    elif mod == 'cp4':
        fname = '/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_cp4_' + scen + '_p25grid.nc'
elif temp_type == 'tw':
    fname = '/nfs/a321/earsch/floods_heatwaves/processed/wetbulb_temp/' + area + '/wb_' + mod + '_daily_'  + scen + '*.nc'
       
        
temp = iris.load_cube(fname)


#add aux coords
temp.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
temp.add_aux_coord(iris.coords.AuxCoord(scen, long_name = 'sim'))

#convert Temp to celsius
temp.convert_units('celsius')


#load threshold 
thres = iris.load_cube(file_path_thres + mod + '_historical_' + hd_thres + '.nc')



#%% drop first half of 1997 and second hafl of 2006
# 1st June 1997 = 150
# 1st June 2006 = 3390

#will give 9 years exactly, starting 1st juen 1997, ending 1st june (exclusive)
#this works for td and Tw - though tw only goes up to 3600 rather than 3660
temp = temp[150:3390,:,:]

#%%

def get_heatwave(cube, thres):
    ''' Calculate locations of hotdays, and total number of hotdays
        
        cube = input temperature cube
        thres = Temperature threshold which defines hotdays (i.e., 95th percnentile)
    '''
    

    
    # will create cube with 1 for each hotday, 0 everything else
    hotdays = copy.deepcopy(cube)
    hotdays.data = np.where(hotdays.data > thres.data, 1.0, 0.0)
    
    # total number of ho days for each gridpoint within time series
    count_hot_days = hotdays.collapsed('time', iris.analysis.SUM)
    

    return [hotdays, count_hot_days]



#%%
#calculate hotdays


output = get_heatwave(temp, thres)

#save data
var_names = ['hotdays', 'ndays']


for i in np.arange(len(output)):

    save_name = save_path + 'per' + hd_thres + '/' + var_names[i] + '_' + mod + '_' + scen + '.nc'

    iris.save(output[i], save_name)


