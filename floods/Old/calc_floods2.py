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
save_path = '/nfs/a321/earsch/floods_heatwaves/processed/wap/pan_africa/floods/'  
save_path_thres = '/nfs/a321/earsch/floods_heatwaves/processed/wap/pan_africa/thres/'  
mod = 'p25'
scen = 'historical'

#P25 
wap = iris.load_cube('/nfs/a321/earsch/floods_heatwaves/processed/wap/west_africa/wap/p25*hist*')
#wap = iris.load_cube('/nfs/a321/earsch/floods_heatwaves/processed/wap/west_africa/wap/p25*rcp85*')

#CP4
#wap = iris.load_cube('/nfs/a321/earsch/floods_heatwaves/processed/wap/west_africa/wap/cp4*hist*')
#wap = iris.load_cube('/nfs/a321/earsch/floods_heatwaves/processed/wap/west_africa/wap/cp4*rcp85*')

#add aux coords
wap.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
wap.add_aux_coord(iris.coords.AuxCoord(scen, long_name = 'sim'))



#%% Extract area
# change save path if run

save_path = '/nfs/a321/earsch/floods_heatwaves/processed/wap/west_africa/floods/' 
save_path_thres = '/nfs/a321/earsch/floods_heatwaves/processed/wap/west_africa/thres/'

min_lat = 3.5
max_lat = 20.0
min_lon = -20.0
max_lon = 16.0 


cons = iris.Constraint(latitude = lambda cell: min_lat < cell < max_lat,
                       longitude = lambda cell: min_lon < cell < max_lon)

wap = wap.extract(cons)



#%% calc threshold


def calc_thres(cube, per, save_path, mod, scen, window = 44., save = True):
    #ignore first days within window as these are forced to 0 - not enoguh data calculate wap
    thres = cube[window:].collapsed('time', iris.analysis.PERCENTILE, percent = per)
    
    if save == True:
        save_name = save_path + mod + '_' + scen + '_' + str(per) + '.nc'
        iris.save(thres, save_name)
        
    return thres


def get_floods(cube, wap_thres):
    ''' Calculate floods events
        Defined as wap > 95th percentile
        Calcualted only for land points
    
        cube = input weighted area percip cube
        wap_thres = wap threshold which defines floods (i.e., 95th percnentile)

    '''
    
    floods = copy.deepcopy(cube)
    floods.data = np.where(floods.data > wap_thres.data, 1.0, 0.0)           

    return floods




#%%
#calculate flods

#thres    
per_95 = calc_thres(wap, 95, save_path_thres, mod, scen)

#floods
floods = get_floods(wap, per_95)
save_name = save_path + 'floods_' + mod + '_' + scen + '.nc'
iris.save(floods, save_name)


