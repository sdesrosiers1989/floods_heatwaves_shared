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
area = 'pan_africa'

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


#%% calc threshold


def calc_thres(cube, per, save_path, mod, scen, save = True):
    thres = cube.collapsed('time', iris.analysis.PERCENTILE, percent = per)
    
    if save == True:
        save_name = save_path + mod + '_' + scen + '_' + str(per) + '.nc'
        iris.save(thres, save_name)
        
    return thres


def get_floods(cube, wap_thres, ls, restrict = 0.0):
    ''' Calculate floods events
        Defined as wap > 95th percentile
        Calcualted only for land points
        
        Find duration, number of flood events, flood days,
        and start and end point of floods
    
        cube = input weighted area percip cube
        wap_thres = wap threshold which defines floods (i.e., 95th percnentile)
        ls = landsea mask

    '''
    
    ## Create cubes to save output data
    
    shape=np.shape(cube)
    nt=shape[0]
    nlat=shape[1]
    nlon=shape[2]

    # total number of flood days and individual flodo events for each gridpoint within time series
    count_flood_days = copy.deepcopy(cube[0])
    count_flood_days.data = np.zeros((nlat,nlon))
    count_flood_event = copy.deepcopy(cube[0])
    count_flood_event.data = np.zeros((nlat,nlon))
    # will create cube with for duration of flood
    # at first day of each flood (0s everything else)
    flood_duration = copy.deepcopy(cube)
    flood_duration.data = np.zeros((nt,nlat,nlon))
    #start and and end day of floods
    start_flood_event = copy.deepcopy(cube)
    start_flood_event.data = np.zeros((nt,nlat,nlon))
    end_flood_event = copy.deepcopy(cube)
    end_flood_event.data = np.zeros((nt,nlat,nlon))
 
    
    #find flood events
    floods = copy.deepcopy(cube)
    floods.data = np.where(floods.data > wap_thres.data, 1.0, 0.0)     
    
    wap_greater_restrict = copy.deepcopy(cube)
    wap_greater_restrict.data = np.where(wap_greater_restrict.data >= restrict, 1.0, 0.0) 

    for y in range(nlat):
        for x in range(nlon):
            if ls[y,x].data > 0.5:
        
                
                f_loc = floods[:, y, x].data
                r_loc = wap_greater_restrict[:,y,x].data
                 
                f=np.zeros(nt+2) # add an element at beginning and end compared to floods
                
                #find locations of floods (both threshold and mm restirction are true)
                fidx = np.where(f_loc + r_loc == 2.0)
           
                f[fidx[0] + 1] = 1 # add 1 to every index after flood defined

                diffs=np.diff(f) 
                # this is 1 where it becomes flood and -1 where it becomes not flood 
                # and 0 where it stays the same
                
                # get indices where changes occur and find length of flood events
                start_flood=np.where(diffs>0)[0]
                end_flood=np.where(diffs<0)[0] # end flood day - day flood ends (not a flood day)
                nstart=len(start_flood)
                #nend=len(end_hot)
                len_flood=end_flood-start_flood
                
               
                #get total num flood days
                n_flood_days = np.sum(len_flood)
                count_flood_days.data[y,x] = n_flood_days
                
                #number of individual flood events
                count_flood_event.data[y,x] = nstart



                # flood_duration:  cube with length of flood at first day of flood
                #  start flood
                #  end flood
                for i in range(nstart):
    
                    #flood duration
                    flood_duration.data[start_flood[i], y,x]=len_flood[i]
                    
                    
                    #assign 1s to start and end dates of heatwaves
                    start_flood_event.data[start_flood[i], y,x] = 1
                    end_flood_event.data[end_flood[i] - 1, y,x] = 1 # last day of flood, 
                                                                    # is a flood day
                    
          
                        
                        
                        

    return [floods, flood_duration, count_flood_days, count_flood_event, start_flood_event, end_flood_event]



                

#%%

#calculate flods

#thres    
per_95 = calc_thres(wap, 95, save_path_thres, mod, scen)

#%%

restrict = 10.0

#floods
output = get_floods(wap, per_95, ls_regrid)


#save data
var_names = ['floods', 'duration', 'ndays', 'nevents', 'start', 'end']


for i in np.arange(len(output)):

    save_name = save_path +  var_names[i] + '_' + mod + '_' + scen + '_r' + str(restrict) +'.nc'

    iris.save(output[i], save_name)
