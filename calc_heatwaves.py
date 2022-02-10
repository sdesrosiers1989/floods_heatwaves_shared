#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Get heatwaves
-number of heatwave days per year
-duration of heatwaves
-itnensity of heatwaves (heatave magnitude index, Russo 2015)
    Sum of intensity on each day of ehatwve, Md(Td)
    Md(Td) = (Td - T25th) / (T75th- T25th) where Td > T25th
    

Defined as:
    3 or more consectuive days with Tav >= 95th percentile Tmean, defined over all days in dataset
        -check also 97.5th percentile

Reclaculate threshold separeately for future and historic data


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
save_path = '/nfs/a321/earsch/floods_heatwaves/processed/heatwaves/pan_africa/td/attributes/'  
save_path_thres = '/nfs/a321/earsch/floods_heatwaves/processed/heatwaves/pan_africa/td/thres/'  
mod = 'p25'
scen = 'historical'

#P25 
temp = iris.load_cube('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_p25_historical.nc')
#temp = iris.load_cube('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_p25_rcp85.nc')

#CP4
#temp = iris.load_cube('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_cp4_historical_p25grid.nc')
#temp = iris.load_cube('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_cp4_rcp85_p25grid.nc')

#add aux coords
temp.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
temp.add_aux_coord(iris.coords.AuxCoord(scen, long_name = 'sim'))

#convert Temp to celsius
temp.convert_units('celsius')

#import landsea mask
cs = temp.coord_system(iris.coord_systems.CoordSystem)

ls = iris.load_cube('/nfs/a277/IMPALA/data/4km/ANCILS/landseamask_ancil_4km_regrid.nc')
ls.coord('longitude').points = ls.coord('longitude').points - 360
ls.coord('longitude').guess_bounds()
ls.coord('latitude').guess_bounds()
temp.coord('longitude').guess_bounds()
temp.coord('latitude').guess_bounds()
ls.coord(axis='x').coord_system = cs
ls.coord(axis='y').coord_system = cs
ls_regrid = ls.regrid(temp, iris.analysis.AreaWeighted())
ls_regrid =ls_regrid[0,0]


#%% Extract area
# change save path if run

save_path = '/nfs/a321/earsch/floods_heatwaves/processed/heatwaves/west_africa/td/attributes/'  
save_path_thres = '/nfs/a321/earsch/floods_heatwaves/processed/heatwaves/west_africa/td/thres/'

min_lat = 3.5
max_lat = 20.0
min_lon = -20.0
max_lon = 16.0 


cons = iris.Constraint(latitude = lambda cell: min_lat < cell < max_lat,
                       longitude = lambda cell: min_lon < cell < max_lon)

temp = temp.extract(cons)

ls_regrid = ls_regrid.extract(cons)


#%% calc threshold


def calc_thres(cube, per, save_path, mod, scen, save = True):
    thres = cube.collapsed('time', iris.analysis.PERCENTILE, percent = per)
    
    if save == True:
        save_name = save_path + mod + '_' + scen + '_' + str(per) + '.nc'
        iris.save(thres, save_name)
        
    return thres



per_25 = calc_thres(temp, 25, save_path_thres, mod, scen)
per_75 = calc_thres(temp, 75, save_path_thres, mod, scen)
per_95 = calc_thres(temp, 95, save_path_thres, mod, scen)
per_975 = calc_thres(temp, 97.5, save_path_thres, mod, scen)

    

#%%

def magnitude(temp):
    ''' Calculate heatwave mangitude ofr each day in dataset as per Russo 2015
    
        After have magnitude, run function 'get_heatwave'to sum mangitudecs and get mangitude 
        of heatwave days only
        
        Md(Td) = (Td - T25th) / (T75th- T25th) where Td > T25th
        
        temp = cube of temperature
        
    '''
    
    per_25 = temp.collapsed('time', iris.analysis.PERCENTILE, percent = 25)
    per_75 = temp.collapsed('time', iris.analysis.PERCENTILE, percent = 75)
    
    denom = per_75 - per_25

    mag = (temp - per_25) / denom
    
    mag.data = np.where(mag.data < 0, 0, mag.data)
    
    return mag


def get_heatwave(cube, hw_thres, nconsecutive_days, mag, ls):
    ''' Calculate duration, number of heatwave days, locations of heatwave,
        start and and day of heatwaves, and 
         and magnitude of heatwaves.
         Calcualted only for land points
    
        ls = landsea mask
        cube = input temperature cube
        hw_thres = Temperature threshold which defines heataves (i.e., 95th percnentile)
        nconsecutive_days =  number consec days to define a heatve
        mag = magnitude of heatwave for each day in dataset
            to calcualte heatwave magnitude, sum for heatwave days only
    '''
    
    ## Create cubes to save output data
    
    shape=np.shape(cube)
    nt=shape[0]
    nlat=shape[1]
    nlon=shape[2]

 
    # total number of ehatawve days for each gridpoint within time series
    count_heatwave_days = copy.deepcopy(cube[0])
    count_heatwave_days.data = np.zeros((nlat,nlon))
    # will create cube with n days of heatwave and duration of heatwave
    # at first day of each heatwave (0s everything else)
    heatwave_duration = copy.deepcopy(cube)
    heatwave_duration.data = np.zeros((nt,nlat,nlon))
    heatwaves = copy.deepcopy(cube)
    heatwaves.data = np.zeros((nt,nlat,nlon))
    start_heatwave = copy.deepcopy(cube)
    start_heatwave.data = np.zeros((nt,nlat,nlon))
    end_heatwave = copy.deepcopy(cube)
    end_heatwave.data = np.zeros((nt,nlat,nlon))
    heatwave_mag = copy.deepcopy(cube)
    heatwave_mag.data = np.zeros((nt,nlat,nlon))
    
    
    for y in range(nlat):
        for x in range(nlon):
            if ls[y,x].data > 0.5:
                hot=np.zeros(nt+2) # add an element at beginning and end compared to temp
                
                # where is the temp greater than or equal to threshold which defines heatwaves
                hotix=np.where(cube[:,y,x].data >= hw_thres[y,x].data)
                
                hot[hotix[0]+1]=1 # add 1 to every index after hot day defined
                
                diffs=np.diff(hot) 
                # this is 1 where it becomes hot and -1 where it becomes cold and 0 where it stays the same
                
                # get indices where changes occur and find length of hot events
                start_hot=np.where(diffs>0)[0]
                end_hot=np.where(diffs<0)[0] # end hot day - day heatwave ends (not a heatwave)
                nstart=len(start_hot)
                #nend=len(end_hot)
                len_heat=end_hot-start_hot
                
               
                #get total num heatwave days
                n_heatwave_days = np.sum(len_heat[len_heat >= nconsecutive_days])
                count_heatwave_days.data[y,x] = n_heatwave_days
                
                
                
                #find heat events longer or equal to ncsonecutive days, and create:
                # heatwaves (location of heatwaves): cube with 1s for each day in heatave
                # heatwaves_duration:  cube with length of heatwave at first day of heatwave
                # heatwave_mag: sum of daily magnitudes for heatwaves, value at first day of ehatwave
                for i in range(nstart):
                    if len_heat[i] >= nconsecutive_days:
                        
                        #heatwave duration
                        heatwave_duration.data[start_hot[i], y,x]=len_heat[i]
                        
                        #heatwave days - get indices assocaited with heatwave
                        #last day of heatwave_days is a heatwave
                        heatwave_days = np.arange(start_hot[i], end_hot[i])
                        
                        #heatwave magnitude
                        mag_sum = mag.data[heatwave_days, y, x].sum()
                        heatwave_mag.data[start_hot[i], y, x] = mag_sum
                        
                        #assign 1s to start and end dates of heatwaves
                        start_heatwave.data[start_hot[i], y,x] = 1
                        end_heatwave.data[end_hot[i] - 1, y,x] = 1 # last day of heatwave, 
                                                                   # is a heatwave day
                        
                        #assign 1 to each day in heatwave (all others 0)
                        for index in heatwave_days:
                            heatwaves.data[index, y, x] = 1.0
                        
                        
                        

    return [heatwaves, heatwave_duration, count_heatwave_days, heatwave_mag, start_heatwave, end_heatwave]




#%% calc daily magnitude 

#magnitude
mag = magnitude(temp)
var_name = 'dailymagnitude'
save_name = save_path + var_name + '_' + mod + '_' + scen + '.nc'

iris.save(mag, save_name)

#%%
#calculate heatwave

#select which percentile threshold to use
per_cube = per_95
per_name = 'per95/'

per_cube = per_975
per_name = 'per975/'

output = get_heatwave(temp, per_cube, 3.0, mag, ls_regrid)

#save data
var_names = ['heatwaves', 'duration', 'ndays', 'heatwavemagnitude', 'startheatwave', 'endheatwave']


for i in np.arange(len(output)):

    save_name = save_path + per_name + var_names[i] + '_' + mod + '_' + scen + '.nc'

    iris.save(output[i], save_name)


