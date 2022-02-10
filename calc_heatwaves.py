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
   
#P25 
temp = iris.load_cube('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_p25_historical.nc')
#temp = iris.load_cube('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_p25_rcp85.nc')

#CP4
#temp = iris.load_cube('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_cp4_historical_p25grid.nc')
#temp = iris.load_cube('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_cp4_rcp85_p25grid.nc')

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

save_path = '/nfs/a321/earsch/floods_heatwaves/processed/heatwaves/pan_africa/td/thres/'

mod = 'p25'
scen = 'historical'

per_25 = calc_thres(temp, 25, save_path, mod, scen)
per_75 = calc_thres(temp, 75, save_path, mod, scen)
per_95 = calc_thres(temp, 95, save_path, mod, scen)
per_975 = calc_thres(temp, 97.5, save_path, mod, scen)

    

#%%


def get_heatwave(cube, hw_thres, nconsecutive_days):
    #hw_thres = cube with temp that corresponds to threshold define heatwave
    
    shape=np.shape(cube)
    nt=shape[0]
    nlat=shape[1]
    nlon=shape[2]

 
    count_heatwave_days = copy.deepcopy(cube[0])
    count_heatwave_days.data = np.zeros((nlat,nlon))
    heatwave_duration = copy.deepcopy(cube[0])
    heatwave_duration.data = np.zeros((nlat,nlon))
    heatwaves = copy.deepcopy(cube)
    heatwaves.data = np.zeros((nt,nlat,nlon))
    
    for y in range(nlat):
        for x in range(nlon):
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
            
            
            #get avg heatwave duration
            heatwave_dur = np.mean(len_heat[len_heat >= nconsecutive_days])
            heatwave_duration.data[y,x] = heatwave_dur
            
            #get total numb heatwave days
            n_heatwave_days = np.sum(len_heat[len_heat >= nconsecutive_days])
            count_heatwave_days.data[y,x] = n_heatwave_days
            
            
            
            #find heat events longer or equal to ncsonecutive days, and crate
            # cube with 1s for each day in heatave
            #print(len_heat)
            for i in range(nstart):
                if len_heat[i] >= nconsecutive_days:
                    #heatwave[start_hot[i]:end_hot[i],y,x]=1
                    #length_heatwave[start_hot[i], y,x]=len_heat[i]
                    #print(len_heat[i])
                    #print(count_heatwave_days[y,x], nstart)
                    
                    #get indices assocaited with heatwave
                    heatwave_days = np.arange(start_hot[i], end_hot[i])
                    
                    for index in heatwave_days:
                        heatwaves.data[index, y, x] = 1.0
                    
                   

    return heatwaves, heatwave_duration, count_heatwave_days




#%% calc heatwaves
    
hw, hw_d, count_hw = get_heatwave(temp, per_95, 3.0)

#save data

save_path = '/nfs/a321/earsch/floods_heatwaves/processed/heatwaves/pan_africa/td/attributes/'    

mod = 'p25'
scen = 'historical'

#duration
var_name = 'duration'
save_name = save_path + var_name + '_' + mod + '_' + scen + '.nc'

iris.save(hw_d, save_name)

#coutndays
var_name = 'ndays'
save_name = save_path + var_name + '_' + mod + '_' + scen + '.nc'

iris.save(count_hw, save_name)

#heatwaves
var_name = 'heatwaves'
save_name = save_path + var_name + '_' + mod + '_' + scen + '.nc'

iris.save(hw, save_name)


#%% Calc hetawve magnitude

def heatwave_mag(cube, heatwave_cub, per_25, per_75):
    
    #calc mag
    
    #multiple by days actually heatwave
    
    #sum up for heatwave days