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

#Import my functions
import sys
sys.path.append('/nfs/see-fs-02_users/earsch/Documents/Leeds/Repos/Tanga/Plot_functions')


proj = ccrs.PlateCarree(central_longitude = 38)



#%% User inputs - model, temp_type and scenario

# set up save details based on laoded data
temp_type = 'td'
#temp_type = 'tw' 
mod = 'cp4'
#scen = 'rcp85'
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

#import landsea mask
cs = temp.coord_system(iris.coord_systems.CoordSystem)
ls = iris.load_cube('/nfs/a277/IMPALA/data/4km/ANCILS/landseamask_ancil_4km_regrid.nc')
ls.coord('longitude').points = ls.coord('longitude').points - 360
ls.coord('longitude').guess_bounds()
ls.coord('latitude').guess_bounds()
try:
	temp.coord('longitude').guess_bounds()
	temp.coord('latitude').guess_bounds()
except:
	print('has bounds')
ls.coord(axis='x').coord_system = cs
ls.coord(axis='y').coord_system = cs
ls_regrid = ls.regrid(temp, iris.analysis.AreaWeighted())
ls_regrid =ls_regrid[0,0]

#%% drop first half of 1997 and second hafl of 2006
# 1st June 1997 = 150
# 1st June 2006 = 3390

#will give 9 years exactly, starting 1st juen 1997, ending 1st june (exclusive)
#this works for td and Tw - though tw only goes up to 3600 rather than 3660
temp = temp[150:3390,:,:]

#%%

def get_heatwave(cube, thres, ls):
    ''' Calculate locations of hotdays, and total number of hotdays
        
        cube = input temperature cube
        thres = Temperature threshold which defines hotdays (i.e., 95th percnentile)
    '''
    

    
    # will create cube with 1 for each hotday, 0 everything else
    hotdays = copy.deepcopy(cube)
    hotdays.data = np.where(hotdays.data > thres.data, 1.0, 0.0)
    
    # total number of ho days for each gridpoint within time series
    count_hot_days = hotdays.collapsed('time', iris.analysis.SUM)
    
    #start and end of multiday events
    shape=np.shape(cube)
    nt=shape[0]
    nlat=shape[1]
    nlon=shape[2]
    
    count_hd_event = copy.deepcopy(cube[0])
    count_hd_event.data = np.zeros((nlat,nlon))
    start_hd_event = copy.deepcopy(cube)
    start_hd_event.data = np.zeros((nt,nlat,nlon))
    end_hd_event = copy.deepcopy(cube)
    end_hd_event.data = np.zeros((nt,nlat,nlon))
    hd_duration = copy.deepcopy(cube)
    hd_duration.data = np.zeros((nt,nlat,nlon))

    for y in range(nlat):
        for x in range(nlon):
            if ls[y,x].data > 0.5:
        
                
                h_loc = hotdays[:, y, x].data
                 
                h=np.zeros(nt+2) # add an element at beginning and end compared to floods
                
                hidx = np.where(h_loc == 1.0)
           
                h[hidx[0] + 1] = 1 # add 1 to every index after flood defined

                diffs=np.diff(h) 
                # this is 1 where it becomes flood and -1 where it becomes not flood 
                # and 0 where it stays the same
                
                # get indices where changes occur and find length of flood events
                start_hd=np.where(diffs>0)[0]
                end_hd=np.where(diffs<0)[0] # end flood day - day flood ends (not a flood day)
                nstart=len(start_hd)
                #nend=len(end_hot)
                len_hd=end_hd-start_hd
                
               #number of individual flood events
                count_hd_event.data[y,x] = nstart

                # hd_duration:  cube with length of hotday events at first hotday

                for i in range(nstart):
    
                    #flood duration
                    hd_duration.data[start_hd[i], y,x]=len_hd[i]
                    
                    
                    #assign 1s to start and end dates of heatwaves
                    start_hd_event.data[start_hd[i], y,x] = 1
                    end_hd_event.data[end_hd[i] - 1, y,x] = 1 # last day of hotday event, 
                                                                    # is a hot day
                    
          
                        
                        
                        

    return [hotdays, hd_duration, count_hot_days, count_hd_event, start_hd_event, end_hd_event]





#%%
#calculate hotdays



if area == 'west_africa':
    output = get_heatwave(temp, thres, ls_regrid)

    #save data
    var_names =  ['hotdays', 'duration', 'ndays', 'nevents', 'start', 'end']
    
    
    for i in np.arange(len(output)):
    
        save_name = save_path + 'per' + hd_thres + '/' + var_names[i] + '_' + mod + '_' + scen + '.nc'
    
        iris.save(output[i], save_name)
else:
    #for pan africa split into parts (need to do all timesteps at once for wap due to window)
    print(mod, scen)
    temp_dims = temp.shape
    
    #print(pr_dims)
    
    start_idx = np.arange(0, temp_dims[1], 50)
    end_idx = start_idx + 50
    end_idx[-1] = temp_dims[1]
    
    n_parts = len(start_idx)
    for k in np.arange(0, n_parts):
        print('Starting hotdays part ', k, start_idx[k], end_idx[k])
        new_t = temp[:, start_idx[k]:end_idx[k], :]
        new_ls = ls_regrid[start_idx[k]:end_idx[k], :]
        new_per = thres[start_idx[k]:end_idx[k], :]
        
        output = get_heatwave(new_t, new_per, new_ls)
                        
        print('Saving hotdays part ', k)
        
        
        #save data
        var_names =  ['hotdays', 'duration', 'ndays', 'nevents', 'start', 'end']
        
        
        for i in np.arange(len(output)):
                            
            save_name = save_path + 'per' + hd_thres + '/partial_files/' + var_names[i] + '_' + mod + '_' + scen + '_part' + str(k) + '.nc'
    
            iris.save(output[i], save_name)
