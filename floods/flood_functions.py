#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Floods functions - functions for calculating flodos

Created on Wed Feb 23 08:15:58 2022

@author: earsch
"""

# set wd and import packages
import iris

import copy

import numpy as np


#define functions

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
    floods_wap_thres = copy.deepcopy(cube)
    floods_wap_thres.data = np.where(floods_wap_thres.data > wap_thres.data, 1.0, 0.0)     
    
    wap_greater_restrict = copy.deepcopy(cube)
    wap_greater_restrict.data = np.where(wap_greater_restrict.data >= restrict, 1.0, 0.0) 
    
    floods = copy.deepcopy(floods_wap_thres)
    floods.data = floods.data + wap_greater_restrict.data
    floods.data = np.where(floods.data == 2.0, 1.0, 0.0)

    for y in range(nlat):
        for x in range(nlon):
            if ls[y,x].data > 0.5:
        
                
                f_loc = floods[:, y, x].data
                 
                f=np.zeros(nt+2) # add an element at beginning and end compared to floods
                
                fidx = np.where(f_loc == 1.0)
           
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



                