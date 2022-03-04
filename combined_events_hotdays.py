#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Get combiend floods and ehatwave events
-number of combined events (same day)
-number of floods before start of heatwves (3, 5, 7 days)
-number of floods after end of heatwves (3, 5, 7 days) 


Created on Wed Mar 17 09:28:31 2021

@author: earsch
"""


# set wd and import packages
import iris
import iris.coord_categorisation

import numpy as np

import copy




#%% Load data

# set up save details for loading data
loc = 'pan_africa'
mod = 'p25'
scen = 'historical'
#scen = 'rcp85'
hd_thres = 'per95' # threshold for defining heatwaves
temp_var = 'td' #dry or wet bulb heatwaves
restrict = str(10.0)

save_path = '/nfs/a321/earsch/floods_heatwaves/processed/combined_events/%s/hotdays/%s/%s/restrict/partial_files/' % (loc, temp_var, hd_thres)  


##flood events
#1 = flood event
fname = '/nfs/a321/earsch/floods_heatwaves/processed/wap/%s/floods/restrict/floods_%s_%s_r%s.nc' % (loc, mod, scen, restrict)
floods = iris.load_cube(fname)

fname = '/nfs/a321/earsch/floods_heatwaves/processed/wap/%s/floods/restrict/start_%s_%s_r%s.nc' % (loc, mod, scen, restrict)
f_start = iris.load_cube(fname)

fname = '/nfs/a321/earsch/floods_heatwaves/processed/wap/%s/floods/restrict/end_%s_%s_r%s.nc' % (loc, mod, scen, restrict)
f_end = iris.load_cube(fname)

##hotdays
fname = '/nfs/a321/earsch/floods_heatwaves/processed/hotdays/%s/%s/%s/hotdays_%s_%s.nc' \
% (loc, temp_var, hd_thres, mod, scen)
hw = iris.load_cube(fname)

#start ehatwave
fname = '/nfs/a321/earsch/floods_heatwaves/processed/hotdays/%s/%s/%s/start_%s_%s.nc' \
% (loc, temp_var, hd_thres, mod, scen)
hw_start = iris.load_cube(fname)

#end heatwave
# 1= last day of heatwave (is a heatwave day)
fname = '/nfs/a321/earsch/floods_heatwaves/processed/hotdays/%s/%s/%s/end_%s_%s.nc' \
% (loc, temp_var, hd_thres, mod, scen)
hw_end = iris.load_cube(fname)


#import landsea mask
cs = floods.coord_system(iris.coord_systems.CoordSystem)

ls = iris.load_cube('/nfs/a277/IMPALA/data/4km/ANCILS/landseamask_ancil_4km_regrid.nc')
ls.coord('longitude').points = ls.coord('longitude').points - 360
ls.coord('longitude').guess_bounds()
ls.coord('latitude').guess_bounds()
ls.coord(axis='x').coord_system = cs
ls.coord(axis='y').coord_system = cs
ls_regrid = ls.regrid(floods, iris.analysis.AreaWeighted())
ls_regrid =ls_regrid[0,0]


#%%

def metrics(events, nt):
    '''takes array 1s, 0s, with 1st where combined events are
    Calculates:
        ndays - number of days of combined floods/heatwaves
        nevents - number of individual combined events
        duration - duration of combined events
        
    '''
    event_duration = np.zeros(events.shape)
    start_event_array = np.zeros(events.shape)
    end_event_array = np.zeros(events.shape)
    
    e=np.zeros(nt+2) # add an element at beginning and end compared to events
                
    #find locations of combined evenst
    eidx = np.where(events == 1.0)
    
    e[eidx[0] + 1] = 1 # add 1 to every index after flood defined

    diffs=np.diff(e) 
    # this is 1 where it becomes flood and -1 where it becomes not flood 
    # and 0 where it stays the same
    
    # get indices where changes occur and find length of flood events
    start_event=np.where(diffs>0)[0]
    end_event=np.where(diffs<0)[0] # end flood day - day flood ends (not a flood day)
    nstart=len(start_event)
    #nend=len(end_hot)
    len_event=end_event-start_event
    
   
    #get total num combined event days
    n_event_days = np.sum(len_event)
    count_event_days = n_event_days
    
    #number of individual flood events
    count_events = nstart



    # event_duration:  cube with length of flood at first day of flood
    #  start flood
    #  end flood
    for i in range(nstart):

        #flood duration
        event_duration[start_event[i]]=len_event[i]
        
        #assign 1s to start and end dates of combined events
        start_event_array[start_event[i]] = 1
        if end_event[i] < nt: # otherwise outside array, a dn event ends outsied time frame
            end_event_array[end_event[i]] = 1 # last day of combined event 
                                                            # is an event day day
                     

    return [count_event_days, count_events, event_duration, start_event_array, end_event_array]


def create_save_cubes(incube, n2d, n3d):
        '''
        Set up cubes for saving output, based on dimensions incube
        n2d = number of 2d cubes (lat, lon) cubes desired
        n3d = number of 3d cubes (timelat, lon) cubes desired
        '''
        shape=np.shape(incube)
        nt=shape[0]
        nlat=shape[1]
        nlon=shape[2]
        
    
        output_list =  iris.cube.CubeList()
        for i in np.arange(n2d):
            cube = copy.deepcopy(incube[0])
            cube.data = np.zeros((nlat,nlon))
            output_list.append(cube)
            
        
        for i in np.arange(n3d):
            cube = copy.deepcopy(incube)
            cube.data = np.zeros((nt, nlat,nlon))
            output_list.append(cube)
        return output_list
    
def new_hw_locs(hws, ndays, shift_dir):
    ''' 
    Will create new array with 1s for each ndays before or after
    start/end heatwave
    
    hws = array with 1/0 for start of ehatwaves (or even)
    ndays = number of days before/after heatwave
    
    if looking for events before the start of heatwaves
    shift_dir = -1.0 (moving array towards left)
    if lookign fore events after the end of heatwaves
    shift_dir = 1.0 (moving array towards right)
    '''

    # shift start back by ndays - will create new array for each day before heatwave event,
    # with 1 for each day before the event, up to ndays
    new_list = []
    for z in np.arange(1, ndays +1):
        new_arr = np.roll(hws, z * shift_dir)
        new_list.append(new_arr)
    
    # add arrays together, so get one array, with 1s for each day before heatwve event 
    # for ndays
    new = new_list[0]
    for z in np.arange(1, len(new_list)):
        new = new + new_list[z]
        
    #make sure only 1s and zeros
    #possible might get ovlerapping heatawves, and when add get a 2 or more
    new = np.where(new >= 1, 1, 0)
    
    return new

#%%

def get_combined(floods, f_start, f_end, hw, hw_start, hw_end, ls, ndays = [3, 5, 7]):
    ''' Find where floods occur - during a heatwave, afer or before a heatwave
    
        ls = landsea mask
        floods = 0,1 flood event cube
        hw_start = 0,1 first day of heatwave cube
        hw_end = 0,1, last day of ehatwave cube
        
        ndays =  number days ebfore and after heatwave to look for floods
       
    '''
    
    ## Create cubes to save output data
    
    shape=np.shape(floods)
    nt=shape[0]
    nlat=shape[1]
    nlon=shape[2]

     
    # same day events
    sameday = copy.deepcopy(floods)
    sameday.data = np.zeros((nt, nlat,nlon))
    #output for saving savemday metrics
    sameday_output = create_save_cubes(floods, 2, 3)
    
    
    #floods before heatwaves - list, with cube for each ndays
    fbefore = iris.cube.CubeList()
    fbefore_output = iris.cube.CubeList()
    for i in np.arange(len(ndays)):
        new_cube = copy.deepcopy(floods)
        new_cube.data = np.zeros((nt, nlat,nlon))
        fbefore.append(new_cube)
        
        new_cubelist = create_save_cubes(floods, 2, 3) 
        fbefore_output.append(new_cubelist)
    
    #flodos after heatwaves
    fafter = iris.cube.CubeList()
    fafter_output = iris.cube.CubeList()
    for i in np.arange(len(ndays)):
        new_cube = copy.deepcopy(floods)
        new_cube.data = np.zeros((nt, nlat,nlon))
        fafter.append(new_cube)
        
        new_cubelist = create_save_cubes(floods, 2, 3) 
        fafter_output.append(new_cubelist)
       
    
    #find sameday
    sameday.data = floods.data + hw.data
    sameday.data = np.where(sameday.data == 2.0, 1.0, 0.0)
    
    for y in range(nlat):
        for x in range(nlon):
            if ls[y,x].data > 0.5:
                
                f_start_loc = f_start.data[:,y,x]
                f_end_loc = f_end.data[:,y,x]
                
                hw_start_loc = hw_start.data[:,y,x]
                hw_end_loc = hw_end.data[:,y,x]
                
                ## Calculate sameday metrics
                sameday_metrics = metrics(sameday[:,y,x].data, nt)
                #assign metrics to save cubes
                #order: count_event_days, count_events, event_duration, start_event_array, end_event_array
                for i in np.arange(len(sameday_metrics)):
                    if i < 2:
                        sameday_output[i].data[y,x] = sameday_metrics[i]
                    else:
                        sameday_output[i].data[:,y,x] = sameday_metrics[i]
                
                
                ## Find flood events occur before/after heatawve
                for z in np.arange(len(ndays)):
                    #flood events before heatwave
                    #end of flood within ndays of start of heatwve
                    #will not include sameday events (or shouldn't)
                    n =  ndays[z]
                    daysbefore = new_hw_locs(hw_start_loc, n, -1)
                    
                    f_before_hw = daysbefore + f_end_loc
                    f_before_hw = np.where(f_before_hw == 2.0, 1.0, 0.0)
                    
                    fbefore[z].data[:,y,x] = f_before_hw
                    
                    #flodo events after heatave
                    #start of flood within ndays of end of heatwve
                    #will not include sameday events (or shouldn't)
                    n =  ndays[z]
                    daysafter = new_hw_locs(hw_end_loc, n, 1)
                    
                    f_after_hw = daysafter + f_start_loc
                    f_after_hw = np.where(f_after_hw == 2.0, 1.0, 0.0)
                    
                    fafter[z].data[:,y,x] = f_after_hw
                
                
                ## Calculate before/after  metrics
                for z in np.arange(len(ndays)):
                    fbefore_metrics = metrics(fbefore[z][:,y,x].data, nt)
                    #assign metrics to save cubes
                    #order: count_event_days, count_events, event_duration, start_event_array, end_event_array
                    for i in np.arange(len(fbefore_metrics)):
                        if i < 2:
                            fbefore_output[z][i].data[y,x] = fbefore_metrics[i]
                        else:
                            fbefore_output[z][i].data[:,y,x] = fbefore_metrics[i]
                            
                    fafter_metrics = metrics(fafter[z][:,y,x].data, nt)
                    #assign metrics to save cubes
                    #order: count_event_days, count_events, event_duration, start_event_array, end_event_array
                    for i in np.arange(len(fafter_metrics)):
                        if i < 2:
                            fafter_output[z][i].data[y,x] = fafter_metrics[i]
                        else:
                            fafter_output[z][i].data[:,y,x] = fafter_metrics[i]
                    
                    

    return [sameday, sameday_output, fbefore, fbefore_output, fafter, fafter_output]



#%%
#calculate combined events

ndays = [1, 2, 3, 4, 5, 6, 7, 15]



#for pan africa split into parts (need to do all timesteps at once for wap due to window)
flood_dims = floods.shape

#print(pr_dims)

start_idx = np.arange(0, flood_dims[1], 50)
end_idx = start_idx + 50
end_idx[-1] = flood_dims[1]

n_parts = len(start_idx)
for k in np.arange(0, n_parts):
    print('Starting combined events part ', k, start_idx[k], end_idx[k])
    new_floods = floods[:, start_idx[k]:end_idx[k], :]
    new_fstart = f_start[:, start_idx[k]:end_idx[k], :]
    new_fend = f_end[:, start_idx[k]:end_idx[k], :]
    new_hw = hw[:, start_idx[k]:end_idx[k], :]
    new_hstart = hw_start[:, start_idx[k]:end_idx[k], :]
    new_hend = hw_end[:, start_idx[k]:end_idx[k], :]
    new_ls = ls_regrid[start_idx[k]:end_idx[k], :]
    
    output = get_combined(new_floods, new_fstart, new_fend, new_hw, new_hstart, new_hend, new_ls, ndays)
            
    print('Saving combined events part ', k)
    
    
    #save data
    var_names = ['ndays', 'nevents', 'duration', 'start', 'end']
    
    
    ##save sameday
    save_name = save_path +  'sameday_combinedevents_' + mod + '_' + scen + '.nc'
    iris.save(output[0], save_name)
    #save saemday metrics
    sameday_mets = output[1]
        
    
    for i in np.arange(len(sameday_mets)):
    
        save_name = save_path +  'sameday_' + var_names[i] + '_' + mod + '_' + scen + '_r' + str(restrict) + '.nc'
    
        iris.save(sameday_mets[i], save_name)
    
    
    # save fbefore/fafter
    fbefore = output[2]
    fafter = output[4]
    
    for n in np.arange(len(ndays)):
        
        save_name = save_path +  'fbefore_' + str(ndays[n]) + '_combinedevents_' + mod + '_' + scen + '_r' + str(restrict) + '.nc'
        iris.save(fbefore[n], save_name)
        save_name = save_path +  'fafter_' + str(ndays[n]) + '_combinedevents_' + mod + '_' + scen + '_r' + str(restrict) + '.nc'
        iris.save(fafter[n], save_name)
    
    #save before/after mets
    fbefore_mets = output[3]
    fafter_mets = output[5]
    for n in np.arange(len(ndays)):
        
        for i in np.arange(len(fbefore_mets[n])):
            save_name = save_path +  'fbefore_' + str(ndays[n]) + '_' + var_names[i] + '_' + mod + '_' + scen + '_r' + str(restrict) + '.nc'
            iris.save(fbefore_mets[n][i], save_name)
            save_name = save_path +  'fafter_' + str(ndays[n]) + '_' + var_names[i] + '_' + mod + '_' + scen + '_r' + str(restrict) + '.nc'
            iris.save(fafter_mets[n][i], save_name)



