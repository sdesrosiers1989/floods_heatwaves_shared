#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''Data Preparation  - Put togethet pan-africa flood files
Date: 1/03/2022
Author: Sarah C

            
'''

#%% Import packages

import matplotlib.pyplot as plt 

import iris
import iris.cube
from cf_units import Unit

from iris.util import equalise_attributes
import iris.coord_categorisation


import glob




#%% Import data
#update model name and scenario name in filenames, and saved path at end of script


mod = 'p25'
scen = 'historical'
scen = 'rcp85'
temp_var = 'tw' 
loc = 'pan_africa'
hd_thres = 'per95' # threshold for defining heatwaves
#dry or wet bulb heatwaves
restrict = str(10.0)

file_path = '/nfs/a321/earsch/floods_heatwaves/processed/combined_events/%s/hotdays/%s/%s/restrict/partial_files/' % (loc, temp_var, hd_thres)  
save_path = '/nfs/a321/earsch/floods_heatwaves/processed/combined_events/%s/hotdays/%s/%s/restrict/' % (loc, temp_var, hd_thres)  


var_names = ['combinedevents', 'ndays', 'nevents', 'duration', 'start', 'end']

cb_timing = ['sameday', 'fbefore', 'fafter']

ndays_list = [1, 2, 3, 4, 5, 6, 7, 15]

def concat_files(model, scen, file_path, save_path, var, timing, ndays):

    m_scen = model + '_' + scen 
    tim_days_var = timing + '_' + str(ndays) + '_' + var
    
    if timing == 'sameday':
        tim_days_var = timing + '_' + var

    filenames = glob.glob(file_path + tim_days_var + '_{}*part?.nc'.format(m_scen))
 
    files = iris.cube.CubeList()
    for file in filenames:
        x = iris.load_cube(file)
        files.append(x)    
    
    print('Loading:')
    for f in filenames:
        print(f)

    equalise_attributes(files)
    
    for cube in files:
        cube.coord('model').units = Unit('unknown')
        cube.coord('sim').units = Unit('unknown')
    
    concat_files = files.concatenate_cube()
    
    

    if timing == 'sameday':
        save_name = save_path + timing + '_' + var + '_' + model + '_' + scen + '.nc'
    else:
        save_name = save_path + timing + '_' + str(ndays) + '_' + var + '_' + model + '_' + scen + '.nc'
    print('Saving as ', save_name)
    
    iris.save(concat_files, save_name)


#%%

for timing in cb_timing:
    for var in var_names:
        if timing == 'sameday':
            concat_files(mod, scen, file_path, save_path, var, timing, 1)
        else:
            for nd in ndays_list:
                concat_files(mod, scen, file_path, save_path, var, timing, nd)