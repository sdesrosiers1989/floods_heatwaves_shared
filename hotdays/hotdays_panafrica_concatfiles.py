#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''Data Preparation  - Put togethet pan-africa flood files
Date: 1/03/2022
Author: Sarah C

            
'''

#%% Import packages

import matplotlib.pyplot as plt 

import iris

from iris.util import equalise_attributes
import iris.coord_categorisation


import glob




#%% Import data
#update model name and scenario name in filenames, and saved path at end of script

model = 'cp4'
scen = 'rcp85'
temp_type = 'td'
#scen = 'historical'



var_names =  ['hotdays', 'duration', 'ndays', 'nevents', 'start', 'end']

file_path = '/nfs/a321/earsch/floods_heatwaves/processed/hotdays/pan_africa/' + temp_type + '/per95/partial_files/'
save_path = '/nfs/a321/earsch/floods_heatwaves/processed/hotdays/pan_africa/' + temp_type + '/per95/'

def concat_files(model, scen, file_path, save_path, var):

    m_scen = model + '_' + scen 

    filenames = glob.glob(file_path + var + '*{}*part?.nc'.format(m_scen))
    files = iris.cube.CubeList()
    for file in filenames:
        x = iris.load_cube(file)
        files.append(x)    
    
    print('Loading:')
    for f in filenames:
        print(f)

    equalise_attributes(files)
    concat_files = files.concatenate_cube()


    save_name = save_path + var + '_' + model + '_' + scen + '.nc'
    
    print('Saving as ', save_name)
    
    iris.save(concat_files, save_name)


#%%
for var in var_names:
    concat_files(model, scen, file_path, save_path, var)
