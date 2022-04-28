# -*- coding: utf-8 -*-
"""
Master code made for running the full processing chain for all retrackers at all thresholds with all DEMs etc. For the Paper 1. 
CryoSat2 Baseline D data, and ICESat2 ATL06 data. With ArcticDEM and Helm. Thresholds at 10%, 20%, 30%, 40%, 50%.

Created by: Kat Sejan 1st October 2021.
Last edited on: 8th October 2021 by Kat Sejan

"""
#%%
"PART 0 - Data download"

#%%
'Download CryoSat2'

#%%
'Download ICESat2'


#%%
"PART 1 - Data preperation"

#%%
'Prepare DEMs'

"""
1) Mosaic and resample the DEMs in QGIS. Method used in Paper 1:  
2) Compute the slope and aspect rasters in QGIS. Method used in Paper 1:  
3) In python 2.7 run the script: 'dem_preperation.py'

"""

#%%
'Prefilter and tide correct ICESat2'
"EDDIT THE BELOW" "+ add the tide correction code to the prematching code"

import os
from pathlib2 import Path
from kryopy import IS2_pre_matching
import multiprocessing as mp
import gc
import csv

nr_threads = 6
cryo_dir = Path('/Users/kat/DATA/2018_2020_CvsI/CryoSat/')
inice = Path('/Users/kat/DATA/2018_2020_CvsI/ICESAT/')
outdirectory = Path('/Users/kat/DATA/2018_2020_CvsI/ICESAT_usefull/')
buffer_dist = 6000
bucket_size = 8

cryo_file_list = [os.path.join(cryo_dir, file) for file in os.listdir(cryo_dir) if file.endswith('.nc')]
ice_f = [os.path.join(inice,file_name) for file_name in os.listdir(inice) if file_name.endswith('.nc')]

usefull_icefiles = []

running = True

while running:        
    if __name__ == '__main__':
        pool = mp.Pool(processes=nr_threads)
        for cryo_file_path in cryo_file_list:
            result = pool.apply_async(IS2_pre_matching.IS2_parallel_pre_match,  args=(cryo_file_path, ice_f, buffer_dist, bucket_size))
            usefull_icefiles.append(result)
        pool.close()
        pool.join()
    del pool
    gc.collect()


#write names of all matched files into a csv so we don't need to run this again (might be usefull)
uif_list = [p.get() for p in usefull_icefiles]
uif = zip(uif_list)
    
out_name = 'ICESAT2_usefull.csv'
out_data_file = os.path.join(outdirectory,out_name)
    
with open(out_data_file, "w") as f:
   writer = csv.writer(f)
   writer.writerow(["file"])
   for row in uif:
       writer.writerow(row)     


#%%
"PART 2 - Processing"

#%%

'Run the retracker function through the CryoSat2 data'
from kryopy import retrack
import multiprocessing as mp
from pathlib2 import Path
import kryopy
import time
import os

t0 = time.time()

alphain = [0.1, 0.2, 0.3, 0.4, 0.5]

indirectory = Path('/Users/kat/DATA/2018_2020_CvsI/CryoSat/')
outdirectory = '/Users/kat/DATA/2018_2020_CvsI/CryoSat_retrack/CryoSat_ocog/'
#alpha = 0.3
extension = '.nc'
files = [file for file in os.listdir(indirectory) if file.endswith(extension)]
notfiles = [file for file in os.listdir(indirectory) if file.endswith('retracked.nc')]
file_list = [os.path.join(indirectory, file_name) for file_name in files if file_name not in notfiles]


#specify that you want to allow only 6 threads to be used by python (not all 8)
pool = mp.Pool(processes=6)
for file in file_list:
    pool.apply_async(retrack.ocog_retrack, args=(file, outdirectory, alphain))
pool.close()
pool.join()

t1 = time.time()
total_time = (t1-t0)/60
print('All files retracked on 6 CPUs in: ' + str(total_time/60) + ' h') 


#%%

'Run the slope correction function through the CryoSat2 data'
from kryopy import slope
import multiprocessing as mp
from pathlib2 import Path
import time
import os
import shutil
import gc

t0 = time.time()

"Parameters"
nr_threads = 6
directory = Path('/Users/kat/DATA/2018_2020_CvsI/CryoSat_slope/CryoSat_tfmra_beta5/') #run for 3 retrackers, folders: 'CryoSat_ocog', 'CryoSat_tfmra', 'CryoSat_tfmra_beta5'

dems = ['gimp', 'arctic', 'helm']
alpha = [0.1, 0.2, 0.3, 0.4, 0.5]

extension = '.nc'
files = [file for file in os.listdir(directory) if file.endswith(extension)]
file_list = [os.path.join(directory, file_name) for file_name in files]
    
for dem in dems:
    print(dem)
    dem_x, dem_y, elev_rast, slope_rast, aspect_rast = slope.dem_open(dem)

    pool = mp.Pool(processes=nr_threads)
    for i,file in enumerate(file_list):
        pool.apply_async(slope.slope_correction,  args=(file, dem_x, dem_y, elev_rast, slope_rast, aspect_rast, alpha, dem))
    pool.close()
    pool.join()

    
t1 = time.time()
total_time = (t1-t0)/60
print('All files slope corrected on 6 CPUs in with both GIMP DEM and ArcticDEM in: ' + str(total_time) + ' min, or ' + str(total_time/60) + ' hours.') 

#%%
'Matching L1b retracked and slope corrected data'
"EDDIT THE BELOW!!!"

import multiprocessing as mp
from pathlib2 import Path
import kryopy
import time
import os
from kryopy import matching_L1b
import gc

"Parameters"
matching_distance = 1000 #in m 
day_bucket = 5 #in days (number of days before and after the date)
nr_threads = 6
main_directory = Path('/Users/kat/DATA/2018_2020_CvsI/CryoSat_retracked/') #directory in which the inand the out folders are
in_folder_name = 'ocog_slope_0_5_1km'
out_folder_name = 'ocog_matched_0_5_1km_' + str(matching_distance) + 'm'
ice_sat_directory = Path('/Users/kat/DATA/2018_2020_CvsI/ICESAT/')

"Create a new folder for the results of this slope correction and copy over the data so that it can be ammended with new results."
out_directory = Path(os.path.join(main_directory, out_folder_name))
in_directory = Path(os.path.join(main_directory, in_folder_name))


    
#for i,month in enumerate(mlist):
inice = ice_sat_directory
inice_before = ice_sat_directory
inice_after = ice_sat_directory
outdirectory = out_directory
outdirectory.mkdir(parents=True, exist_ok=True)

cryo_files_list = [file_name for file_name in os.listdir(in_directory) if file_name.endswith('.nc')]
cryo_paths_list = [os.path.join(in_directory,f) for f in cryo_files_list] 
#cryo_file_path = cryo_paths_list[2]
#bucket_size = 5
#buffer_dist = 1000

t0 = time.time()    
if __name__ == '__main__':
    #specify that you want to allow only 6 threads to be used by python (not all 8)
    pool = mp.Pool(processes=nr_threads)
    for cryo_file_path in cryo_paths_list:
        pool.apply_async(matching_L1b.cryo_parallel_match, args=(cryo_file_path, inice, inice_before, inice_after, outdirectory, matching_distance, day_bucket))
    pool.close()
    pool.join()
del pool
gc.collect()
t1 = time.time()
total_time = (t1-t0)/60
print('All files matched to ICESat2 on 6 CPUs in: ' + str(total_time) + ' min') 
  

#%%
'Matching L2 data'
"EDDIT THE BELOW!!!"

import multiprocessing as mp
from pathlib2 import Path
import kryopy
import time
import os
from kryopy import matching_L2

"Parameters"
matching_distance = 1000 #in m 
day_bucket = 5 #in days (number of days before and after the date)
nr_threads = 6
main_directory = Path('/Users/kat/DATA/2018_2020_CvsI/') #directory in which the inand the out folders are
in_folder_name = 'CryoSat_L2i'
out_folder_name = 'CryoSat_L2i_matched_' + str(matching_distance) + 'm'
ice_sat_directory = Path('/Users/kat/DATA/2018_2020_CvsI/ICESAT_usefull/')

"Create a new folder for the results of this slope correction and copy over the data so that it can be ammended with new results."
out_directory = Path(os.path.join(main_directory, out_folder_name))
in_directory = Path(os.path.join(main_directory, in_folder_name))



inice = ice_sat_directory
inice_before = ice_sat_directory
inice_after = ice_sat_directory
outdirectory = out_directory
outdirectory.mkdir(parents=True, exist_ok=True)

cryo_files_list = [file_name for file_name in os.listdir(in_directory) if file_name.endswith('.nc')]
cryo_paths_list = [os.path.join(in_directory,f) for f in cryo_files_list] 


t0 = time.time()    
if __name__ == '__main__':
    #specify that you want to allow only 6 threads to be used by python (not all 8)
    pool = mp.Pool(processes=nr_threads)
    for cryo_file_path in cryo_paths_list:
        pool.apply_async(matching_L2.cryo_parallel_match, args=(cryo_file_path, inice, inice_before, inice_after, outdirectory, matching_distance, day_bucket))
    pool.close()
    pool.join()
t1 = time.time()
total_time = (t1-t0)/60
print('All files matched to ICESat2 on 6 CPUs in: ' + str(total_time) + ' min') 


#%%
"THE END"
