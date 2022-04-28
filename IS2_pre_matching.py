"This code is for reducing the number of ICESat2 files by prematching to CryoSat2."



def track_poly(file_points_list, buffer_dist):
    
    from scipy.spatial import ConvexHull
    from shapely import geometry
    import numpy as np
    
    #Calculate the convex hull
    hull = ConvexHull(file_points_list) #in the code do this for cryo_points_gimp_scipy and cryo_points_arctic_scipy
    vertices = hull.vertices
    polygon_points = []
    for vertice in vertices:
        polygon_points.append(file_points_list[vertice])
    polygon_points = np.array(polygon_points)

    #build a polygon from the points including a buffer
    #build polygon
    poly = geometry.Polygon([[p[0], p[1]] for p in polygon_points])
    
    #add buffer
    poly_buff = geometry.Polygon(poly.buffer(buffer_dist).exterior)
    
    
    return poly_buff


#%%
import numpy as np
import os
from datetime import datetime, timedelta
from netCDF4 import Dataset
from pathlib import Path
from pyproj import Transformer
import csv
import pandas as pd
import shutil

cryo_dir = Path('/Users/kat/DATA/2018_2020_CvsI/CryoSat/')
inice = Path('/Users/kat/DATA/2018_2020_CvsI/ICESAT/')
outdirectory = Path('/Users/kat/DATA/2018_2020_CvsI/')
buffer_dist = 6000
bucket_size = 8


#%%
from shapely.geometry import box 
in_path = '/Users/kat/DATA/2018_2020_CvsI/ICESAT/ICESAT2 Files bounds.csv'

data = pd.read_table(in_path, delimiter=',', names=["file", "minx", "miny", "maxx", "maxy"], dtype=str, error_bad_lines=False, header=1)

fbox = []
fileb = []
for fi,ff in enumerate(data['file']):
    fbox.append(box(float(data['minx'][fi]), float(data['miny'][fi]), float(data['maxx'][fi]), float(data['maxy'][fi])))
    fileb.append(ff)

box_list = [fileb, fbox]
box_data = pd.concat([pd.Series(x) for x in box_list], axis=1)
df_fbox = pd.DataFrame(box_data.values, columns = ['file', 'box'])


#%%

cryo_file_list = [os.path.join(cryo_dir, file) for file in os.listdir(cryo_dir) if file.endswith('.nc')]
ice_f = [file_name for file_name in os.listdir(inice) if file_name.endswith('.nc')]

usefull_icefiles = []

for iic,cryo_file_path in enumerate(cryo_file_list):
    print('Checking for file: ' + str(iic+1) + '/' + str(len(cryo_file_list)))
    
    "Extract date and make a time bucket"
    cfile_end = cryo_file_path.split('1B_')[1]
    file_out_name = cfile_end[0:31:1]
    file_date = cfile_end[0:8:1]
    cdate = datetime.strptime(file_date, '%Y%m%d')
    year = datetime.strptime(file_date, '%Y%m%d').year
    month = datetime.strptime(file_date, '%Y%m%d').month
    
    bucket_dates = []
    bucket_dates.append(cdate)
    for n in range(1,bucket_size+1):
        bucket_dates.append(cdate - timedelta(days=n))
        bucket_dates.append(cdate + timedelta(days=n))
    bucket_dates = sorted(bucket_dates)
    
    month_dates = [x for x in bucket_dates if x.month == month]
    
    ifile_month = [d.strftime('%Y%m%d') for d in month_dates]
    
    ice_files_list = []
    for s in ifile_month:
        ice_files_list.append([file_name for file_name in ice_f if s in file_name])

    ice_files = [x for sub in ice_files_list for x in sub]
    
      
    "Read in CryoSat-2 file data"
    data_file = Dataset(cryo_file_path, "r", format="NETCDF4")

    clon = data_file.variables['lon_20_ku'][:].tolist()
    clat= data_file.variables['lat_20_ku'][:].tolist()
        
    data_file.close()     
    
    #Converting the wavelet location to polar stereographic projection to match the DEM extend
    inProj = 'epsg:4326'
    outProj = 'epsg:3413'
    polar_transformer = Transformer.from_crs(inProj, outProj)
    #below hashtaged seems to be depricated: 
    #inProj = Proj(init='epsg:4326')
    #outProj = Proj(init='epsg:3413')
    cx = []
    cy = []
    nr_points = len(clon)
    for i3 in range(nr_points):
        x_trans, y_trans  = polar_transformer.transform(clat[i3],clon[i3])
        #x_trans, y_trans = transform(inProj,outProj,nadir_long[i3],nadir_lat[i3])
        cx.append(x_trans)
        cy.append(y_trans)

    
    cryo_points = []
    for i,p in enumerate(cx):
        cryo_points.append([cx[i], cy[i]])
    
    cryo_points = np.array(cryo_points)
    
    "Create CryoSat-2 polygon"
    cryo_poly = track_poly(cryo_points, buffer_dist)


    for indx, path_to_ifile in enumerate(ice_files):
        bby = df_fbox.loc[df_fbox['file']==path_to_ifile]['box'].values
        if len(bby) != 0:
            boxby = bby[0]
        if cryo_poly.intersects(boxby) == True:
            usefull_icefiles.append(path_to_ifile)


#write names of all matched files into a csv so we don't need to run this again (might be usefull)
uif = zip(list(set(usefull_icefiles)))
    
out_name = 'ICESAT2_usefull.csv'
out_data_file = os.path.join(outdirectory,out_name)
    
with open(out_data_file, "w") as f:
   writer = csv.writer(f)
   writer.writerow(["file"])
   for row in uif:
       writer.writerow(row)     

#%%
files_keep = list(set(usefull_icefiles))

files_delete = [file.replace('.nc', '') for file in ice_f if file not in files_keep]

odir = Path('/Users/kat/DATA/2018_2020_CvsI/ICESAT_usefull/')

for ik, fike in enumerate(files_keep):
    dir_keep = os.path.join(inice,fike)
    new_dir = os.path.join(odir,fike)
    shutil.copy(dir_keep, new_dir)
    
