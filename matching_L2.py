# -*- coding: utf-8 -*-
"""
Matching CryoSat-2 to ICESat2 points with kdtree method.

Created by: Kat Sejan 30th March 2020.

"""

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



def intersect_matches_kdtree(cryo_points, icesat_points, cryo_elev, ice_elev, cryo_date, ice_date, longitudes, latitudes, buffer_dist):
    #cryo_poly must be either cryo_poly_gimp or cryo_poly_arctic
    
    
    from scipy import spatial
    import numpy as np


    #using double tree query method to find all points within buffer_dist:
    
    ice_tree = spatial.cKDTree(np.array(icesat_points), leafsize=32, balanced_tree=False)
    
    cryo_tree = spatial.cKDTree(np.array(cryo_points), leafsize=32, balanced_tree=False)

    ice_match_index = cryo_tree.query_ball_tree(ice_tree, buffer_dist)

    
    #average out the resulting icesta elevation:

    #cryo points lists
    cryo_elev_matched = [cryo_elev[indxi] for indxi, l in enumerate(ice_match_index) if l] #output for no matches is empy list, if l is used to find only True lists i.e. full lists
    #len(cryo_elev_matched)
    cryo_x_matched = [cryo_points[indxi][0] for indxi, l in enumerate(ice_match_index) if l] 
    #len(cryo_x_matched)
    cryo_y_matched = [cryo_points[indxi][1] for indxi, l in enumerate(ice_match_index) if l]
    #len(cryo_x_matched)
    cryo_date_matched = [cryo_date for i in cryo_x_matched]
    longitude_matched = [longitudes[indx] for indx, l in enumerate(ice_match_index) if l]
    latitude_matched = [latitudes[indx] for indx, l in enumerate(ice_match_index) if l]
    
    #ice
    ice_elev_matched = []
    ice_date_matched = []
    ice_x_m = []
    ice_y_m = []
    for l in ice_match_index:
        if l:
            ice_indices_mached = l
            ice_elevs_m = [ice_elev[i] for i in ice_indices_mached]
            ice_elevs_m = [np.nan if x > 3.0e+30 else x for x in ice_elevs_m] #usually icesat fill value for elevation is 3.4e+38
            ice_elev_matched.append(np.nanmean(ice_elevs_m))
            ice_date_matched.append(ice_date)
            ice_x = np.nanmean([icesat_points[i][0] for i in ice_indices_mached])
            ice_y = np.nanmean([icesat_points[i][1] for i in ice_indices_mached])
            ice_x_m.append(ice_x)
            ice_y_m.append(ice_y)

    return  cryo_elev_matched, cryo_x_matched, cryo_y_matched, longitude_matched, latitude_matched, cryo_date_matched, ice_elev_matched, ice_date_matched, ice_x_m, ice_y_m,

def points_greenland(longitudes, latitudes):
    'This function checks all of the longs and lats against the Greenland shape file, and if points are not in Greenland they are not retracked.'
    
    from shapely.geometry import Point
    from pathlib import Path
    from cartopy.io.shapereader import Reader as Reader
    
    """Reading in the rougth outline shapefile of Greenland."""
    shfile_path = str(Path('/Users/kat/DATA/Greenland_shapefile/Greenland_outline_buffed.shp'))
    greenland_shape = Reader(shfile_path).records()
    
    """Extracting the polygon of Greenland from the shapefile using cartopy. 
    The polygon is represented below as the geometry."""
    greenland = next(greenland_shape)
    g_geometry = greenland.geometry

    
    in_index = []
    out_index = []
    
    for i,l in enumerate(longitudes):
               
        llpoint = Point(longitudes[i], latitudes[i])
        
        #if the line matches the shapefile do nothing, if it does not match then remove the file:
        if g_geometry.contains(llpoint) == True: #checking if the file matches the shapefile
            in_index.append(i)
        else:
            out_index.append(i)
    
    return in_index, out_index 


 
def cryo_parallel_match(cryo_file_path, inice, inice_before, inice_after, outdirectory, buffer_dist, bucket_size):
    #inice_before, inice_after, are directories to data from month before and after the analyesed month
    
    
    import numpy as np
    import os
    from datetime import datetime, timedelta
    from netCDF4 import Dataset
    from pyproj import Transformer
        
    "Paths"
    month_ice_files = [file_name for file_name in os.listdir(inice) if file_name.endswith('.nc')]
    before_ice_files = [file_name for file_name in os.listdir(inice_before) if file_name.endswith('.nc')]
    after_ice_files = [file_name for file_name in os.listdir(inice_after) if file_name.endswith('.nc')]

    "ICESat2 file list"
    print('Load paths of ICESat2 data that match the bucket range of the CryoSat-2 file:' + str(cryo_file_path))
    
    cfile_end = cryo_file_path.split('2__')[1]
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
    
    before_month_dates = [x for x in bucket_dates if x.month < month]
    after_month_dates = [x for x in bucket_dates if x.month > month]
    month_dates = [x for x in bucket_dates if x.month == month]
    
    ifile_month = [d.strftime('%Y%m%d') for d in month_dates]
    ifile_before = [d.strftime('%Y%m%d') for d in before_month_dates]
    ifile_after = [d.strftime('%Y%m%d') for d in after_month_dates]
    
    ice_files_list = []
    for s in ifile_before:
        ice_files_list.append([os.path.join(inice_before,file_name) for file_name in before_ice_files if s in file_name])
    for s in ifile_month:
        ice_files_list.append([os.path.join(inice,file_name) for file_name in month_ice_files if s in file_name])
    for s in ifile_after:
        ice_files_list.append([os.path.join(inice_after,file_name) for file_name in after_ice_files if s in file_name])

    ice_files = [x for sub in ice_files_list for x in sub]
    
    "Read in CryoSat-2 file data"
    print('Create polygon of the CryoSat-2 file being matched.')
    data_file = Dataset(cryo_file_path, "r", format="NETCDF4")

    ###
    x_l = data_file.variables['lon_poca_20_ku'][:]
    y_l = data_file.variables['lat_poca_20_ku'][:]
    elevation = data_file.variables['height_1_20_ku'][:] #OCOG retracker and slope corrected is _3_ ; Land retracker is _2_ ; Ocean retracker is _1_ ;
    
    longitudes = data_file.variables['lon_20_ku'][:]#read in the longitude values for the points
    latitudes = data_file.variables['lat_20_ku'][:] #read in the latitude values for the points
    
    data_file.close()     
    
    #in_index, out_index = points_greenland(longitudes, latitudes)
    
    #x_l = [x_l[il] for il,ll in enumerate(x_l) if il in in_index]
    #y_l = [y_l[il] for il,ll in enumerate(y_l) if il in in_index]
    #elevation = [elevation[il] for il,ll in enumerate(elevation) if il in in_index]
    
    if len(x_l) > 0:
        
        ###Convert the long and lat to polar coordinates    
        inProj = 'epsg:4326'
        outProj = 'epsg:3413'
        polar_transformer = Transformer.from_crs(inProj, outProj)
        #below hashtaged seems to be depricated: 
        #inProj = Proj(init='epsg:4326')
        #outProj = Proj(init='epsg:3413')
        x_p = []
        y_p = []
        nr_points = len(x_l)
        for i in range(nr_points):
            x_trans, y_trans  = polar_transformer.transform(y_l[i],x_l[i])
            #x_trans, y_trans  = polar_transformer.transform(nadir_lat[i3],nadir_long[i3])
            #x_trans, y_trans = transform(inProj,outProj,nadir_long[i3],nadir_lat[i3])
            x_p.append(x_trans)
            y_p.append(y_trans)
    
        
        
        ###Create geo points from polar x and polar y
        cryo_points = []
        for i,p in enumerate(x_p):
            cryo_points.append([x_p[i], y_p[i]])
    
    
        "Create CryoSat-2 polygon"
        cryo_poly = track_poly(cryo_points, buffer_dist)
        
        "Run through ICESat2 files, create polygons and check for intersection and find matching points."
        
        cryo_long= []
        cryo_lat = []
        cryo_elev = []
        lon = []
        lat = []
        cryo_date = []
        ice_elev = []
        ice_date = []
        ice_x = []
        ice_y = []
        
        for indx, path_to_ifile in enumerate(ice_files):
            #read in data points
            "This code is for testing only"
    # =============================================================================
    #         indx = 38
    #         path_to_ifile = ice_files[indx]
    #         "End of testing code"
    # =============================================================================
            
            icesat_file = Dataset(path_to_ifile, "r", format="NETCDF4")        
    
            x_icesat = icesat_file.variables['polar_x'][:].tolist()  
            y_icesat = icesat_file.variables['polar_y'][:].tolist() 
         
            elev_icesat = icesat_file.variables['elevation'][:].tolist() 
            
            file_end = path_to_ifile.split('ATL06_')[1]
            idate = file_end[0:8:1]
            
            icesat_file.close()
            
            icesat_points = []
            for i,p in enumerate(x_icesat):
                icesat_points.append([x_icesat[i],y_icesat[i]])
            
            #Create polygon
            ice_poly = track_poly(icesat_points, buffer_dist)
            
            #Crude check if polygon within min/max of other polygon
            if cryo_poly.intersects(ice_poly) == False:
                #print('Does not intersect for GIMP.')
                pass
            else: 
                #print (str(indx))
                cryo_elev_matched, cryo_x_matched, cryo_y_matched, longitude_matched, latitude_matched, cryo_date_matched, ice_elev_matched, ice_date_matched, ice_x_matched, ice_y_matched = intersect_matches_kdtree(cryo_points, icesat_points, elevation, elev_icesat, file_date, idate, longitudes, latitudes, buffer_dist) 
                                        
    
                cryo_long.append(cryo_x_matched)
                cryo_lat.append(cryo_y_matched)
                cryo_elev.append(cryo_elev_matched)
                lon.append(longitude_matched)
                lat.append(latitude_matched)
                cryo_date.append(cryo_date_matched)
                ice_elev.append(ice_elev_matched)
                ice_date.append(cryo_date_matched)
                ice_x.append(ice_x_matched)
                ice_y.append(ice_y_matched)
                
                
    
        "Prepare output for saving"   
        long_cryo = [np.float32(x) for sub in cryo_long for x in sub]
        lat_cryo = [np.float32(x) for sub in cryo_lat for x in sub]
        elev_cryo = [np.float32(x) for sub in cryo_elev for x in sub]
        lon_c = [np.float32(x) for sub in lon for x in sub]
        lat_c = [np.float32(x) for sub in lat for x in sub]
        date_cryo = [x for sub in cryo_date for x in sub]
        elev_ice = [np.float32(x) for sub in ice_elev for x in sub]
        date_ice = [x for sub in ice_date for x in sub]
        long_ice = [np.float32(x) for sub in ice_x for x in sub]
        lat_ice = [np.float32(x) for sub in ice_y for x in sub]
    
    
        dcg = [int(i) for i in date_cryo]
        dig = [int(i) for i in date_ice]
    
    
        length_output = len(elev_cryo)
    
        
        
        while length_output!=0: 
            print('All the matches found! Saving file ...')    
            "Save this matched points in an output file"
            #save the data into a .nc format:
            #print('Saving the data to a .nc file...')        
            out_name = str(file_out_name) + '_Month_' + str(month) + '_' + str(year) + '_' + str(bucket_size) + 'day_bucket.nc'
            out_data_file = os.path.join(outdirectory,out_name)
            out_file = Dataset(out_data_file, 'w', format='NETCDF4_CLASSIC')
            
            #creating dimentions for the data:
            points = out_file.createDimension('points', length_output)
            
            #assigning variables for the data with above dimentions:
            date_cryo = out_file.createVariable('date_cryo', np.int32, ('points',))
            x_cryo = out_file.createVariable('px_cryo', np.float32, ('points',))
            y_cryo = out_file.createVariable('py_cryo', np.float32, ('points',))
            e_cryo = out_file.createVariable('elevation_cryo', np.float32, ('points',))
            lon_orig = out_file.createVariable('original_longitude', np.float32, ('points',))
            lat_orig = out_file.createVariable('original_latitude', np.float32, ('points',))
            date_ice = out_file.createVariable('date_ice', np.int32, ('points',))
            e_ice = out_file.createVariable('elevation_ice', np.float32, ('points',))
            x_ice = out_file.createVariable('px_ice', np.float32, ('points',))
            y_ice = out_file.createVariable('py_ice', np.float32, ('points',))
    
                   
            #assign atributes to the file (global) and to the variables (variable):
            out_file.description = 'This file contains matched CryoSat-2 data points and IceSat-2 data points. Coordinates have been convarted to stereographic polar projection from the decimal degree latitude and longitude.'
            out_file.source = 'Original data from an IceSat-2 (NASA). CryoSat-2 data computed witha retracker by KM Sejan 2019.'
            
            date_cryo.units = 'yyyymmdd'
            x_cryo.units = 'meters'
            y_cryo.units = 'meters'
            e_cryo.units = 'meters'
            lon_orig.units = 'decimal degrees'
            lat_orig.units = 'decimal degrees'
            date_ice.units = 'yyyymmdd'
            e_ice.units = 'meters'
            x_ice.units = 'meters'
            y_ice.units = 'meters'
    
    
            date_cryo[:] = dcg
            x_cryo[:] = long_cryo
            y_cryo[:] = lat_cryo
            e_cryo[:] = elev_cryo
            lon_orig[:] = lon_c
            lat_orig[:] = lat_c
            date_ice[:] = dig
            e_ice[:] = elev_ice
            x_ice[:] = long_ice
            y_ice[:] = lat_ice
    
        
          
            out_file.close()
        
            print('...SAVED.')
            break
        print('Month ' + str(month) + ' ' + str(year) + ' has ' + str(length_output) + ' matches.')


