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



def intersect_matches_kdtree(cryo_points, icesat_points, cryo_elev, ice_elev, cryo_date, ice_date, buffer_dist, acdc_flag, le_width, te_slope, backscatter, slope):
    #cryo_poly must be either cryo_poly_gimp or cryo_poly_arctic
    
    "This below is for testing"
    # cryo_points = cryo_points_gimp
    # cryo_elev = elevation_gimp
    # ice_elev = elev_icesat
    # cryo_date = file_date
    # ice_date = idate
    # slope = slope_gimp
    #buffer_dist = 500
    "End of testing code"
    
    from scipy import spatial
    import numpy as np


    #using double tree query method to find all points within buffer_dist:
    
    ice_tree = spatial.cKDTree(np.array(icesat_points), leafsize=32, balanced_tree=False)
    
    cryo_tree = spatial.cKDTree(np.array(cryo_points), leafsize=32, balanced_tree=False)

    ice_match_index = cryo_tree.query_ball_tree(ice_tree, buffer_dist)

    
    #average out the resulting icesta elevation:

    #cryo points lists
    cryo_elev_matched = [cryo_elev[indx] for indx, l in enumerate(ice_match_index) if l] #output for no matches is empy list, if l is used to find only True lists i.e. full lists
    #len(cryo_elev_matched)
    cryo_x_matched = [cryo_points[indx][0] for indx, l in enumerate(ice_match_index) if l] 
    #len(cryo_x_matched)
    cryo_y_matched = [cryo_points[indx][1] for indx, l in enumerate(ice_match_index) if l]
    #len(cryo_x_matched)
    cryo_date_matched = [cryo_date for i in cryo_x_matched]
    acdc_flag_matched = [acdc_flag[indx] for indx, l in enumerate(ice_match_index) if l]
    #len(cryo_x_matched)
    le_width_matched = [le_width[indx] for indx, l in enumerate(ice_match_index) if l]
    te_slope_matched = [te_slope[indx] for indx, l in enumerate(ice_match_index) if l]
    backscatter_matched =[backscatter[indx] for indx, l in enumerate(ice_match_index) if l]
    slope_matched =[slope[indx] for indx, l in enumerate(ice_match_index) if l]
    
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
    
    
    return  cryo_elev_matched, cryo_x_matched, cryo_y_matched, cryo_date_matched, ice_elev_matched, ice_x_m, ice_y_m, ice_date_matched, acdc_flag_matched, le_width_matched, te_slope_matched, backscatter_matched, slope_matched


 
def cryo_parallel_match(cryo_file_path, inice, inice_before, inice_after, outdirectory, buffer_dist, bucket_size):
    #inice_before, inice_after, are directories to data from month before and after the analyesed month
    
    
    import numpy as np
    import os
    from datetime import datetime, timedelta
    from netCDF4 import Dataset

        
    "Paths"
    month_ice_files = [file_name for file_name in os.listdir(inice) if file_name.endswith('.nc')]
    before_ice_files = [file_name for file_name in os.listdir(inice_before) if file_name.endswith('.nc')]
    after_ice_files = [file_name for file_name in os.listdir(inice_after) if file_name.endswith('.nc')]

    "ICESat2 file list"
    print('Load paths of ICESat2 data that match the bucket range of the CryoSat-2 file:' + str(cryo_file_path))
    
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
    
    before_month_dates = [x for x in bucket_dates if x.month < month]
    after_month_dates = [x for x in bucket_dates if x.month > month]
    month_dates = [x for x in bucket_dates if x.month == month]
    
    ifile_month = ['_'+d.strftime('%Y%m%d') for d in month_dates]
    ifile_before = ['_'+d.strftime('%Y%m%d') for d in before_month_dates]
    ifile_after = ['_'+d.strftime('%Y%m%d') for d in after_month_dates]
    
    ice_files_list = []
    for s in ifile_before:
        ice_files_list.append([os.path.join(inice_before,file_name) for file_name in before_ice_files if s in file_name])
    for s in ifile_month:
        ice_files_list.append([os.path.join(inice,file_name) for file_name in month_ice_files if s in file_name])
    for s in ifile_after:
        ice_files_list.append([os.path.join(inice_after,file_name) for file_name in after_ice_files if s in file_name])

    ice_files = [x for sub in ice_files_list for x in sub]
    
    
   
    "Read in CryoSat-2 file data"
    """
    1) unhashtag the arctic dem
    2) make sure variable names are correct
    3) add other needed variables (wavelet, all the betas, all the coordinates etc.) to read in and loop through
    4) save the other variables to new file
    
    """
    
    print('Create polygon of the CryoSat-2 file being matched.')
    data_file = Dataset(cryo_file_path, "r", format="NETCDF4")

    ###GIMP
    x_p_gimp = data_file.variables['ecc_x_' + str(threshold) + '_gimp'][:].tolist()
    y_p_gimp = data_file.variables['ecc_y_' + str(threshold) + '_gimp'][:].tolist()
    print('number of gimp points: ' + str(len(x_p_gimp)))
    elevation_gimp = data_file.variables['ecc_elevation_' + str(threshold) + '_gimp'][:].tolist()

    dem_elev_gimp = data_file.variables['nearest_elevation' + '_gimp'][:].tolist()
    aspect_gimp = data_file.variables['nearest_aspect_radian' + '_gimp'][:].tolist()
    slope_gimp = data_file.variables['nearest_slope_radian' + '_gimp'][:].tolist()
    
    # ###ArcticDEM
    x_p_arctic = data_file.variables['ecc_x_' + str(threshold) + '_arctic'][:].tolist() 
    y_p_arctic = data_file.variables['ecc_y_' + str(threshold) + '_arctic'][:].tolist() 
    print('number of arctic points: ' + str(len(x_p_arctic)))
    elevation_arctic = data_file.variables['ecc_elevation_' + str(threshold) + '_arctic'][:].tolist() 
    
    dem_elev_arctic = data_file.variables['nearest_elevation' + '_arctic'][:].tolist()
    aspect_arctic = data_file.variables['nearest_aspect_radian' + '_arctic'][:].tolist()
    slope_arctic = data_file.variables['nearest_slope_radian' + '_arctic'][:].tolist()

    ###Helm DEM
    x_p_helm = data_file.variables['ecc_x_' + str(threshold) + '_helm'][:].tolist() 
    y_p_helm = data_file.variables['ecc_y_' + str(threshold) + '_helm'][:].tolist() 
    print('number of helm points: ' + str(len(x_p_helm)))
    elevation_helm = data_file.variables['ecc_elevation_' + str(threshold) + '_helm'][:].tolist() 
    
    dem_elev_helm = data_file.variables['nearest_elevation' + '_helm'][:].tolist()
    aspect_helm = data_file.variables['nearest_aspect_radian' + '_helm'][:].tolist()
    slope_helm = data_file.variables['nearest_slope_radian' + '_helm'][:].tolist()

    
    ###Common parameters
    longitudes = data_file.variables['lon_20_ku'][:].astype('float32').tolist() #read in the longitude values for the points
    latitudes = data_file.variables['lat_20_ku'][:].astype('float32').tolist() #read in the latitude values for the points
    power_waveform = data_file.variables['pwr_waveform_20_ku'][:]  #ectract the waveforms for each measurement
    
    max_power_wavelet = data_file.variables['max_power_wavelet'][:].tolist()
    mean_power_wavelet = data_file.variables['mean_power_wavelet'][:].tolist()
    noise_wavelet = data_file.variables['noise_wavelet'][:].tolist()
    signal_to_noise = data_file.variables['signal_to_noise'][:].tolist()
    zero_gate_wavelet = data_file.variables['zero_gate_wavelet'][:].tolist()
    zero_gate_flags = data_file.variables['zero_gate_flag'][:].tolist()
    peak_gate_wavelet = data_file.variables['peak_gate_wavelet'][:].tolist()
    peak_gate_flags = data_file.variables['peak_gate_flag'][:].tolist()
    leading_edge_width = data_file.variables['leading_edge_width'][:].tolist()
    trailing_edge_slope = data_file.variables['trailing_edge_slope'][:].tolist()
    backscatter_wavelet = data_file.variables['backscatter_wavelet'][:].tolist()
    pulse_peakiness = data_file.variables['pulse_peakiness'][:].tolist()
    max_amplitude_wavelet = data_file.variables['max_amplitude_wavelet'][:].tolist()
    corrections_sum_wavelet = data_file.variables['corrections_sum_wavelet'][:].tolist()

    
    data_file.close()     
    
    cryo_points_gimp = []
    cryo_points_arctic = []
    cryo_points_helm = []
    for i,p in enumerate(x_p_gimp):
        cryo_points_gimp.append([x_p_gimp[i], y_p_gimp[i]])
        cryo_points_arctic.append([x_p_arctic[i], y_p_arctic[i]])
        cryo_points_helm.append([x_p_helm[i], y_p_helm[i]])
    
    cryo_points_gimp = np.array(cryo_points_gimp)
    
    "Create CryoSat-2 polygon"
    cryo_poly_gimp = track_poly(cryo_points_gimp, buffer_dist)
    cryo_poly_arctic= track_poly(cryo_points_arctic, buffer_dist)
    cryo_poly_helm= track_poly(cryo_points_helm, buffer_dist)
        
    
    "Run through ICESat2 files, create polygons and check for intersection and find matching points."
    
    cryo_long_gimp = []
    cryo_lat_gimp = []
    cryo_elev_gimp = []
    cryo_date_gimp = []
    ice_elev_gimp = []
    ice_x_gimp = []
    ice_y_gimp = []
    ice_date_gimp = []
    acdc_flag_gimp = []
    le_width_gimp = []
    te_slope_gimp = []
    backscatter_gimp = []
    slope_gimp_matched = []
    
    # cryo_long_arctic = []
    # cryo_lat_arctic = []
    # cryo_elev_arctic = []
    # cryo_date_arctic = []
    # ice_elev_arctic = []
    # ice_x_arctic = []
    # ice_y_arctic = []
    # ice_date_arctic = []
    # acdc_flag_arctic = []
    # le_width_arctic = []
    # te_slope_arctic = []
    # backscatter_arctic = []
    # slope_arctic_matched = []

    
    cryo_long_helm = []
    cryo_lat_helm = []
    cryo_elev_helm = []
    cryo_date_helm = []
    ice_elev_helm = []
    ice_x_helm = []
    ice_y_helm = []
    ice_date_helm = []
    acdc_flag_helm = []
    le_width_helm = []
    te_slope_helm = []
    backscatter_helm = []
    slope_helm_matched = []

    
    
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
        if cryo_poly_gimp.intersects(ice_poly) == False:
            #print('Does not intersect for GIMP.')
            pass
        else: 
            #print (str(indx))
            cryo_elev_matched_gimp, cryo_x_matched_gimp, cryo_y_matched_gimp, cryo_date_matched_gimp, ice_elev_matched_gimp, ice_x_matched_gimp, ice_y_matched_gimp, ice_date_matched_gimp, acdc_flag_matched_gimp, le_width_matched_gimp, te_slope_matched_gimp, backscatter_matched_gimp, slope_matched_gimp = intersect_matches_kdtree(cryo_points_gimp, icesat_points, elevation_gimp, elev_icesat, file_date, idate, buffer_dist, acdc_flag, le_width,te_slope, backscatter, slope_gimp) 
                                    

            cryo_long_gimp.append(cryo_x_matched_gimp)
            cryo_lat_gimp.append(cryo_y_matched_gimp)
            cryo_elev_gimp.append(cryo_elev_matched_gimp)
            cryo_date_gimp.append(cryo_date_matched_gimp)
            ice_elev_gimp.append(ice_elev_matched_gimp)
            ice_x_gimp.append(ice_x_matched_gimp)
            ice_y_gimp.append(ice_y_matched_gimp)
            ice_date_gimp.append(ice_date_matched_gimp)
            acdc_flag_gimp.append(acdc_flag_matched_gimp)
            le_width_gimp.append(le_width_matched_gimp)
            te_slope_gimp.append(te_slope_matched_gimp)
            backscatter_gimp.append(backscatter_matched_gimp)
            slope_gimp_matched.append(slope_matched_gimp)
            
            
        # if cryo_poly_arctic.intersects(ice_poly) == False:
        #     #print('Does not intersect for Arctic.')
        #     pass
        # else:
            
        #     cryo_elev_matched_arctic, cryo_x_matched_arctic, cryo_y_matched_arctic, cryo_date_matched_arctic, ice_elev_matched_arctic, ice_x_matched_arctic, ice_y_matched_arctic, ice_date_matched_arctic, acdc_flag_matched_arctic, le_width_matched_arctic,te_slope_matched_arctic, backscatter_matched_arctic, slope_matched_arctic = intersect_matches_kdtree(cryo_points_arctic, icesat_points, elevation_arctic, elev_icesat, file_date, idate, buffer_dist, acdc_flag, le_width, te_slope, backscatter, slope_arctic) 
            
        #     cryo_long_arctic.append(cryo_x_matched_arctic)
        #     cryo_lat_arctic.append(cryo_y_matched_arctic)
        #     cryo_elev_arctic.append(cryo_elev_matched_arctic)
        #     cryo_date_arctic.append(cryo_date_matched_arctic)
        #     ice_elev_arctic.append(ice_elev_matched_arctic)
        #     ice_x_arctic.append(ice_x_matched_arctic)
        #     ice_y_arctic.append(ice_y_matched_arctic)
        #     ice_date_arctic.append(ice_date_matched_arctic)
        #     acdc_flag_arctic.append(acdc_flag_matched_arctic)    
        #     le_width_arctic.append(le_width_matched_arctic)
        #     te_slope_arctic.append(te_slope_matched_arctic)
        #     backscatter_arctic.append(backscatter_matched_arctic)
        #     slope_arctic_matched.append(slope_matched_arctic)

        
        
        if cryo_poly_helm.intersects(ice_poly) == False:
            #print('Does not intersect for Arctic.')
            pass
        else:
            
            cryo_elev_matched_helm, cryo_x_matched_helm, cryo_y_matched_helm, cryo_date_matched_helm, ice_elev_matched_helm, ice_x_matched_helm, ice_y_matched_helm, ice_date_matched_helm, acdc_flag_matched_helm, le_width_matched_helm,te_slope_matched_helm, backscatter_matched_helm, slope_matched_helm = intersect_matches_kdtree(cryo_points_helm, icesat_points, elevation_helm, elev_icesat, file_date, idate, buffer_dist, acdc_flag, le_width,te_slope, backscatter, slope_helm) 
            
            cryo_long_helm.append(cryo_x_matched_helm)
            cryo_lat_helm.append(cryo_y_matched_helm)
            cryo_elev_helm.append(cryo_elev_matched_helm)
            cryo_date_helm.append(cryo_date_matched_helm)
            ice_elev_helm.append(ice_elev_matched_helm)
            ice_x_helm.append(ice_x_matched_helm)
            ice_y_helm.append(ice_y_matched_helm)
            ice_date_helm.append(ice_date_matched_helm)
            acdc_flag_helm.append(acdc_flag_matched_helm)
            le_width_helm.append(le_width_matched_helm)
            te_slope_helm.append(te_slope_matched_helm)
            backscatter_helm.append(backscatter_matched_helm)
            slope_helm_matched.append(slope_matched_helm)

        
        
    "Prepare output for saving"   
    long_cryo_gimp = [np.float32(x) for sub in cryo_long_gimp for x in sub]
    lat_cryo_gimp = [np.float32(x) for sub in cryo_lat_gimp for x in sub]
    elev_cryo_gimp = [np.float32(x) for sub in cryo_elev_gimp for x in sub]
    date_cryo_gimp = [np.float32(x) for sub in cryo_date_gimp for x in sub]
    long_ice_gimp = [np.float32(x) for sub in ice_x_gimp for x in sub]
    lat_ice_gimp = [np.float32(x) for sub in ice_y_gimp for x in sub]
    elev_ice_gimp = [np.float32(x) for sub in ice_elev_gimp for x in sub]
    date_ice_gimp = [np.float32(x) for sub in ice_date_gimp for x in sub]
    flag_acdc_gimp = [np.float32(x) for sub in acdc_flag_gimp for x in sub]
    le_width_gimp = [np.float32(x) for sub in le_width_gimp for x in sub]
    te_slope_gimp = [np.float32(x) for sub in te_slope_gimp for x in sub]
    backscatter_gimp = [np.float32(x) for sub in backscatter_gimp for x in sub]
    slope_gimp_matched = [np.float32(x) for sub in slope_gimp_matched for x in sub]
    
    # long_cryo_arctic = [np.float32(x) for sub in cryo_long_arctic for x in sub]
    # lat_cryo_arctic = [np.float32(x) for sub in cryo_lat_arctic for x in sub]
    # elev_cryo_arctic = [np.float32(x) for sub in cryo_elev_arctic for x in sub]
    # date_cryo_arctic = [np.float32(x) for sub in cryo_date_arctic for x in sub]
    # long_ice_arctic = [np.float32(x) for sub in ice_x_arctic for x in sub]
    # lat_ice_arctic = [np.float32(x) for sub in ice_y_arctic for x in sub]
    # elev_ice_arctic = [np.float32(x) for sub in ice_elev_arctic for x in sub]
    # date_ice_arctic = [np.float32(x) for sub in ice_date_arctic for x in sub]
    # flag_acdc_arctic = [np.float32(x) for sub in acdc_flag_arctic for x in sub]
    # le_width_arctic = [np.float32(x) for sub in le_width_arctic for x in sub]
    # te_slope_arctic = [np.float32(x) for sub in te_slope_arctic for x in sub]
    # backscatter_arctic = [np.float32(x) for sub in backscatter_arctic for x in sub]
    # slope_arctic_matched = [np.float32(x) for sub in slope_arctic_matched for x in sub]
    
    long_cryo_helm = [np.float32(x) for sub in cryo_long_helm for x in sub]
    lat_cryo_helm = [np.float32(x) for sub in cryo_lat_helm for x in sub]
    elev_cryo_helm = [np.float32(x) for sub in cryo_elev_helm for x in sub]
    date_cryo_helm = [np.float32(x) for sub in cryo_date_helm for x in sub]
    long_ice_helm = [np.float32(x) for sub in ice_x_helm for x in sub]
    lat_ice_helm = [np.float32(x) for sub in ice_y_helm for x in sub]
    elev_ice_helm = [np.float32(x) for sub in ice_elev_helm for x in sub]
    date_ice_helm = [np.float32(x) for sub in ice_date_helm for x in sub]
    flag_acdc_helm = [np.float32(x) for sub in acdc_flag_helm for x in sub]
    le_width_helm = [np.float32(x) for sub in le_width_helm for x in sub]
    te_slope_helm = [np.float32(x) for sub in te_slope_helm for x in sub]
    backscatter_helm = [np.float32(x) for sub in backscatter_helm for x in sub]
    slope_helm_matched = [np.float32(x) for sub in slope_helm_matched for x in sub]

    dcg = [i for i in date_cryo_gimp]
    dig = [i for i in date_ice_gimp]
    # dca = [i for i in date_cryo_arctic]
    # dia = [i for i in date_ice_arctic]
    dch = [i for i in date_cryo_helm]
    dih = [i for i in date_ice_helm]

    g_length_output = len(elev_cryo_gimp)
#    a_length_output = len(elev_cryo_arctic)
    h_length_output = len(elev_cryo_helm)
    
    while g_length_output!=0 or h_length_output!=0: 
        print('All the matches found! Saving file ...')    
        "Save this matched points in an output file"
        #save the data into a .nc format:
        #print('Saving the data to a .nc file...')        
        out_name = str(file_out_name) + '_Month_' + str(month) + '_' + str(year) + '_' + str(bucket_size) + 'day_bucket.nc'
        out_data_file = os.path.join(outdirectory,out_name)
        out_file = Dataset(out_data_file, 'w', format='NETCDF4_CLASSIC')
        
        #creating dimentions for the data:
        gpoints = out_file.createDimension('gpoints', g_length_output)
#        apoints = out_file.createDimension('apoints', a_length_output)
        hpoints = out_file.createDimension('hpoints', h_length_output)
        
        #assigning variables for the data with above dimentions:
        date_cryo_g = out_file.createVariable('date_cryo_g', np.int32, ('gpoints',))
        x_cryo_g = out_file.createVariable('px_cryo_g', np.float32, ('gpoints',))
        y_cryo_g = out_file.createVariable('py_cryo_g', np.float32, ('gpoints',))
        elev_cryo_g = out_file.createVariable('elevation_cryo_g', np.float32, ('gpoints',))
        date_ice_g = out_file.createVariable('date_ice_g', np.int32, ('gpoints',))
        x_ice_g = out_file.createVariable('px_ice_g', np.float32, ('gpoints',))
        y_ice_g = out_file.createVariable('py_ice_g', np.float32, ('gpoints',))
        elev_ice_g = out_file.createVariable('elevation_ice_g', np.float32, ('gpoints',))
        slope_g = out_file.createVariable('slope_g', np.float32, ('gpoints',))
        
        # date_cryo_a = out_file.createVariable('date_cryo_a', np.int32, ('apoints',))
        # x_cryo_a = out_file.createVariable('px_cryo_a', np.float32, ('apoints',))
        # y_cryo_a = out_file.createVariable('py_cryo_a', np.float32, ('apoints',))
        # elev_cryo_a = out_file.createVariable('elevation_cryo_a', np.float32, ('apoints',))
        # date_ice_a = out_file.createVariable('date_ice_a', np.int32, ('apoints',))
        # x_ice_a = out_file.createVariable('px_ice_a', np.float32, ('apoints',))
        # y_ice_a = out_file.createVariable('py_ice_a', np.float32, ('apoints',))
        # elev_ice_a = out_file.createVariable('elevation_ice_a', np.float32, ('apoints',))
        # slope_a = out_file.createVariable('slope_a', np.float32, ('apoints',))
        
        # ad_flag = out_file.createVariable('acdc_flag', np.float32, ('gpoints',))
        
        date_cryo_h = out_file.createVariable('date_cryo_h', np.int32, ('hpoints',))
        x_cryo_h = out_file.createVariable('px_cryo_h', np.float32, ('hpoints',))
        y_cryo_h = out_file.createVariable('py_cryo_h', np.float32, ('hpoints',))
        elev_cryo_h = out_file.createVariable('elevation_cryo_h', np.float32, ('hpoints',))
        date_ice_h = out_file.createVariable('date_ice_h', np.int32, ('hpoints',))
        x_ice_h = out_file.createVariable('px_ice_h', np.float32, ('hpoints',))
        y_ice_h = out_file.createVariable('py_ice_h', np.float32, ('hpoints',))
        elev_ice_h = out_file.createVariable('elevation_ice_h', np.float32, ('hpoints',))
        slope_h = out_file.createVariable('slope_h', np.float32, ('hpoints',))

        leading_edge_width = out_file.createVariable('leading_edge_width', np.float32, ('gpoints',))
        trailing_edge_slope = out_file.createVariable('trailing_edge_slope', np.float32, ('gpoints',))
        backscatter_wavelet = out_file.createVariable('backscatter_wavelet', np.float32, ('gpoints',))
        
        #assign atributes to the file (global) and to the variables (variable):
        out_file.description = 'This file contains matched CryoSat-2 data points and IceSat-2 data points. Coordinates have been convarted to stereographic polar projection from the decimal degree latitude and longitude.'
        out_file.source = 'Original data from an IceSat-2 (NASA). CryoSat-2 data computed witha retracker by KM Sejan 2019.'
        
        date_cryo_g.units = 'yyyymmdd'
        x_cryo_g.units = 'meters'
        y_cryo_g.units = 'meters'
        elev_cryo_g.units = 'meters'
        date_ice_g.units = 'yyyymmdd'
        x_ice_g.units = 'meters'
        y_ice_g.units = 'meters'
        elev_ice_g.units = 'meters'
        slope_g.units = 'radians'
        # date_cryo_a.units = 'yyyymmdd'
        # x_cryo_a.units = 'meters'
        # y_cryo_a.units = 'meters'
        # elev_cryo_a.units = 'meters'
        # date_ice_a.units = 'yyyymmdd'
        # x_ice_a.units = 'meters'
        # y_ice_a.units = 'meters'
        # elev_ice_a.units = 'meters'
        # slope_a.units = 'radians'
        date_cryo_h.units = 'yyyymmdd'
        x_cryo_h.units = 'meters'
        y_cryo_h.units = 'meters'
        elev_cryo_h.units = 'meters'
        date_ice_h.units = 'yyyymmdd'
        x_ice_h.units = 'meters'
        y_ice_h.units = 'meters'
        elev_ice_h.units = 'meters'
        slope_h.units = 'radians'
        leading_edge_width.units = 'meters'
        trailing_edge_slope.units = 'unitless, but theoretically it represents the change of return power over one bin.'
        backscatter_wavelet.units = 'count'
        
        date_cryo_g[:] = dcg
        x_cryo_g[:] = long_cryo_gimp
        y_cryo_g[:] = lat_cryo_gimp
        elev_cryo_g[:] = elev_cryo_gimp
        date_ice_g[:] = dig
        x_ice_g[:] = long_ice_gimp
        y_ice_g[:] = lat_ice_gimp
        elev_ice_g[:] = elev_ice_gimp
        slope_g[:] = slope_gimp_matched
        
        # date_cryo_a[:] = dca
        # x_cryo_a[:] = long_cryo_arctic
        # y_cryo_a[:] = lat_cryo_arctic
        # elev_cryo_a[:] = elev_cryo_arctic
        # date_ice_a[:] = dia
        # x_ice_a[:] = long_ice_arctic
        # y_ice_a[:] = lat_ice_arctic
        # elev_ice_a[:] = elev_ice_arctic
        # slope_a[:] = slope_arctic_matched
        
        # ad_flag[:] = flag_acdc_gimp
        
        date_cryo_h[:] = dch
        x_cryo_h[:] = long_cryo_helm
        y_cryo_h[:] = lat_cryo_helm
        elev_cryo_h[:] = elev_cryo_helm
        date_ice_h[:] = dih
        elev_ice_h[:] = elev_ice_helm
        x_ice_h[:] = long_ice_helm
        y_ice_h[:] = lat_ice_helm
        slope_h[:] = slope_helm_matched
        
        leading_edge_width[:] = le_width_gimp
        trailing_edge_slope[:] = te_slope_gimp
        backscatter_wavelet[:] = backscatter_gimp

        
        out_file.close()
        
        print('...SAVED.')
        break
    print('Month ' + str(month) + ' ' + str(year) + ' has ' + str(g_length_output) + ' matches for GIMP and ' + str(h_length_output) + ' matches for Helm DEM.')


