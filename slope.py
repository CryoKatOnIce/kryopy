# -*- coding: utf-8 -*-
"""
Slope correction for the CryoSat2 Baseline D data.

Created by: Kat Sejan 21st February 2020.
Last edited by: Kat Sejan 8th October 2021.

"""

def coordinates_geo_to_polar(nadir_lat, nadir_long):
    
    from pyproj import Transformer
    
    inProj = 'epsg:4326'
    outProj = 'epsg:3413'
    polar_transformer = Transformer.from_crs(inProj, outProj)
    
    x = []
    y = []

    for i3,l in enumerate(nadir_long):
        x_trans, y_trans  = polar_transformer.transform(nadir_lat[i3], nadir_long[i3])
        x.append(x_trans)
        y.append(y_trans)

    return x, y



def dem_open(dem = 'gimp'):

    import os
    from pathlib import Path
    from netCDF4 import Dataset


    """Read in the tiles bounds"""
    if dem == 'gimp': 
        dem_path = Path('/Users/kat/DATA/DEM_GIMP_processed/')
        dem_file = 'GIMP_dem_mosaic_1km_sae.nc'

    elif dem == 'arctic':
        dem_path = Path('/Users/kat/DATA/ArcticDEM_processed/')
        dem_file = 'ArcticDEM_mosaic_1km_sae.nc'

    elif dem == 'helm':
        dem_path = Path('/Users/kat/DATA/Helm_DEM_processed/')
        dem_file = 'DEM_GRE_CS_20130826_sae.nc'

    else:
        print('ERROR: Input DEM not specified!')
    
    """Get DEM raster values."""
    path_to_file = os.path.join(dem_path, dem_file)
    
    dem_file = Dataset(path_to_file, "r", format="NETCDF4")
    dem_x = dem_file.variables['polar_x'][:]
    dem_y = dem_file.variables['polar_y'][:]
    
    elev_rast = dem_file.variables['elevation'][:]
    slope_rast = dem_file.variables['slope'][:]
    aspect_rast = dem_file.variables['aspect'][:]

    dem_file.close()
    
    return dem_x, dem_y, elev_rast, slope_rast, aspect_rast


def dem_values(x, y, dem_x, dem_y, elev_rast, slope_rast, aspect_rast):
    
    import numpy as np
    import math


    """Find nearest values for each point."""
    nearest_elevation = []
    nearest_x = []
    nearest_y = []
    nearest_aspect = []
    nearest_slope = []    
    for i,p in enumerate(x):           
        #find closest x and y coordinates in the array:
        idx_x = (np.abs(dem_x - x[i])).argmin()
        idx_y = (np.abs(dem_y - y[i])).argmin()
        
        elev_val = elev_rast[idx_y, idx_x]
        slp_val = slope_rast[idx_y, idx_x] #not abs() if using ascending -pi and descending as is
        asp_val = aspect_rast[idx_y, idx_x]
    
        nearest_x.append(dem_x[idx_x])
        nearest_y.append(dem_y[idx_y])
        nearest_elevation.append(elev_val)        
        nearest_slope.append(math.radians(slp_val))
        nearest_aspect.append(math.radians(asp_val))
    

    return nearest_x, nearest_y, nearest_elevation, nearest_slope, nearest_aspect


def hurkmans_correction(x, y, nadir_range, sat_alt, nearest_slope, nearest_aspect):
    #Hurkmans method needs anticlockwise aspect, therfore previously generated aspect is anticlockwise
    
    import numpy as np
    
    hurkmans_displacement = []
    hurkmans_dx = []
    hurkmans_dy = []
    hurkmans_x = []
    hurkmans_y = []
    hurkmans_range = []
    
    for i,r in enumerate(nadir_range):        
        
        if r == -9999:
            hurkmans_displacement.append(np.float32(-9999))
            hurkmans_dx.append(np.float32(-9999))
            hurkmans_dy.append(np.float32(-9999))
            hurkmans_x.append(np.float32(-9999999999))
            hurkmans_y.append(np.float32(-9999999999)) 
            hurkmans_range.append(np.float32(-9999))
        else:        
            #calculated the displacement and diffrence in x and y:
            displacement = sat_alt[i] * np.sin(nearest_slope[i]) * np.cos(nearest_slope[i]) 
            dx = displacement * np.sin(nearest_aspect[i] - np.pi)
            dy = displacement * np.cos(nearest_aspect[i] - np.pi)
                        
            #calculating the new coordinates
            x_corr_smooth = x[i] + dx # these should be added as the dx and dy already have a positive or negative sign based on the displacement direction opposite to the aspect
            y_corr_smooth = y[i] + dy
            
            #Range correction
            range_corr = nadir_range[i] * np.cos(nearest_slope[i]) #in [m] #second method by Hurkmans, not the direct method of Nilsson
            
            #Saving calculculated values
            hurkmans_displacement.append(np.float32(displacement))
            hurkmans_dx.append(np.float32(dx))
            hurkmans_dy.append(np.float32(dy))
            hurkmans_x.append(np.float32(x_corr_smooth))
            hurkmans_y.append(np.float32(y_corr_smooth)) 
            hurkmans_range.append(np.float32(range_corr))

    return hurkmans_x, hurkmans_y, hurkmans_dx, hurkmans_dy, hurkmans_displacement, hurkmans_range

def earth_curvature_correction(nadir_long, nadir_lat, nadir_range, sat_alt, nearest_x, nearest_y, nearest_slope, nearest_aspect, dem = 'gimp'):
    #specify the DEM to use: 'gimp' or 'arctic'
    #the smoothing is a radius over which the slope is smoothed in km, value of 0.03 is equal to no smoothing.
    
    import numpy as np
    import math
    from pyproj import Transformer
    import pyproj
    import pymap3d as pm

        
    """Apply the slope correction and curvature correction with the found slope and aspect values."""    
    #specifying constants:
    """
    a = np.float(6378137.0) #Earth's equatorial radius [m]
    b = np.float(6356752.3) #Earth's polar radius [m]
    e = np.sqrt(1-(b**2/a**2)) #[m] Bamber has:(a**2-b**2)/a**2
    """
    #as in matlab functions for wllipsoid WGS84
    a=6378137.0;
    finv=298.257223563;
    f=1/finv
    b=a*(1-f)
    e=1-(1-f)**2

    
    "Earth curvature correction"
    #comupte the ascending/decending flag:
    acdc = ['D' if i==0 else 'D' if nadir_lat[i-1]-nadir_lat[i]>0 else 'A' for i,l in enumerate(nadir_lat)]
    acdc_flag = all(flag=='D' for flag in acdc) #True if decending, false if adcending 
    #ad_flag = [1 if i==0 else 1 if nadir_lat[i-1]-nadir_lat[i]>0 else 0 for i,l in enumerate(nadir_lat)]
    ad_flag = [1 if acdc_flag == True else 0 for i in acdc]
    
    inProj = 'epsg:4326'
    outProj = 'epsg:3413'
    polar_transformer = Transformer.from_crs(inProj, outProj)
    
    inProj = 'epsg:3413'
    outProj = 'epsg:4326'
    geo_transformer = Transformer.from_crs(inProj, outProj)

    geodesic = pyproj.Geod(ellps='WGS84')
    
    ecc_elevation = []
    ecc_correction = []
    ecc_range = []
    ecc_x = []
    ecc_y = []
    for i,r in enumerate(nadir_range):
        
        if r == -9999:
            ecc_correction.append(np.float32(-9999))
            ecc_range.append(np.float32(-9999))
            ecc_elevation.append(np.float32(-9999))
            ecc_x.append(np.float32(-9999999999)) #[m] 
            ecc_y.append(np.float32(-9999999999)) #[m]

        else:
            """
            #compute N and M
            lat_rad = math.radians(nadir_lat[i])
            lon_rad = math.radians(nadir_long[i])
            alpha = nearest_aspect[i] - np.pi #-nearest_aspect[i] + (np.pi/2) - np.pi #correcting for East = 0 with +pi/2 and aniclockwise to clockwise with -angle
                        
            r_m = a*(1-e**2)/((1-e**2*(np.sin(lat_rad)**2))**(3/2)) #[m] 
            r_n = (a*np.cos(lat_rad))/np.sqrt((1-e**2*(np.sin(lat_rad)**2))) #[m] Bamber does not include *np.cos(lat_rad) in the numerator
            
            #compute radius of curvature
            r_alpha = (r_m * r_n) / ((r_m*np.cos(lat_rad)*(np.sin(alpha)**2)) + (r_n*(np.cos(alpha)**2))) #Bamber does not include *np.cos(lat_rad) in the denominator

            #compute correction
            curv_corr = ((nadir_range[i]*np.sin(nearest_slope[i]))**2)/(2.0*r_alpha) #[m] 
                    
            #Corrected Elevation
            range_corr_curv = nadir_range[i]*np.cos(nearest_slope[i]) + curv_corr #[m]
            elevation = sat_alt[i] - range_corr_curv #[m]
            lon_smooth = math.degrees(lon_rad +((nadir_range[i]*np.sin(nearest_slope[i])*np.sin(alpha))/r_n)) #longitude [dd]
            lat_smooth = math.degrees(lat_rad +((nadir_range[i]*np.sin(nearest_slope[i])*np.cos(alpha))/r_m)) #latitude [dd]
            
            #change degrees to polar coordinates for the final corrected coordinates
            x_smooth, y_smooth  = polar_transformer.transform(lat_smooth,lon_smooth)
            
            current_method = dict( {'long' : nadir_long[i], 'lat': nadir_lat[i], 'long_corr' :lon_smooth, 'lat_corr' :lat_smooth, 
                                    'elevation' : elevation, 'alpha': alpha, 'r_alpha': r_alpha, 'Curv_correction': curv_corr})
                                    
            
            #Saving calculculated values
            ecc_correction.append(np.float32(curv_corr))
            ecc_range.append(np.float32(range_corr_curv))
            ecc_elevation.append(np.float32(elevation))
            ecc_x.append(np.float32(x_smooth)) #[m] 
            ecc_y.append(np.float32(y_smooth)) #[m]
            
            """
            aspect = nearest_aspect[i] - np.pi
            
            DEMx = nearest_x[i] + np.sin(aspect)
            DEMy = nearest_y[i] + np.cos(aspect)
            
            lon1, lat1 = geo_transformer.transform(nearest_x[i], nearest_y[i])
            lon2, lat2 = geo_transformer.transform(DEMx, DEMy)
            
            fwd_azimuth,back_azimuth,distance = geodesic.inv(lon1, lat1, lon2, lat2)
            
            azimuth = fwd_azimuth 
            
            i_azimuth = np.cos(azimuth)+1j*np.sin(azimuth)  #(np.cos(azimuth)+1j*np.sin(azimuth)).real
            
            alpha = math.degrees(np.angle(i_azimuth)) % 360 #aA_s
            #wrapTo360 >>> % 360
            #rad2deg >>> math.degrees
            #angle >>> np.angle
            #interp2 >>> it is just finding the nearest point to our xy in a grid of azimuth values (nearest_azimuth)
            #we have however already used 
            
            
            #sA_s >>> nearest_slope[i]
            
            lat_rad = math.radians(nadir_lat[i])
            lon_rad = math.radians(nadir_long[i])
                        
            r_m = a*(1-e**2)/((1-e**2*(np.sin(lat_rad)**2))**(3/2)) #[m] 
            r_n = a/np.sqrt((1-e**2*(np.sin(lat_rad)**2))) #[m] Bamber does not include *np.cos(lat_rad) in the numerator
            
            #compute radius of curvature
            r_alpha = (r_m * r_n) / ((r_m*(np.sin(alpha)**2)) + (r_n*(np.cos(alpha)**2))) #Bamber does not include *np.cos(lat_rad) in the denominator, but this foes not affect the resulting r_alpha
            
            
            R_s = r_alpha + sat_alt[i] #??? is CS.GEO.H >>> sat_alt[i] ?
            
            curv_corr = math.degrees(math.asin(nadir_range[i]*np.sin(nearest_slope[i]) / R_s))
            
            elevation = (R_s * np.sin(nearest_slope[i] - curv_corr)/np.sin(nearest_slope[i])) - r_alpha
            
            #matlab method to get corrected lona nd lat:
            arclength = r_alpha*curv_corr
            lonc, latc, back_azimuth2 = geodesic.fwd(nadir_long[i], nadir_lat[i], alpha, arclength)
            
            dart_method = dict( {'long' : nadir_long[i], 'lat': nadir_lat[i], 'long_corr' :lonc, 'lat_corr' :latc, 
                                    'elevation' : elevation, 'alpha': alpha, 'r_alpha': r_alpha, 'Curv_correction': curv_corr})

            
            #bamber method in matlab for corrected lona and lat:
            dx = r_alpha * curv_corr * np.cos(alpha)
            dy = r_alpha * curv_corr * np.sin(alpha)
            
            #convert nadir coordinates to cartesian:
            v=a/np.sqrt(1 - e * np.sin(nadir_lat[i]) * np.sin(nadir_lat[i]))
            nadir_x = (v + nadir_range[i]) * np.cos(nadir_lat[i]) * np.cos(nadir_long[i])
            nadir_y = (v + nadir_range[i]) * np.cos(nadir_lat[i]) * np.sin(nadir_long[i])
            nadir_z = (v * (1 - e) + nadir_range[i]) * np.sin(nadir_lat[i])
            
            x_c = nadir_x + dx
            y_c = nadir_y + dy
            
            #convert coordinates from cartesian to geo:
            # Latitude and height convergence criteria
            elat = 1e-12
            eht = 1e-5
            # Initial values for iteration
            p = np.sqrt(x_c * x_c + y_c * y_c)
            lat = math.atan2(nadir_z, p*(1-e))
            h=0
            dh=1
            dlat=1
            # Iterate until lat & h converge to elat & eht
            while dlat>elat or dh>eht:
                lat0 = lat
                h0 = h
                v = a/np.sqrt(1 - e * np.sin(lat) * np.sin(lat))
                h = p/np.cos(lat) - v
                lat = math.atan2(nadir_z, p * (1-e * v/(v+nadir_z)))
                dlat = abs(lat - lat0)
                dh = abs(h - h0)
            
            lon_c = math.degrees(math.atan2(y_c, x_c))
            lat_c = math.degrees(dlat)
            range_corr_curv = dh
            
            #change degrees to polar coordinates for the final corrected coordinates
            x_smooth, y_smooth  = polar_transformer.transform(lat_c,lon_c)
            
            bamber_method = dict( {'long' : nadir_long[i], 'lat': nadir_lat[i], 'long_corr' :lon_c, 'lat_corr' :lat_c, 
                                    'elevation' : elevation, 'alpha': alpha, 'r_alpha': r_alpha, 'Curv_correction': curv_corr})

            
            #Saving calculculated values
            ecc_correction.append(np.float32(curv_corr))
            ecc_range.append(np.float32(range_corr_curv))
            ecc_elevation.append(np.float32(elevation))
            ecc_x.append(np.float32(x_smooth)) #[m] 
            ecc_y.append(np.float32(y_smooth)) #[m]
            """
            
    
    return ecc_x, ecc_y, ecc_range, ecc_elevation, ecc_correction, ad_flag


def slope_correction(infile, dem_x, dem_y, elev_rast, slope_rast, aspect_rast, alpha = [0.3], dem = 'gimp'):
    #writes output to the same file as already retracked data
    #specify the DEM to use: 'gimp' or 'arctic'
    
    from netCDF4 import Dataset
    import numpy as np     
    from shapely.geometry import Point
    from shapely.geometry import box
       
    
    try:
        "Reading in the data from file:"
        data_file = Dataset(infile, mode="r+")
        #print('File read in, file name: ' + str(infile))
        
        nadir_long = data_file.variables['lon_20_ku'][:]
        nadir_lat = data_file.variables['lat_20_ku'][:]
        sat_alt = data_file.variables['alt_20_ku'][:]
        
        
        #Converting the wavelet location to polar stereographic projection to match the DEM extend
        x, y = coordinates_geo_to_polar(nadir_lat, nadir_long)        
                
        nearest_x, nearest_y, nearest_elevation, nearest_slope, nearest_aspect = dem_values(x,y, dem_x, dem_y, elev_rast, slope_rast, aspect_rast)

        
        if dem == 'gimp': 
            extension = '_gimp'
            description = ' using a GIMP DEM (Version 1, collected between 30 June 1999 to 04 September 2002).'
            source = 'KM Sejan code for slope correction in retracker. Calculations based on GIMP DEM (Version 1, collected between 30 June 1999 to 04 September 2002).'
        elif dem == 'arctic':
            extension = '_arctic'
            description = ' using ArcticDEM (Release 7, data collected from 2015 to 2017). Replacing the non value value -9999 with an average of serrounding cells.'
            source = 'KM Sejan code for slope correction in retracker. Calculations based on ArcticDEM (Release 7, data collected from 2015 to 2017).'
        elif dem == 'helm':
            extension = '_helm'
            description = ' using Helm (2014). Replacing the non value value -9999 with an average of serrounding cells.'
            source = 'KM Sejan code for slope correction in retracker. Calculations based on Helm (2014).'


        dem_elev_n_s = data_file.createVariable('nearest_elevation' + str(extension), np.float32, ('time_20_ku',))
        aspect_n_s = data_file.createVariable('nearest_aspect_radian' + str(extension), np.float32, ('time_20_ku',))
        slope_n_s = data_file.createVariable('nearest_slope_radian' + str(extension), np.float32, ('time_20_ku',))

        dem_elev_n_s.description = 'Elevation of the nearest point in the DEM' + str(description) + ' Matched to the measurement (wavelet) point (both points in polar coordinate system).'
        dem_elev_n_s.units = 'meters'
        dem_elev_n_s.source = str(source)
        aspect_n_s.description = 'Slope aspect of the nearest point in the DEM' + str(description) + ' Matched to the measurement (wavelet) point (both points in polar coordinate system).'
        aspect_n_s.units = 'radians'
        aspect_n_s.source = str(source)
        slope_n_s.description = 'Slope of the nearest point in the DEM' + str(description) + ' Matched to the measurement (wavelet) point (both points in polar coordinate system).'
        slope_n_s.units = 'radians'
        slope_n_s.source = str(source)
        
        dem_elev_n_s[:] = nearest_elevation
        aspect_n_s[:] = nearest_aspect
        slope_n_s[:] = nearest_slope


        for ai,aa in enumerate(alpha):
            
            threshold = str(aa).replace('.', '_')

            nadir_range = data_file.variables['range_wavelet_'+threshold][:]
            
                    
            hurkmans_x, hurkmans_y, hurkmans_dx, hurkmans_dy, hurkmans_displacement, hurkmans_range = hurkmans_correction(x, y, nadir_range, sat_alt, nearest_slope, nearest_aspect)

            ecc_x, ecc_y, ecc_range, ecc_elevation, ecc_correction, ad_flag = earth_curvature_correction(nadir_long, nadir_lat, nadir_range, sat_alt, nearest_x, nearest_y, nearest_slope, nearest_aspect, dem = 'gimp')

            displ_s = data_file.createVariable('hurkmans_displacement_' + str(threshold) + str(extension), np.float32, ('time_20_ku',))
            dx_calc_s = data_file.createVariable('hurkmans_dx_' + str(threshold) + str(extension), np.float32, ('time_20_ku',))
            dy_calc_s = data_file.createVariable('hurkmans_dy_' + str(threshold) + str(extension), np.float32, ('time_20_ku',))
            calc_range = data_file.createVariable('hurkmans_range_' + str(threshold) + str(extension), np.float32, ('time_20_ku',))            
            x_s = data_file.createVariable('ecc_x_' + str(threshold) + str(extension), np.float32, ('time_20_ku',))
            y_s = data_file.createVariable('ecc_y_' + str(threshold) + str(extension), np.float32, ('time_20_ku',))
            c_elev = data_file.createVariable('ecc_elevation_' + str(threshold) + str(extension), np.float32, ('time_20_ku',))
            calc_range_curv = data_file.createVariable('ecc_range_' + str(threshold) + str(extension), np.float32, ('time_20_ku',))
            curv_corr = data_file.createVariable('ecc_correction_' + str(threshold) + str(extension), np.float32, ('time_20_ku',))
            acdc_flag = data_file.createVariable('acdc_flag_' + str(threshold) + str(extension), np.int32, ('time_20_ku',))

            displ_s.description = 'Horizontal displacement calculated based on slope at a nearest point to the measurement (wavelet) point, ' + str(description)
            displ_s.units = 'meters'
            displ_s.source = str(source)
            dx_calc_s.description = 'Horizontal displacement in the direction of x coordinate, calculated based on displacement (slope) and aspect at a nearest point to the measurement (wavelet) point,' + str(description)
            dx_calc_s.units = 'meters'
            dx_calc_s.source = str(source)
            dy_calc_s.description = 'Horizontal displacement in the direction of y coordinate, calculated based on displacement (slope) and aspect at a nearest point to the measurement (wavelet) point,' + str(description)
            dy_calc_s.units = 'meters'
            dy_calc_s.source = str(source)
            calc_range.description = 'Range corrected for slope error'
            calc_range.units = 'meters' 
            calc_range.sourse = str(source)
            x_s.description = 'Polar x coordinate of the measurement (wavelet) corrected for slope error' + str(description)
            x_s.units = 'meters'
            x_s.source = str(source) 
            y_s.description = 'Polar y coordinate of the measurement (wavelet) corrected for slope error' + str(description)
            y_s.units = 'meters'
            y_s.source = str(source)
            c_elev.description = 'Elevation calculated after appluing slope and curvature correction to the range'
            c_elev.units = 'meters'
            c_elev.source = str(source)
            calc_range_curv.description = 'Range corrected for slope error and for Earths curvature'
            calc_range_curv.units = 'meters' 
            calc_range_curv.sourse = str(source)
            curv_corr.description = 'Earths curvature correction applied to range' + str(description)
            curv_corr.units = 'meters'
            curv_corr.source = str(source)
            acdc_flag.description = 'descending tracks are flagged 1, ascending tracks are flagged 0.'
            acdc_flag.units = '0 or 1'
            acdc_flag.source = str(source)

            calc_range[:] = hurkmans_range
            displ_s[:] = hurkmans_displacement
            dx_calc_s[:] = hurkmans_dx
            dy_calc_s[:] =  hurkmans_dy
            x_s[:] = ecc_x
            y_s[:] = ecc_y
            c_elev[:] = ecc_elevation
            calc_range_curv[:] = ecc_range
            curv_corr[:] = ecc_correction
            acdc_flag[:] = ad_flag


        data_file.close()
        
    
    except (IOError, RuntimeError, ValueError, IndexError, UnboundLocalError) as e:
        print(e)
        
