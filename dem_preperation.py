# -*- coding: utf-8 -*-
"""
For Python 2.7 with richdem and gdal.
This is a code for computing slope, aspect and curvature from GIMP DEM with 30m resolution. Also, for cimputing point coordniates (decimal degrees).

Created by: Kat Sejan 23rd September 2020. Last edidted on 11th October 2021
"""


def dem_tif_to_netcdf(directory, name, slope = True, aspect = True):
    #directory where dem or dem tiles are; string path
    #the name of the dem geotif or the common name of the tiles; string of name
    
    import gdal
    import gdalconst
    import os
    import numpy as np
    from pyproj import Proj, transform
    from netCDF4 import Dataset
    import time
    import csv

    dem_tiles = [file_name for file_name in os.listdir(directory) if file_name.endswith(name + '.tif')]

    tiles = []
    minx_tiles = []
    miny_tiles = []
    maxx_tiles = []
    maxy_tiles = []
    
    for indx,d in enumerate(dem_tiles):
    
        tiles.append(d)
                
        #path to the DEM:
        elevation_file = str(d)
        elevation_dem_path = os.path.join(directory, elevation_file)

        #read the DEM:
        elev_dem = gdal.Open(elevation_dem_path)

        #extract the elevation values:
        elev_band = elev_dem.GetRasterBand(1)
        dem_elev = elev_band.ReadAsArray() #this array is 15000 in length => the y is the length of the array

        #create tiles info for csv file
        width2 = elev_dem.RasterXSize
        height2 = elev_dem.RasterYSize
        gt2 = elev_dem.GetGeoTransform()
    
        minx = gt2[0]
        miny = gt2[3] + width2*gt2[4] + height2*gt2[5] 
        maxx = gt2[0] + width2*gt2[1] + height2*gt2[2]
        maxy = gt2[3] 
            
        minx_tiles.append(minx)
        miny_tiles.append(miny)
        maxx_tiles.append(maxx)
        maxy_tiles.append(maxy)
            
        #create a list of x and y cordinate points for all rasters:
        x_coord = np.arange(minx, maxx, gt2[1])
        y_coord = np.flip(np.arange(miny, maxy, -gt2[5])) #the y coordinate list must be fliped to go from maxy to miny because we want to have coordinates staring from top left corner(the same as data)
 

        #save the data as netCDF file so that the file is smaller (saving processing time):
        if slope==True and aspect==True:            
            out_name = name + '_sae.nc'
        elif slope==True and aspect==False:
            out_name = name + '_se.nc'
        elif slope==False and aspect==True:
            out_name = name + '_ae.nc'
        else:
            out_name = name + '_e.nc'
        
        out_data_file = os.path.join(directory,out_name)
        data_file = Dataset(out_data_file, 'w', format='NETCDF4_CLASSIC')
        
        #creating dimentions for the data:
        polar_x = data_file.createDimension('polar_x', width2)
        polar_y = data_file.createDimension('polar_y', height2)    
    
        #assigning variables for the above dimentions:
        polar_x = data_file.createVariable('polar_x', np.float32, ('polar_x',))
        polar_y = data_file.createVariable('polar_y', np.float32, ('polar_y',))
    
        #create the 2-d variables of actual data (elevation, slope, aspect)
        elevation = data_file.createVariable('elevation', np.float32, ('polar_y','polar_x'))

        #assign atributes to the file (global) and to the variables (variable):
        data_file.description = 'This file is an output of the analysis of DEM ' + str(name) + '. The coordinates have been converted to decimal degree latitude and longitude from the stereographic polar projection, and the slope and aspect have been computed from elevation.'
        data_file.history = 'Created ' + time.ctime(time.time()) 
        data_file.source = 'Computed with python code by KM Sejan, Utrecht University. Original elevation and location data from DEM ' + str(name) + '.'
        polar_x.units = 'meters'
        polar_y.units = 'meters'
        elevation.units = 'meters'

        #assign the data lists into the data file variables:
        polar_x[:] = x_coord
        polar_y[:] = y_coord
        elevation[:] = dem_elev


        while slope==True:
        
            slope_file = str(d).split('.')[0] + '_slope.tif'
            slope_dem_path = os.path.join(directory, slope_file)
            
            #read the DEM:
            slope_dem = gdal.Open(slope_dem_path)
    
            #extract the values:
            slope_band = slope_dem.GetRasterBand(1)
            dem_slope = slope_band.ReadAsArray() #this array is 15000 in length => the y is the length of the array
    
            slope = data_file.createVariable('slope', np.float32, ('polar_y','polar_x'))
            slope.units = 'degrees' 
            slope[:] = dem_slope

        
        while aspect == True:
            aspect_file = str(d).split('.')[0] + '_aspect.tif'
            aspect_dem_path = os.path.join(directory, aspect_file)
            
            #read the DEM:
            aspect_dem = gdal.Open(aspect_dem_path)
                    
            #extract the values:
            aspect_band = aspect_dem.GetRasterBand(1)
            dem_aspect = aspect_band.ReadAsArray() #this array is 15000 in length => the y is the length of the array        
          
            aspect = data_file.createVariable('aspect', np.float32, ('polar_y','polar_x'))
            aspect.units = 'degrees' #trigonometric angle     
            
            aspect[:] = dem_aspect
        
        data_file.close()
                    
    #write the csv file with tile info
    bounds = zip(tiles, minx_tiles, miny_tiles, maxx_tiles, maxy_tiles)
        
    out_name = name + '_tiles_bounds.csv'
    out_data_file = os.path.join(directory,out_name)
        
    with open(out_data_file, "w") as f:
       writer = csv.writer(f)
       writer.writerow(["tiles", "minx", "miny", "maxx", "maxy"])
       for row in bounds:
           writer.writerow(row) 

#%%

'Prepare DEMs'

#GIMP
gdirectory = '/Users/kat/DATA/DEM_GIMP_processed/'
gname = 'GIMP_dem_mosaic_1km'

dem_tif_to_netcdf(gdirectory, gname, slope = True, aspect = True)

#Helm
hdirectory = '/Users/kat/DATA/Helm_DEM_processed/'
hname = 'DEM_GRE_CS_20130826'

dem_tif_to_netcdf(hdirectory, hname, slope = True, aspect = True)

#Arctic
adirectory = '/Users/kat/DATA/ArcticDEM_processed/'
aname = 'ArcticDEM_mosaic_1km'

dem_tif_to_netcdf(adirectory, aname, slope = True, aspect = True)

#%%
"THE END"