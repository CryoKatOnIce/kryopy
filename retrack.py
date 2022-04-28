# -*- coding: utf-8 -*-
"""
Retracking the CryoSat2 Baseline D data.

Created by: Kat Sejan 14th February 2020.
Last edited on: 29th September 2021 by Kat Sejan

"""

def find_valid_window(array, ltresh, htresh, window_size):
    #finds valid window for given value threshold and window length
    #if such window does not exists, looks for shorter windows
    
    import numpy as np
    
    hmask = htresh < array #create mask to find which values of the array are above threshold
    lmask = array > ltresh #create mask to find which values of the array are below threshold
    mask = hmask != lmask
    mask_extended = np.r_[False,mask,False] #just puts false at the begining and at the end
    change_indices = np.flatnonzero(mask_extended[:-1] != mask_extended[1:]) #checking if previous element mask value is different than current element mask value, if so, gives index
    windows_starts = change_indices[::2]
    windows_ends = change_indices[1::2]
    windows_lengths = windows_ends-windows_starts # Take every second element, even then odd, and substract from each other to get length of same mask windows
    if any(windows_lengths >= window_size) == True:
        valid_window = windows_starts[(windows_lengths >= window_size).argmax()] #taking first found window higher than specified size from the list of starting indices
        window_length = window_size
    elif any(windows_lengths >= window_size-int(0.2*window_size)) == True:
        valid_window = windows_starts[(windows_lengths >= window_size-int(0.2*window_size)).argmax()] 
        window_length = window_size-int(0.2*window_size)
    elif any(windows_lengths >= window_size-int(0.4*window_size)) == True:
        valid_window = windows_starts[(windows_lengths >= window_size-int(0.4*window_size)).argmax()]
        window_length = window_size-int(0.4*window_size)
    elif any(windows_lengths >= window_size-int(0.6*window_size)) == True:
        valid_window = windows_starts[(windows_lengths >= window_size-int(0.6*window_size)).argmax()]
        window_length = window_size-int(0.6*window_size)
    elif any(windows_lengths >= window_size-int(0.8*window_size)) == True:
        valid_window = windows_starts[(windows_lengths >= window_size-int(0.8*window_size)).argmax()]
        window_length = window_size-int(0.8*window_size)
    else:
        valid_window = np.abs(array).argmin()
        window_length = 1
    return  valid_window, window_length


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


def file_copy_greenland(path_to_infile, out_directory):
    
    'Library Dependencies'
    from netCDF4 import Dataset
    import os
    import numpy as np
    
    data_infile = Dataset(path_to_infile, "r", format="NETCDF4")

    #extract coordinates of the meauserments form the file:
    longitudes = data_infile.variables['lon_20_ku'][:].astype('float32').tolist() #read in the longitude values for the points
    latitudes = data_infile.variables['lat_20_ku'][:].astype('float32').tolist() #read in the latitude values for the points
    
    #get index of points within Greenland
    in_index, out_index = points_greenland(longitudes, latitudes)
    
    if len(in_index) != 0:
        #get output path
        file_name = path_to_infile.split('/')[-1]
        path_to_outfile = os.path.join(out_directory, file_name)
        
        #output file
        data_outfile = Dataset(path_to_outfile, "w", format="NETCDF4")
        
        #Copy dimensions
        for dname, the_dim in data_infile.dimensions.items():
            if len(the_dim) == len(in_index+out_index):
                len_dim = len(in_index)
            else:
                len_dim = len(the_dim)                
            data_outfile.createDimension(dname, len_dim)
            
                
        # Copy variables
        for v_name, varin in data_infile.variables.items():
            outVar = data_outfile.createVariable(v_name, varin.datatype, varin.dimensions, zlib=True)
            
            # Copy variable attributes
            outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
            
            variablein = varin[:].filled()
            
            if '_FillValue' in varin.ncattrs():
                variableout = np.ma.masked_values(np.array([variablein[il] for il,ll in enumerate(variablein) if il in in_index]), varin[:].fill_value)
            else:
                variableout = np.array([variablein[il] for il,ll in enumerate(variablein) if il in in_index])
            
            outVar[:] = variableout
            
        # close the output file
        data_outfile.close()
    
    data_infile.close()

    return path_to_outfile

def read_L1b(path_to_outfile):
    
    'Library Dependencies'
    from netCDF4 import Dataset
    
    data_file = Dataset(path_to_outfile, "r", format="NETCDF4")

    #extract coordinates of the meauserments form the file:
    longitudes = data_file.variables['lon_20_ku'][:].astype('float32').tolist() #read in the longitude values for the points
    latitudes = data_file.variables['lat_20_ku'][:].astype('float32').tolist() #read in the latitude values for the points
    power_waveform = data_file.variables['pwr_waveform_20_ku'][:]  #ectract the waveforms for each measurement
    
    #extract measurement time and flags:
    time_sat = data_file.variables['time_20_ku'][:] #needed?
    flag_corrections = data_file.variables['flag_cor_status_01'][:] #value 4095 indicates that all of the corrections have been used: 
    time_corrections = data_file.variables['time_cor_01'][:]
    """ Corrections in the flag corrections are as follow in this order: 
        model_dry_called model_wet_called inv_bar_called hf_fluctuations_called 
        iono_gim_called iono_model_called ocean_tide_called ocean_tide_equil_called 
        load_tide_called solid_earth_called pole_tide_called surface_type_called"""
    
    #extracting the corrections values:
    dry_tropo_corr = data_file.variables['mod_dry_tropo_cor_01'][:] #[m] scale_factor: 0.001
    wet_tropo_corr = data_file.variables['mod_wet_tropo_cor_01'][:] #[m] scale_factor: 0.001
    ionosph_corr = data_file.variables['iono_cor_01'][:] #[m] scale_factor: 0.001
    ionosph_gim_corr = data_file.variables['iono_cor_gim_01'][:] #[m] scale_factor: 0.001
    ocean_corr = data_file.variables['ocean_tide_01'][:] #[m] scale_factor: 0.001
    ocean_eqi_corr = data_file.variables['ocean_tide_eq_01'][:] #[m] scale_factor: 0.001
    load_corr = data_file.variables['load_tide_01'][:] #[m] scale_factor: 0.001
    solid_corr = data_file.variables['solid_earth_tide_01'][:] #[m] scale_factor: 0.001
    pole_corr = data_file.variables['pole_tide_01'][:] #[m] scale_factor: 0.001
    
    #extracting other important parameters
    doppler_corr = data_file.variables['dop_cor_20_ku'][:] #[m] scale_factor: 0.001
    window_delay = data_file.variables['window_del_20_ku'][:] #[s] scale_factor: 1.e-12        
    noise_floor = data_file.variables['noise_power_20_ku'][:] #[dB]        
    sat_alt = data_file.variables['alt_20_ku'][:].astype('float32').tolist() #satelite altitude for each measurement point
    
    print('Reading the file...')
    data_file.close()
    
    in_index, out_index = points_greenland(longitudes, latitudes)
    
    "Here we need to add the file copy loop together with the loop for shortening dimentions. "
    
    longitudes = [longitudes[il] for il,ll in enumerate(longitudes) if il in in_index]
    latitudes = [latitudes[il] for il,ll in enumerate(latitudes) if il in in_index]
    power_waveform = [power_waveform[il] for il,ll in enumerate(power_waveform) if il in in_index]
    time_sat = [time_sat[il] for il,ll in enumerate(time_sat) if il in in_index]
    flag_corrections = [flag_corrections[il] for il,ll in enumerate(flag_corrections) if il in in_index]
    time_corrections = [time_corrections[il] for il,ll in enumerate(time_corrections) if il in in_index]
    dry_tropo_corr = [dry_tropo_corr[il] for il,ll in enumerate(dry_tropo_corr) if il in in_index]
    wet_tropo_corr = [wet_tropo_corr[il] for il,ll in enumerate(wet_tropo_corr) if il in in_index]
    ionosph_corr = [ionosph_corr[il] for il,ll in enumerate(ionosph_corr) if il in in_index]
    ionosph_gim_corr = [ionosph_gim_corr[il] for il,ll in enumerate(ionosph_gim_corr) if il in in_index]
    ocean_corr = [ocean_corr[il] for il,ll in enumerate(ocean_corr) if il in in_index]
    ocean_eqi_corr = [ocean_eqi_corr[il] for il,ll in enumerate(ocean_eqi_corr) if il in in_index]
    load_corr = [load_corr[il] for il,ll in enumerate(load_corr) if il in in_index]
    solid_corr = [solid_corr[il] for il,ll in enumerate(solid_corr) if il in in_index]
    pole_corr = [pole_corr[il] for il,ll in enumerate(pole_corr) if il in in_index]
    doppler_corr = [doppler_corr[il] for il,ll in enumerate(doppler_corr) if il in in_index]
    window_delay = [window_delay[il] for il,ll in enumerate(window_delay) if il in in_index]
    noise_floor = [noise_floor[il] for il,ll in enumerate(noise_floor) if il in in_index]
    sat_alt = [sat_alt[il] for il,ll in enumerate(sat_alt) if il in in_index]
    
    
    corrections = [dry_tropo_corr, wet_tropo_corr, ionosph_corr, ionosph_gim_corr, ocean_corr, ocean_eqi_corr, load_corr, solid_corr, pole_corr, doppler_corr]
    corr_labels = ['dry_tropo_corr', 'wet_tropo_corr', 'ionosph_corr', 'ionosph_gim_corr', 'ocean_corr', 'ocean_eqi_corr', 'load_corr', 'solid_corr', 'pole_corr', 'doppler_corr']

        
    return longitudes, latitudes, power_waveform, time_sat, flag_corrections, time_corrections, corrections, corr_labels, window_delay, noise_floor, sat_alt

def interpolate_corrections(power_waveform, time_sat, flag_corrections, time_corrections, corrections, corr_labels):
    
    'Library Dependencies'
    import pandas as pd
    import numpy as np

    
    #Checking the length of correction lists and interpolating for all wavelets, but first making sure that all of the corrections have been applied for a given measurement point:
    flag_check = [True if n == 4095 else False for n in flag_corrections] #value 4095 means alll corrections have been applied for this measurement 
    flag_action = all(flag_check)
    
    if flag_action == True:
        
        flag_ind = [indx for indx,n in enumerate(time_sat) if n in time_corrections] #creating the list of idex for each correction flag based on the time stamps (flags are given for the next 20 measurements usually, but not guaranteed so use of time stamps is necessary)
        length_check = [True if len(name)!=len(power_waveform) and len(name)==len(flag_corrections) else False for name in corrections] #getting the corrections that have to be interpolated (not all corrections are given for every 20 measurements, some are given for each measurement)
        short_corrections = [name for i,name in enumerate(corrections) if length_check[i] == True] #corrections to be interpolated
        short_corr_labels = [label for i,label in enumerate(corr_labels) if length_check[i] == True]
        long_corrections = [name for i,name in enumerate(corrections) if length_check[i] == False] #corrections that don't need interpolation
        long_corr_labels = [label for i,label in enumerate(corr_labels) if length_check[i] == False]
    
        #create a data frame to structure the data and interpolate all in one go:
        df_corrections = pd.concat([pd.Series(x, dtype ='float64') for x in short_corrections], axis=1)
        flag_index = pd.Index(flag_ind, name = 'flag_index')
        df_cor = pd.DataFrame(df_corrections.values, index = flag_index, columns = short_corr_labels)                                                                   
        
        n_wavelets = [0, len(power_waveform)]
        
        new_index = pd.Index(np.arange(n_wavelets[0], n_wavelets[1],1), name = 'interp') 
        expand_corr = df_cor.set_index(flag_index).reindex(new_index)
        
        #i_corr_df = expand_corr.fillna(value = None, method = 'ffill', axis = 0)
        interp_corr_df = expand_corr.fillna(value = None, method = 'ffill', axis = 0)
        
        #Add corrections that did not need to be interpolated to the data frame:
        df_long_corr = pd.concat([pd.Series(x) for x in long_corrections], axis=1)
        for i,n in enumerate(long_corr_labels):
            #i_corr_df[n] = df_long_corr[i]
            interp_corr_df[n] = df_long_corr[i]
            
        "Values should not be sacaled according to the scale_factor."
        
    return interp_corr_df

def wave_stats(wave, smoothing_filter='hann', order = 3, cut_freq = 0.5):
    #wave is the power waveform 
    #smoothing filter is the window filter for signal.firwin() default is hann for retracking 
    #order is the order of the filter, here it is 4th order so =3
    #cut off frequency of the filter is set to 0.5Hz
    
    
    'Library Dependencies'
    import numpy as np
    from scipy import signal

    "Smooth the wavelet by using 4th order zero-phase low-pass filter with cut-off frequency of  0.5 Hz for LRM data:"
    b = signal.firwin(order, cut_freq, width=None, window=smoothing_filter, pass_zero=True, scale=True) #Usinf Hann filter
    a = 1
    smooth_wave = signal.filtfilt( b,a,wave)   
    #smooth_wave = signal.convolve(wave, a, mode='same')
    
    "Find the maximum, mean and minimum power of the wavelet and their indices:"
    #find the maximum and the index
    max_p = max(smooth_wave)
    max_p_index = [indx for indx,i in enumerate(smooth_wave) if i == max_p]
    #find the mean power of the wavelet:
    mean_p = np.mean(smooth_wave)
    #checking min of the wave
    min_p = min(smooth_wave)
    min_p_index = [indx for indx,i in enumerate(smooth_wave) if i == min_p]

    return smooth_wave, max_p, max_p_index, mean_p, min_p


def wave_gates(wave, mean_p, smoothing_filter='hann', order = 3, max_order = 5, cut_freq = 0.3):
    #wave is the power waveform 
    #smoothing filter is the window filter for signal.firwin() default is hann for retracking 
    #order is the order of the filter, here it is 4th order so =3
    #max order is a higher order used to oversmooth the wavelet even more
    #cut off frequency of the filter is set to 0.3Hz for oversmoothing the wavelet
    
    'Library Dependencies'
    from scipy import signal
    
    "Oversmooth the wavelet with 4th order filer to find zero_gate."
    #Note: oversmoothing of the waveform changes the location of the peak. By trial and error decided not to use any lowe than 0.3 Hz cut-off frequency:
    c = signal.firwin(3, 0.3, width=None, window='hann', pass_zero=True, scale=True)
    a = 1
    oversmooth_wave1 = signal.filtfilt( c,a,wave)   
    #oversmooth_wave = signal.convolve(wave, c, mode='same')
    
    "Find the first peak and it's start. The peak is defined as being at least larger than the mean of the wavelet:"
    peaks, peaks_properties = signal.find_peaks(oversmooth_wave1, mean_p)
    if len(peaks) != 0:
        peak_locations = peaks.tolist()
        peaks_prominence = signal.peak_prominences(oversmooth_wave1, peak_locations)
        peaks_lbases = peaks_prominence[1]
        zero_gate = peaks_lbases[0] #this is the start of the leading edge peak
    else:
        zero_gate = 0
    
    "Oversmooth the wavelet with 8th order filter to find peak_gate." 
    #Note: oversmoothing of the waveform changes the location of the peak. By trial and error decided not to use any lowe than 0.3 Hz cut-off frequency:
    c = signal.firwin(5, 0.3, width=None, window='hann', pass_zero=True, scale=True)
    a = 1
    oversmooth_wave2 = signal.filtfilt( c,a,wave)   
    #oversmooth_wave = signal.convolve(wave, c, mode='same')
    
    "Find the first peak and it's end. The peak is defined as being at least larger than the mean of the wavelet:"
    peaks, peaks_properties = signal.find_peaks(oversmooth_wave2, mean_p)
    if len(peaks) !=0:
        peak_locations = peaks.tolist()
        for l in peak_locations:
            if l > zero_gate:
                peak_gate = l  #this is the top of the leading edge, the maxima of this peak
                break
        peaks_heights = peaks_properties['peak_heights']
        peak_gate_height = peaks_heights[0]
    else:
        peak_gate = zero_gate
        peak_gate_height = oversmooth_wave2[zero_gate]
        peaks_heights = [max(oversmooth_wave2)]

    return zero_gate, peak_gate, peak_gate_height, peaks_heights

def interpolate_edge(lew, zero_gate, peak_gate):
    
    'Library Dependencies'
    import pandas as pd
    import numpy as np
    import math

    "Interpolate the leading edge to get more samples adding 10 samples between every point pair, but keep the original indexing (range bins) only adding decimals to the index for interpolated values:"
    lindex = np.arange(math.floor(zero_gate), math.ceil(peak_gate)+1,1).tolist()
    data_df = [lindex, lew]
    df_wave3 = pd.concat([pd.Series(x, dtype='float64') for x in data_df], axis=1)
    df_w = pd.DataFrame(df_wave3.values, columns = ['sample', 'wavelet'])
    
    points_range = [lindex[0], lindex[-1]]
    
    df_w.set_index('sample')
    new_index = pd.Index([round(x,1) for x in np.arange(points_range[0], points_range[1]+1,0.1)], name = 'interp') 
    expand_wave = df_w.set_index('sample').reindex(new_index)
    
    interp_wave_df = expand_wave.interpolate(method = 'linear', axis = 0)
    interp_wave = interp_wave_df['wavelet'].to_list() 
    
    return interp_wave, new_index


def output_retracked_file(path_to_outfile, out_data, alpha):
    
    'Library Dependencies'
    import numpy as np
    from netCDF4 import Dataset
      
    
    "Write an output netCDF4 file with all the calculated data"
    try:
        print('... saving netCDF4 file...')

        data_file = Dataset(path_to_outfile, 'r+', format='NETCDF4_CLASSIC')

        max_power_wavelet = data_file.createVariable('max_power_wavelet', np.float32, ('time_20_ku',))
        mean_power_wavelet = data_file.createVariable('mean_power_wavelet', np.float32, ('time_20_ku',))
        noise_wavelet = data_file.createVariable('noise_wavelet', np.float32, ('time_20_ku',))
        signal_to_noise = data_file.createVariable('signal_to_noise', np.float32, ('time_20_ku',))
        zero_gate_wavelet = data_file.createVariable('zero_gate_wavelet', np.float32, ('time_20_ku',))
        zero_gate_flags = data_file.createVariable('zero_gate_flag', np.float32, ('time_20_ku',))
        peak_gate_wavelet = data_file.createVariable('peak_gate_wavelet', np.float32, ('time_20_ku',))
        peak_gate_flags = data_file.createVariable('peak_gate_flag', np.float32, ('time_20_ku',))
        leading_edge_width = data_file.createVariable('leading_edge_width', np.float32, ('time_20_ku',))
        trailing_edge_slope = data_file.createVariable('trailing_edge_slope', np.float32, ('time_20_ku',))
        backscatter_wavelet = data_file.createVariable('backscatter_wavelet', np.float32, ('time_20_ku',))
        pulse_peakiness = data_file.createVariable('pulse_peakiness', np.float32, ('time_20_ku',))
        max_amplitude_wavelet = data_file.createVariable('max_amplitude_wavelet', np.float32, ('time_20_ku',))
        corrections_sum_wavelet = data_file.createVariable('corrections_sum_wavelet', np.float32, ('time_20_ku',))


        max_power_wavelet.description = 'Maximum power of the wavelet.'
        max_power_wavelet.units = 'count'
        max_power_wavelet.source = 'KM Sejan threshold retracker.'
        mean_power_wavelet.description = 'Mean power of the wavelet.'
        mean_power_wavelet.units = 'count'
        mean_power_wavelet.source = 'KM Sejanthreshold retracker.'
        noise_wavelet.description = 'Noise level in the wavelet, background noise.'
        noise_wavelet.units = 'count'
        noise_wavelet.source = 'KM Sejan threshold retracker.'
        signal_to_noise.description = 'Signal to noise ratio.'
        signal_to_noise.units = 'dB'
        signal_to_noise.source = 'KM Sejan threshold retracker.'
        zero_gate_wavelet.description = 'Zero gate, the begining of the leading edge.'
        zero_gate_wavelet.units = 'range_bin_number'
        zero_gate_wavelet.source = 'KM Sejan threshold retracker.'
        zero_gate_flags.description = 'Flag assigned based on zero gate location, the begining of the leading edge. If zero gate is below the 5th bin rage the flag equal to 1 is raised, else flag is equal 0.'
        zero_gate_flags.units = '0 or 1'
        zero_gate_flags.source = 'KM Sejan threshold retracker.'
        leading_edge_width.description = 'This is the leading edge width computed using the peak gate minus the zero gate, times the bin size [m]'
        leading_edge_width.units = 'meters'
        leading_edge_width.source = 'KM Sejan threshold retracker.'
        trailing_edge_slope.description = 'This is the trailing edge slope computed using the wave between peak gate and the end, and using linear regression.'
        trailing_edge_slope.units = 'unitless, but theoretically it represents the change of return power over one bin.'
        trailing_edge_slope.source = 'KM Sejan threshold retracker.'
        peak_gate_wavelet.description = 'Peak gate, the peak of the leading edge.'
        peak_gate_wavelet.units = 'range_bin_number'
        peak_gate_wavelet.source = 'KM Sejan threshold retracker.'
        peak_gate_flags.description = 'Flag assigned based on peak gate height, the peak of the leading edge. If peak gate is lower than the the maximum peak of the waveform, then flag equal to 1 is raised, else flag is equal 0.'
        peak_gate_flags.units = '0 or 1'
        peak_gate_flags.source = 'KM Sejan threshold retracker.'
        backscatter_wavelet.description = 'Backscatter coefficient, the integrated power of the wavelet.'
        backscatter_wavelet.units = 'count'
        backscatter_wavelet.source = 'KM Sejan threshold retracker.'
        pulse_peakiness.description = 'Pulse peakiness, waveform peakiness.'
        pulse_peakiness.units = ''
        pulse_peakiness.source = 'KM Sejan threshold retracker.'
        max_amplitude_wavelet.description = 'Maximum aplitude of the wavelet.'
        max_amplitude_wavelet.units = 'count'
        max_amplitude_wavelet.source = 'KM Sejan threshold retracker.'
        corrections_sum_wavelet.description = 'Sum of all of the corrections.'
        corrections_sum_wavelet.units = 'meters'
        corrections_sum_wavelet.source = 'KM Sejan based on ESA correction values.'

        
        #assign the data lists into the data file variables:
        max_power_wavelet[:] = out_data[0]
        mean_power_wavelet[:] = out_data[1]
        noise_wavelet[:] = out_data[2]
        signal_to_noise[:] = out_data[3]
        zero_gate_wavelet[:] = out_data[4]
        zero_gate_flags[:] = out_data[5]        
        peak_gate_wavelet[:] = out_data[6]
        peak_gate_flags[:] = out_data[7]
        leading_edge_width[:] = out_data[8]
        trailing_edge_slope[:] = out_data[9]        
        backscatter_wavelet[:] = out_data[10]
        pulse_peakiness[:] = out_data[11]
        max_amplitude_wavelet[:] = out_data[12]
        corrections_sum_wavelet[:] = out_data[13]

        
        for aai,aa in enumerate(alpha):
            
            extension = str(aa).replace('.', '_')
            
            #extract indivifual aa data
            out_data_th = list(out_data[14:])
            out_data_thres = [[sub[aai] for sub in out_data_th[0]], [sub[aai] for sub in out_data_th[1]], [sub[aai] for sub in out_data_th[2]], [sub[aai] for sub in out_data_th[3]]]

            threshold_power_wavelet = data_file.createVariable('threshold_power_wavelet_'+extension, np.float32, ('time_20_ku',))
            threshold_index_wavelet = data_file.createVariable('threshold_index_wavelet_'+extension, np.float32, ('time_20_ku',))
            range_wavelet = data_file.createVariable('range_wavelet_'+extension, np.float32, ('time_20_ku',))
            elevation_wavelet = data_file.createVariable('elevation_wavelet_'+extension, np.float32, ('time_20_ku',))
            
            threshold_power_wavelet.description = 'Threshold power for the threshold of ' + extension + '. '
            threshold_power_wavelet.units = 'count'
            threshold_power_wavelet.source = 'KM Sejan threshold retracker.'
            threshold_index_wavelet.description = 'The index of the threshold power for the retracker at threshold of ' + extension + '. '
            threshold_index_wavelet.units = 'range_bin_number'
            threshold_index_wavelet.source = 'KM Sejan threshold retracker.'
            range_wavelet.description = 'Calculated range from the satleite to the surface at threshold of ' + extension + '. '
            range_wavelet.units = 'm'
            range_wavelet.source = 'KM Sejan threshold retracker.'
            elevation_wavelet.description = 'Computed elevation at the point of measurement at threshold of ' + extension + '. '
            elevation_wavelet.units = 'm'
            elevation_wavelet.source = 'KM Sejan threshold retracker.'

            threshold_power_wavelet[:] = out_data_thres[0]
            threshold_index_wavelet[:] = out_data_thres[1]
            range_wavelet[:] = out_data_thres[2]
            elevation_wavelet[:] = out_data_thres[3]

        
        #finish writing the file by closing it:
        data_file.close()
        
        print('Output file saved.')
    
    except (IOError, RuntimeError, ValueError, IndexError) as e:
        print(e)


def beta_fit(wave, betas, bounds):
    
    'Library Dependencies'
    import numpy as np
    from scipy import optimize
    from scipy import special
    
    'Making the y function as a vector'
    
    wave_indx_vector = np.array(range(0,len(wave)))
    Q = np.zeros(len(wave))
    Q = np.array([i-(betas[2]+0.5*betas[3]) if i >= (betas[2] + 0.5*betas[3]) else 0 for i,q in enumerate(Q)])
    
    
    def model(t, coeffs):
        #y = B1 + B2*(1+B5*Q)*(1/2+ 1/2*(special.erf(np.sqrt(1/2)*((t-B3)/B4))))
       return coeffs[0] + (coeffs[1]*(1+coeffs[4]*Q))*(0.5+(np.float64(1/2)* special.erf(np.sqrt(1/2)*((t-coeffs[2])/coeffs[3]))))
    
    def residuals(coeffs, wavelet, t):
        return wavelet - model(t, coeffs)
    

    fit_output = optimize.least_squares(residuals, betas, args= (wave, wave_indx_vector), bounds= bounds)

    return fit_output


###############################################################################

def tfmra_retrack(path_to_infile, out_directory, alpha = [0.3]):
    #alpha is percentage of the threshold, it's a value between 0 and 1. We use 30% as standard.
    'This function applies retracker to one input file (all waveforms) and returns lists of all calculated data for every waveform within the file.'
    
    'Library Dependencies'
    import numpy as np
    from scipy import stats
    import math
    
    try:
        "Filter points only falling within Greenland and save them into a copy of a original file."
        path_to_outfile = file_copy_greenland(path_to_infile, out_directory)
        
        "Reading in the data from file:"
        longitudes, latitudes, power_waveform, time_sat, flag_corrections, time_corrections, corrections, corr_labels, window_delay, noise_floor, sat_alt = read_L1b(path_to_outfile)


        "Interpolating Corrections to the length of the wavelet and applying scale factor."        
        interp_corr_df = interpolate_corrections(power_waveform, time_sat, flag_corrections, time_corrections, corrections, corr_labels)

        "THE RETRACKING SECTION"
        
        
        max_power = []
        mean_power = []
        noises = []
        snrs = []
        zero_gates = []
        zero_gate_flags = []
        peak_gates = []
        peak_gate_flags = []
        le_width = []
        te_slope = []
        backscatter = []
        peakiness = []
        max_amplitude = []
        threshold_power = []
        threshold_index = []
        sum_corrections = []
        range_measurement = []
        elevations = []
       
    
        print('File contains ' + str(len(power_waveform)) + ' waveforms.')
        print('Start retracking...')
        

        for i,n in enumerate (power_waveform):
            
            "Extract the wavelet:"
            wave = power_waveform[i]
            
            #add back the value of 65535 that has been masked
            mask = np.ma.getmaskarray(wave)
            fill_value = 65535
            
            wave[mask] = fill_value
                            
            smooth_wave, max_p, max_p_index, mean_p, min_p = wave_stats(wave, smoothing_filter='hann', order = 3, cut_freq = 0.5 )
            
            zero_gate, peak_gate, peak_gate_height, peaks_heights = wave_gates(wave, mean_p, smoothing_filter='hann', order = 3, max_order = 5, cut_freq = 0.3)
            
            
            "Check if there are any wavelets that have a zero gate sooner than range bin 5, i.e. if there are any possibly cut wavelets or bad wavelets:"
            #if zero_gate < 5: print('Zero gate at < 5, wavelet ' + str(i) + ' must be rejected. The peak gate is at: ' + str(peak_gate) + '. The minimum is at: ' + str(min_p_index[0]) + '.') #we also check where the wavelet minimum is and the leading edge ends to see the shape of the wavelet
            if zero_gate <5: 
                zero_gate_flag = 1
            else: 
                zero_gate_flag = 0
            
            "Check if the peak_gate (i.e. the maxima of leading edge is the highest peak in the waveform. If not, then flag the waveform.)"
            if peak_gate_height < max(peaks_heights):
                peak_gate_flag = 1
            else:
                peak_gate_flag = 0
            
            """Calculate noise based on 5 range bins before the zero gate (5 bin window with gradiens lower than 1000 and higher than -1000, ie. straightes part between begining and xero-gate). 
            In case that the zero gate is before range bin 5 of the wavelet, then take all the bins before zero gate. In case zero gate is at 0, then take zero gate as noise level."""
            if zero_gate == 0:
                noise = np.mean(smooth_wave[zero_gate])
            else:                            
                noise_range = 5
                min_gradient = -1000 #largest negative change in gradient acceptable
                max_gradient = 1000 #largest positive chnage in gradient acceptable
                noise_gate, window_length = find_valid_window(np.gradient(smooth_wave[0:int(zero_gate)+1]), min_gradient, max_gradient, noise_range)
                noise = np.mean(smooth_wave[noise_gate:int(noise_gate + window_length)])

                
            "Compute signal to noise ratio (SNR), and check if it is lower than 5 db:"
            snr = math.log10(max_p/max(noise,1))*10 #[dB]
            snr_check = [False if snr<5 else True] 

                                                            
            "Calculate the backscatter coefficient:"
            integrating_power= np.cumsum(smooth_wave) #integrated power (sum of all the values of wavelet)
            b_coeff = integrating_power[-1] #backscatter coefficient is the integrated power of the wavelet
            
            "Calculate the Pulse Peakiness:"
            pp = max_p/b_coeff
            
            "Extract the leading edge:"
            lew = smooth_wave[math.floor(zero_gate): math.ceil(peak_gate)+1] 
            
            "Find the trailing edge slope"
            tereg = stats.linregress(range(0, len(smooth_wave[peak_gate:])), smooth_wave[peak_gate:]) #extracting here the trailing edge part of the wave and computing linear regression for it.
            tes = tereg[0] #taking the computed slope value. 
            
            "Interpolate the leading edge to get more samples adding 10 samples between every point pair, but keep the original indexing (range bins) only adding decimals to the index for interpolated values:"
            interp_wave, new_index = interpolate_edge(lew, zero_gate, peak_gate)
            
            "Cmpute maximum Amplitude:"
            amp = max(interp_wave)
            
            "Satellite altitude"
            satelite_alt = sat_alt[i]
                       
            "Get the values of corrections for this wavelet."
            corrections_wave = [interp_corr_df.dry_tropo_corr[i], interp_corr_df.wet_tropo_corr[i], interp_corr_df.ionosph_gim_corr[i], interp_corr_df.load_corr[i], interp_corr_df.solid_corr[i], interp_corr_df.pole_corr[i]]

            
            "Calculate range:"
            rb = 0.4684 #This is Range bin sample [m], it is different in SIN= 0.2342 [m]
            c = 299792458 #[m/s]
            w_delay = window_delay[i]
            ref_position = len(wave)/2 #it is still 64 as in baseline B
            corr_sum = sum(corrections_wave)
            
            edge_width = (peak_gate - zero_gate)*rb
            
            "Writing the wavelet data into a list per variable."
            max_power.append(np.float32(round(max_p,4)))
            mean_power.append(np.float32(round(mean_p, 4)))
            noises.append(np.float32(round(noise, 4)))
            snrs.append(np.float32(round(snr, 4)))
            zero_gates.append(np.float32(round(zero_gate, 4)))
            zero_gate_flags.append(int(zero_gate_flag))
            peak_gates.append(np.float32(round(peak_gate, 4)))
            peak_gate_flags.append(int(peak_gate_flag))
            le_width.append(np.float32(round(edge_width, 4)))
            te_slope.append(np.float32(round(tes, 4)))
            backscatter.append(np.float32(round(b_coeff, 4)))
            peakiness.append(np.float32(round(pp, 4)))
            max_amplitude.append(np.float32(round(amp, 4)))
            sum_corrections.append(np.float32(round(corr_sum, 4)))


            
            th_power = []
            th_index = []
            r_measurement = []
            elevs = []

            for ai, a in enumerate(alpha): #where alpha is a list of alphas [0.1, 0.2, 0.3, 0.4, 0.5]
                "Calculate threshold power:"
                
                if noise < amp:
                    thold_power = noise + alpha[ai]*(amp-noise)
                                                                
                    "Find threshold index"
                    indx_thp_interwave = [i for i,m in enumerate(interp_wave) if m >= thold_power][0] #first bin above the treshold in interpolated wave list index
                    value_indx_thp = interp_wave[indx_thp_interwave] #power at the first bin above the treshold
                    indx_thp_below = [i for i,m in enumerate(interp_wave) if m >= thold_power][0] - 1 #last bin below the treshold in interpolated wave list index
                    value_indx_thp_before = interp_wave[indx_thp_below] #Power at the last bin below the threshold
                    indx_thp_deciaml_below = [new_index[i] for i,m in enumerate(interp_wave) if m >= thold_power][0] - 0.1 #last bin below the treshold in interpolated index (decimal)
                    
                    thp_int_index = ((thold_power-value_indx_thp_before)/(value_indx_thp-value_indx_thp_before)) + indx_thp_deciaml_below #previous version (commented above): + indx_val_before
                                            
                    R = 0.5*c*w_delay + rb*(thp_int_index+1-ref_position) + corr_sum
                    elevation = satelite_alt - R
                                        
                    th_power.append(np.float32(round(thold_power, 4)))
                    th_index.append(np.float32(round(thp_int_index, 4)))
                    r_measurement.append(np.float32(round(R, 4)))
                    elevs.append(np.float32(round(elevation, 4)))
                else:
                    th_power.append(-9999)
                    th_index.append(-9999)
                    r_measurement.append(-9999)
                    elevs.append(-9999)

            
            threshold_power.append(th_power)
            threshold_index.append(th_index)
            range_measurement.append(r_measurement)
            elevations.append(elevs)

        
        out_data = [max_power, mean_power, noises, snrs, zero_gates, 
            zero_gate_flags, peak_gates, peak_gate_flags, le_width, 
            te_slope, backscatter, peakiness, max_amplitude, sum_corrections,
            threshold_power, threshold_index, range_measurement, elevations]

        
        output_retracked_file(path_to_outfile, out_data, alpha) #save to apropriate out directory
 
        
        print('... retracking finished.')
        
            
    except (IOError, RuntimeError, ValueError, IndexError) as e:
        print(e)
    

            

def ocog_retrack(path_to_infile, out_directory, alpha = [0.3]):
        #alpha is percentage of the threshold, it's a value between 0 and 1. We use 30% as standard.
    'This function applies retracker to one input file (all waveforms) and returns lists of all calculated data for every waveform within the file.'
    
    'Library Dependencies'
    import numpy as np
    from scipy import stats
    import math
    
    try:
        "Filter points only falling within Greenland and save them into a copy of a original file."
        path_to_outfile = file_copy_greenland(path_to_infile, out_directory)

        
        "Reading in the data from file:"
        longitudes, latitudes, power_waveform, time_sat, flag_corrections, time_corrections, corrections, corr_labels, window_delay, noise_floor, sat_alt = read_L1b(path_to_outfile)


        "Interpolating Corrections to the length of the wavelet and applying scale factor."        
        interp_corr_df = interpolate_corrections(power_waveform, time_sat, flag_corrections, time_corrections, corrections, corr_labels)

        "THE RETRACKING SECTION"
             
        
        max_power = []
        mean_power = []
        noises = []
        snrs = []
        zero_gates = []
        zero_gate_flags = []
        peak_gates = []
        peak_gate_flags = []
        le_width = []
        te_slope = []
        backscatter = []
        peakiness = []
        max_amplitude = []
        threshold_power = []
        threshold_index = []
        sum_corrections = []
        range_measurement = []
        elevations = []
    
        print('File contains ' + str(len(power_waveform)) + ' waveforms.')
        print('Start retracking...')
        
        "THE RETRACKING SECTION"
        for i,n in enumerate (power_waveform):
            
            "Extract the wavelet:"
            wave = power_waveform[i]
            
            #add back the value of 65535 that has been masked
            mask = np.ma.getmaskarray(wave)
            fill_value = 65535
            
            wave[mask] = fill_value
            
         
            smooth_wave, max_p, max_p_index, mean_p, min_p = wave_stats(wave, smoothing_filter='hann', order = 3, cut_freq = 0.5 )
            
            zero_gate, peak_gate, peak_gate_height, peaks_heights = wave_gates(wave, mean_p, smoothing_filter='hann', order = 3, max_order = 5, cut_freq = 0.3)
            
            
            "Check if there are any wavelets that have a zero gate sooner than range bin 5, i.e. if there are any possibly cut wavelets or bad wavelets:"
            #if zero_gate < 5: print('Zero gate at < 5, wavelet ' + str(i) + ' must be rejected. The peak gate is at: ' + str(peak_gate) + '. The minimum is at: ' + str(min_p_index[0]) + '.') #we also check where the wavelet minimum is and the leading edge ends to see the shape of the wavelet
            if zero_gate <5: 
                zero_gate_flag = 1
            else: 
                zero_gate_flag = 0
            
            "Check if the peak_gate (i.e. the maxima of leading edge is the highest peak in the waveform. If not, then flag the waveform.)"
            if peak_gate_height < max(peaks_heights):
                peak_gate_flag = 1
            else:
                peak_gate_flag = 0
            
                        
            
            """Calculate noise based on 5 range bins before the zero gate (5 bin window with gradiens lower than 1000 and higher than -1000, ie. straightes part between begining and xero-gate). 
            In case that the zero gate is before range bin 5 of the wavelet, then take all the bins before zero gate. In case zero gate is at 0, then take zero gate as noise level."""
            if zero_gate == 0:
                noise = np.mean(smooth_wave[zero_gate])
            else:                            
                noise_range = 5
                min_gradient = -1000 #largest negative change in gradient acceptable
                max_gradient = 1000 #largest positive chnage in gradient acceptable
                noise_gate, window_length = find_valid_window(np.gradient(smooth_wave[0:int(zero_gate)+1]), min_gradient, max_gradient, noise_range)
                noise = np.mean(smooth_wave[noise_gate:int(noise_gate + window_length)])

            
            "Compute signal to noise ratio (SNR), and check if it is lower than 5 db:"
            snr = math.log10(max_p/max(noise,1))*10 #[dB] 
            snr_check = [False if snr<5 else True] 
            

            "Calculate the backscatter coefficient:"
            integrating_power= np.cumsum(smooth_wave) #integrated power (sum of all the values of wavelet)
            b_coeff = integrating_power[-1] #backscatter coefficient is the integrated power of the wavelet
            
            "Calculate the Pulse Peakiness:"
            pp = max_p/b_coeff
              
            "Cmpute maximum Amplitude:"
            aln = 4 #number of aliased bins
            center_gravity = (sum([indx*f**2 for indx,f in enumerate(smooth_wave[aln:128-aln])]))/(sum([g**2 for g in smooth_wave[aln:128-aln]]))
            amp = np.sqrt((sum([f**4 for f in smooth_wave[aln:128-aln]]))/(sum([g**2 for g in smooth_wave[aln:128-aln]])))
            wave_width = (sum([f**2 for f in smooth_wave[aln:128-aln]]))**2/(sum([g**4 for g in smooth_wave[aln:128-aln]]))
            
            
            "Find the index/sample number/range bin for the thold_power:"
            if center_gravity > (wave_width/2):
                le_index = round(center_gravity - (wave_width/2), 1)
            else: 
                le_index = round(center_gravity - center_gravity , 1)
            
            "Find the trailing edge slope"
            tereg = stats.linregress(range(0, len(smooth_wave[int(center_gravity):])), smooth_wave[int(center_gravity):]) #extracting here the trailing edge part of the wave and computing linear regression for it.
            tes = tereg[0] #taking the computed slope value. 

            
            "Satellite altitude"
            satelite_alt = sat_alt[i]
                       
            "Get the values of corrections for this wavelet."
            corrections_wave = [interp_corr_df.dry_tropo_corr[i], interp_corr_df.wet_tropo_corr[i], interp_corr_df.ionosph_gim_corr[i], interp_corr_df.load_corr[i], interp_corr_df.solid_corr[i], interp_corr_df.pole_corr[i]]

            
            "Calculate range:"
            rb = 0.4684 #This is Range bin sample [m], it is different in SIN= 0.2342 [m]
            c = 299792458 #[m/s]
            w_delay = window_delay[i]
            ref_position = len(wave)/2 #it is still 64 as in baseline B
            corr_sum = sum(corrections_wave)
            
            edge_width = (wave_width/2)*rb
            
            "Writing the wavelet data into a list per variable."
            max_power.append(np.float32(round(max_p,4)))
            mean_power.append(np.float32(round(mean_p, 4)))
            noises.append(np.float32(round(noise, 4)))
            snrs.append(np.float32(round(snr, 4)))
            zero_gates.append(np.float32(round(zero_gate, 4)))
            zero_gate_flags.append(int(zero_gate_flag))
            peak_gates.append(np.float32(round(peak_gate, 4)))
            peak_gate_flags.append(int(peak_gate_flag))
            le_width.append(np.float32(round(edge_width, 4)))
            te_slope.append(np.float32(round(tes, 4)))
            backscatter.append(np.float32(round(b_coeff, 4)))
            peakiness.append(np.float32(round(pp, 4)))
            max_amplitude.append(np.float32(round(amp, 4)))
            sum_corrections.append(np.float32(round(corr_sum, 4)))
            
            th_power = []
            th_index = []
            r_measurement = []
            elevs = []

            for ai, a in enumerate(alpha): #where alpha is a list of alphas [0.1, 0.2, 0.3, 0.4, 0.5]
                "Calculate threshold power:"
                thold_power = noise + alpha[ai]*(amp-noise)
                
                "Interpolate leading edge"
                if thold_power > max(smooth_wave[math.floor(le_index):math.ceil(center_gravity)+1]):
                    max_p_i = max_p_index[0]
                    if max_p_i < le_index:
                        lew = smooth_wave[math.floor(0):math.ceil(max_p_i)+1]
                        interp_wave, new_index = interpolate_edge(lew, 0, max_p_i)
                    else:
                        lew = smooth_wave[math.floor(le_index):math.ceil(max_p_i)+1]
                        interp_wave, new_index = interpolate_edge(lew, le_index, max_p_i)
                else:
                    lew = smooth_wave[math.floor(le_index):math.ceil(center_gravity)+1] 
                    interp_wave, new_index = interpolate_edge(lew, le_index, center_gravity)
                        
                
                "Find threshold index"
                indx_thp_interwave = [iw for iw,m in enumerate(interp_wave) if m >= thold_power][0] #first bin above the treshold in interpolated wave list index
                value_indx_thp = interp_wave[indx_thp_interwave] #power at the first bin above the treshold
                indx_thp_below = [iw for iw,m in enumerate(interp_wave) if m >= thold_power][0] - 1 #last bin below the treshold in interpolated wave list index
                value_indx_thp_before = interp_wave[indx_thp_below] #Power at the last bin below the threshold
                indx_thp_deciaml_below = [new_index[iw] for iw,m in enumerate(interp_wave) if m >= thold_power][0] - 0.1 #last bin below the treshold in interpolated index (decimal)
                
                
                thp_int_index = ((thold_power-value_indx_thp_before)/(value_indx_thp-value_indx_thp_before)) + indx_thp_deciaml_below #previous version (commented above): + indx_val_before

                                        
                R = 0.5*c*w_delay + rb*(thp_int_index+1-ref_position) + corr_sum
                elevation = satelite_alt - R
                                    
                th_power.append(np.float32(round(thold_power, 4)))
                th_index.append(np.float32(round(thp_int_index, 4)))
                r_measurement.append(np.float32(round(R, 4)))
                elevs.append(np.float32(round(elevation, 4)))
            
            
            threshold_power.append(th_power)
            threshold_index.append(th_index)
            range_measurement.append(r_measurement)
            elevations.append(elevs)
                

        out_data = [max_power, mean_power, noises, snrs, zero_gates, 
            zero_gate_flags, peak_gates, peak_gate_flags, le_width, 
            te_slope, backscatter, peakiness, max_amplitude, sum_corrections,
            threshold_power, threshold_index, range_measurement, elevations]

        
        output_retracked_file(path_to_outfile, out_data, alpha) #save to apropriate out directory

            
        print('... retracking finished.')
        
          
    except (IOError, RuntimeError, ValueError, IndexError) as e:
        print(e)
    



def Beta5_tfmra_retrack(path_to_infile, out_directory, alpha = [0.3]):
        #alpha is percentage of the threshold, it's a value between 0 and 1. We use 30% as standard.
    'This function applies retracker to one input file (all waveforms) and returns lists of all calculated data for every waveform within the file.'
    
    'Library Dependencies'
    import numpy as np
    from scipy import stats
    import math
    
    try:
        "Filter points only falling within Greenland and save them into a copy of a original file."
        path_to_outfile = file_copy_greenland(path_to_infile, out_directory)


        "Reading in the data from file:"
        longitudes, latitudes, power_waveform, time_sat, flag_corrections, time_corrections, corrections, corr_labels, window_delay, noise_floor, sat_alt = read_L1b(path_to_outfile)


        "Interpolating Corrections to the length of the wavelet and applying scale factor."        
        interp_corr_df = interpolate_corrections(power_waveform, time_sat, flag_corrections, time_corrections, corrections, corr_labels)

        "THE RETRACKING SECTION"
            
        max_power = []
        mean_power = []
        noises = []
        snrs = []
        zero_gates = []
        zero_gate_flag = []
        peak_gates = []
        peak_gate_flag = []
        le_width = []
        te_slope = []
        backscatter = []
        peakiness = []
        max_amplitude = []
        threshold_power = []
        threshold_index = []
        sum_corrections = []
        range_measurement = []
        elevations = []
        
    
        print('File contains ' + str(len(power_waveform)) + ' waveforms.')
        print('Start retracking...')
        
        "THE RETRACKING SECTION"
        for i,n in enumerate (power_waveform):
            
            "Extract the wavelet:"
            wave = power_waveform[i]
            
            #add back the value of 65535 that has been masked
            mask = np.ma.getmaskarray(wave)
            fill_value = 65535
            
            wave[mask] = fill_value
            
            
            smooth_wave, max_p, max_p_index, mean_p, min_p = wave_stats(wave, smoothing_filter='hann', order = 3, cut_freq = 0.5 )
            
            zero_gate, peak_gate, peak_gate_height, peaks_heights = wave_gates(wave, mean_p, smoothing_filter='hann', order = 3, max_order = 5, cut_freq = 0.3)
            
                        
            """Calculate noise based on 5 range bins before the zero gate (5 bin window with gradiens lower than 1000 and higher than -1000, ie. straightes part between begining and xero-gate). 
            In case that the zero gate is before range bin 5 of the wavelet, then take all the bins before zero gate. In case zero gate is at 0, then take zero gate as noise level."""
            if zero_gate == 0:
                noise = np.mean(smooth_wave[zero_gate])
            else:                            
                noise_range = 5
                min_gradient = -1000 #largest negative change in gradient acceptable
                max_gradient = 1000 #largest positive chnage in gradient acceptable
                noise_gate, window_length = find_valid_window(np.gradient(smooth_wave[0:int(zero_gate)+1]), min_gradient, max_gradient, noise_range)
                noise = np.mean(smooth_wave[noise_gate:int(noise_gate + window_length)])
            
            
            "Compute signal to noise ratio (SNR), and check if it is lower than 5 db:"
            snr = math.log10(max_p/max(noise,1))*10 #[dB]
            snr_check = [False if snr<5 else True] 
            
            #currently indtead of rejecting the wavefrom, just output the snr in the putput file
            bad_snr = []
            if snr_check==False:
                bad_snr.append(i) 
            
            "Calculate the backscatter coefficient:"
            integrating_power= np.cumsum(smooth_wave) #integrated power (sum of all the values of wavelet)
            b_coeff = integrating_power[-1] #backscatter coefficient is the integrated power of the wavelet
            
            "Calculate the Pulse Peakiness:"
            pp = max_p/b_coeff
            
            "Extract the leading edge:"
            lew = smooth_wave[math.floor(zero_gate): math.ceil(peak_gate)+1] 
            
            "Find the trailing edge slope"
            tereg = stats.linregress(range(0, len(smooth_wave[peak_gate:])), smooth_wave[peak_gate:]) #extracting here the trailing edge part of the wave and computing linear regression for it.
            tes = tereg[0] #taking the computed slope value. 
            
            
            "Interpolate the leading edge to get more samples adding 10 samples between every point pair, but keep the original indexing (range bins) only adding decimals to the index for interpolated values:"
            interp_wave, new_index = interpolate_edge(lew, zero_gate, peak_gate)
            
            "Cmpute maximum Amplitude:"
            amp = max(interp_wave)
            
            edge_width = (peak_gate - zero_gate)
            
            """5-Beta retracker with TFMRA initial betas"""
            
            B1 = np.float32(noise) #noise where nuber of samples used to estimate it depends on the location of the leading edge (above)
            B2 = np.float32(amp-noise) #wavelet amplitude - noise 
            B3 = np.float32(peak_gate - edge_width/2) #leading edge mid point 
            B4 = np.float32(edge_width) #leading edge width 
            B5 = np.float32(tes/amp) #trailing edge slope, set to zero in many cases (D.C. Slobbe, 2016) but here we use simple linear least squares regression (above) to estimate the slope and then divide by amplitude (because B2*B5*t)
            
            
            betas = [B1, B2, B3, B4, B5]
            lb  = [0,  0, zero_gate, 0, -0.2];
            ub  = [amp/2, amp, peak_gate, edge_width, -0.0001];
            bounds = (lb, ub)


            "We apply the fit over 5 iterations with increasing bounds. If still unacceptable, then flag wavelet for removal."
            wave_flag = 0
            fit_output = []
            
            for ii in range(2):
              for attempt in range(5):
                try:
                  fit_output = beta_fit(wave, betas, bounds)
                except:
                  #print(attempt)
                  it = attempt+1
                  lb  = [0,  0, zero_gate-(it*2), 0, -0.2];
                  ub  = [(amp/2) + (amp/(it/10)), amp + (amp/(it/10)), peak_gate + (it*2), edge_width + (it*2), -0.0001];
                  bounds = (lb, ub)
                else:
                  break
              else:
                # we failed all the attempts - deal with the consequences. 
                wave_flag = 1

            "Now we make sure we do not append any wavelets with bad fit"
            if wave_flag == 1:
                max_power.append(-9999)
                mean_power.append(-9999)
                noises.append(-9999)
                snrs.append(-9999)
                zero_gates.append(-9999)
                peak_gates.append(-9999)
                le_width.append(-9999)
                te_slope.append(-9999)
                backscatter.append(-9999)
                peakiness.append(-9999)
                max_amplitude.append(-9999)
                sum_corrections.append(-9999)
                zero_gate_flag.append(1)
                peak_gate_flag.append(1)
                
                th_power = []
                th_index = []
                r_measurement = []
                elevs = []
                for ai, a in enumerate(alpha):
                    th_power.append(-9999)
                    th_index.append(-9999)
                    r_measurement.append(-9999)
                    elevs.append(-9999)                
                
                threshold_power.append(th_power)
                threshold_index.append(th_index)
                range_measurement.append(r_measurement)
                elevations.append(elevs)

                
            else:

                noise_B1 = fit_output['x'][0]
                amplitude_B2 = fit_output['x'][1]
                mid_le_B3 = fit_output['x'][2]
                edge_width_B4 = fit_output['x'][3]
                tes_B5 = fit_output['x'][4]*fit_output['x'][2]
                
                snr_new = math.log10(max_p/noise_B1)*10 #[dB]
                zero_gate_new = mid_le_B3 - (edge_width_B4/2)
                peak_gate_new = mid_le_B3 + (edge_width_B4/2)
                               
                "Check if there are any wavelets that have a zero gate sooner than range bin 5, i.e. if there are any possibly cut wavelets or bad wavelets:"
                #if zero_gate < 5: print('Zero gate at < 5, wavelet ' + str(i) + ' must be rejected. The peak gate is at: ' + str(peak_gate) + '. The minimum is at: ' + str(min_p_index[0]) + '.') #we also check where the wavelet minimum is and the leading edge ends to see the shape of the wavelet
                if zero_gate <5: 
                    zero_gate_flag.append(1)
                else: 
                    zero_gate_flag.append(0)
                
                "Check if the peak_gate (i.e. the maxima of leading edge is the highest peak in the waveform. If not, then flag the waveform.)"
                if peak_gate_height < max(peaks_heights):
                    peak_gate_flag.append(1)
                else:
                    peak_gate_flag.append(0)
    
                satelite_alt = sat_alt[i]
                
                "Get the values of corrections for this wavelet."
                corrections_wave = [interp_corr_df.dry_tropo_corr[i], interp_corr_df.wet_tropo_corr[i], interp_corr_df.ionosph_gim_corr[i], interp_corr_df.load_corr[i], interp_corr_df.solid_corr[i], interp_corr_df.pole_corr[i]]
    
                "Calculate range:"
                rb = 0.4684 #This is Range bin sample [m], it is different in SIN= 0.2342 [m]
                c = 299792458 #[m/s]
                w_delay = window_delay[i]
                ref_position = len(wave)/2 #it is still 64 as in baseline B
                corr_sum = sum(corrections_wave)
                
                le_width_m = edge_width_B4*rb #leading edge width in meters
     
                "Writing the wavelet data into a list per variable."
                max_power.append(np.float32(round(max_p,4)))
                mean_power.append(np.float32(round(mean_p, 4)))
                noises.append(np.float32(round(noise_B1, 4)))
                snrs.append(np.float32(round(snr_new, 4)))
                zero_gates.append(np.float32(round(zero_gate_new, 4)))
                peak_gates.append(np.float32(round(peak_gate_new, 4)))
                le_width.append(np.float32(round(le_width_m, 4)))
                te_slope.append(np.float32(round(tes_B5, 4)))
                backscatter.append(np.float32(round(b_coeff, 4)))
                peakiness.append(np.float32(round(pp, 4)))
                max_amplitude.append(np.float32(round(amplitude_B2, 4)))
                sum_corrections.append(np.float32(round(corr_sum, 4)))
                
                
                th_power = []
                th_index = []
                r_measurement = []
                elevs = []
    
                for ai, a in enumerate(alpha): #where alpha is a list of alphas [0.1, 0.2, 0.3, 0.4, 0.5] 
                    thp_int_index = fit_output['x'][2] - ((fit_output['x'][3]*(1-alpha[ai]))/2)
                    
                    if thp_int_index > peak_gate:
                        interp_thp_indx = int(thp_int_index)
                        thold_power = wave[interp_thp_indx]
                    elif thp_int_index < 0:
                        interp_thp_indx = 0
                        thold_power = wave[interp_thp_indx]
                    elif thp_int_index < zero_gate:
                        interp_thp_indx = int(thp_int_index)
                        thold_power = wave[interp_thp_indx]
                    else:
                        interp_thp_indx = int(round(thp_int_index - math.floor(zero_gate), 1) * 10)
                        thold_power = interp_wave[interp_thp_indx]
    
                        
                    R = 0.5*c*w_delay + rb*(thp_int_index+1-ref_position) + corr_sum
                        
                    "Substract from satelite altitude:"
                    elevation = satelite_alt - R
                         
                    th_power.append(np.float32(round(thold_power, 4)))
                    th_index.append(np.float32(round(thp_int_index, 4)))
                    r_measurement.append(np.float32(round(R, 4)))
                    elevs.append(np.float32(round(elevation, 4)))
                
                
                threshold_power.append(th_power)
                threshold_index.append(th_index)
                range_measurement.append(r_measurement)
                elevations.append(elevs)


        out_data = [max_power, mean_power, noises, snrs, zero_gates, 
            zero_gate_flag, peak_gates, peak_gate_flag, le_width, 
            te_slope, backscatter, peakiness, max_amplitude, sum_corrections,
            threshold_power, threshold_index, range_measurement, elevations]

        
        output_retracked_file(path_to_outfile, out_data, alpha) #save to apropriate out directory
                
        
        print('... retracking finished.')        
            
    except (IOError, RuntimeError, ValueError, IndexError) as e:
        print(e)
        


def Beta5_ocog_retrack(path_to_infile, out_directory, alpha = [0.3]):
    #alpha is percentage of the threshold, it's a value between 0 and 1. We use 30% as standard.
    'This function applies retracker to one input file (all waveforms) and returns lists of all calculated data for every waveform within the file.'
    
    'Library Dependencies'
    import numpy as np
    from scipy import stats
    import math
    
    try:
        "Filter points only falling within Greenland and save them into a copy of a original file."
        path_to_outfile = file_copy_greenland(path_to_infile, out_directory)


        "Reading in the data from file:"
        longitudes, latitudes, power_waveform, time_sat, flag_corrections, time_corrections, corrections, corr_labels, window_delay, noise_floor, sat_alt = read_L1b(path_to_outfile)


        "Interpolating Corrections to the length of the wavelet and applying scale factor."        
        interp_corr_df = interpolate_corrections(power_waveform, time_sat, flag_corrections, time_corrections, corrections, corr_labels)

        "THE RETRACKING SECTION"
        
        max_power = []
        mean_power = []
        noises = []
        snrs = []
        zero_gates = []
        zero_gate_flag = []
        peak_gates = []
        peak_gate_flag = []
        le_width = []
        te_slope = []
        backscatter = []
        peakiness = []
        max_amplitude = []
        threshold_power = []
        threshold_index = []
        sum_corrections = []
        range_measurement = []
        elevations = []
         
    
        print('File contains ' + str(len(power_waveform)) + ' waveforms.')
        print('Start retracking...')
        
        "THE RETRACKING SECTION"
        for i,n in enumerate (power_waveform):
            
            "Extract the wavelet:"
            wave = power_waveform[i]
            
            #add back the value of 65535 that has been masked
            mask = np.ma.getmaskarray(wave)
            fill_value = 65535
            
            wave[mask] = fill_value
            
            smooth_wave, max_p, max_p_index, mean_p, min_p = wave_stats(wave, smoothing_filter='hann', order = 3, cut_freq = 0.5 )
            
            zero_gate, peak_gate, peak_gate_height, peaks_heights = wave_gates(wave, mean_p, smoothing_filter='hann', order = 3, max_order = 5, cut_freq = 0.3)
            
                        
            """Calculate noise based on 5 range bins before the zero gate (5 bin window with gradiens lower than 1000 and higher than -1000, ie. straightes part between begining and xero-gate). 
            In case that the zero gate is before range bin 5 of the wavelet, then take all the bins before zero gate. In case zero gate is at 0, then take zero gate as noise level."""
            if zero_gate == 0:
                noise = np.mean(smooth_wave[zero_gate])
            else:                            
                noise_range = 5
                min_gradient = -1000 #largest negative change in gradient acceptable
                max_gradient = 1000 #largest positive chnage in gradient acceptable
                noise_gate, window_length = find_valid_window(np.gradient(smooth_wave[0:int(zero_gate)+1]), min_gradient, max_gradient, noise_range)
                noise = np.mean(smooth_wave[noise_gate:int(noise_gate + window_length)])

            
            "Compute signal to noise ratio (SNR), and check if it is lower than 5 db:"
            snr = math.log10(max_p/max(noise,1))*10 #[dB]
            snr_check = [False if snr<5 else True] 
            
            #currently indtead of rejecting the wavefrom, just output the snr in the putput file
            bad_snr = []
            if snr_check==False:
                bad_snr.append(i) 
            
            
            "Calculate the backscatter coefficient:"
            integrating_power= np.cumsum(smooth_wave) #integrated power (sum of all the values of wavelet)
            b_coeff = integrating_power[-1] #backscatter coefficient is the integrated power of the wavelet
            
            "Calculate the Pulse Peakiness:"
            pp = max_p/b_coeff
              
            "Cmpute maximum Amplitude:"
            center_gravity = (sum([indx*f**2 for indx,f in enumerate(smooth_wave)]))/(sum([g**2 for g in smooth_wave]))
            amp = np.sqrt((sum([f**4 for f in smooth_wave]))/(sum([g**2 for g in smooth_wave])))
            wave_width = (sum([f**2 for f in smooth_wave]))**2/(sum([g**4 for g in smooth_wave]))
            
            edge_width = wave_width/2
            
            "Find the index/sample number/range bin for the thold_power:"
            if center_gravity > (wave_width/2):
                le_index = round(center_gravity - (wave_width/2), 1)
            else: 
                le_index = round(center_gravity - center_gravity , 1)
            

            "Interpolate the leading edge"    
            lew = smooth_wave[math.floor(le_index):math.ceil(center_gravity)+1] 
            
            interp_wave, new_index = interpolate_edge(lew, le_index, center_gravity)

                        
            tereg = stats.linregress(range(0, len(smooth_wave[int(center_gravity):])), smooth_wave[int(center_gravity):]) #extracting here the trailing edge part of the wave and computing linear regression for it.
            tes = tereg[0] #taking the computed slope value. 
            

            """5-Beta retracker model building with OCOG initial beetas"""
            
            B1 = np.float32(noise) #noise where nuber of samples used to estimate it depends on the location of the leading edge (above)
            B2 = np.float32(amp-noise) #wavelet amplitude - noise 
            B3 = np.float32(center_gravity - edge_width/2) #leading edge mid point 
            B4 = np.float32(edge_width) #leading edge width 
            B5 = np.float32(tes/amp) #trailing edge slope, set to zero in many cases (D.C. Slobbe, 2016) but here we use simple linear least squares regression (above) to estimate the slope and then divide by amplitude (because B2*B5*t)
            
           
            betas = [B1, B2, B3, B4, B5]
            lb  = [0,  0, center_gravity - edge_width, 0, -0.2];
            ub  = [amp/2, amp, center_gravity, edge_width, -0.0001]; #B5 values in lb ub from (D.C. Slobbe, 2016) -0.2 to -0.0001, but in some cases the tes might be positive so the B5 is positive, we assume that B5 would never be larger than the slope of the leadng edge, which usually is lower than 0.1
            bounds = (lb, ub)
           
            
            "We apply the fit over 5 iterations with increasing bounds. If still unacceptable, then flag wavelet for removal."
            wave_flag = 0
            fit_output = []
            
            for ii in range(2):
              for attempt in range(5):
                try:
                  fit_output = beta_fit(wave, betas, bounds)
                except:
                  #print(attempt)
                  it = attempt+1
                  lb  = [0,  0, (center_gravity - edge_width)-(it*2), 0, -0.2];
                  ub  = [(amp/2) + (amp/(it/10)), amp + (amp/(it/10)), center_gravity + (it*2), edge_width + (it*2), -0.0001];
                  bounds = (lb, ub)
                else:
                  break
              else:
                # we failed all the attempts - deal with the consequences. 
                wave_flag = 1


            noise_B1 = fit_output['x'][0]
            amplitude_B2 = fit_output['x'][1]
            mid_le_B3 = fit_output['x'][2]
            edge_width_B4 = fit_output['x'][3]
            tes_B5 = fit_output['x'][4]*fit_output['x'][2]
            
            snr_new = math.log10(max_p/noise_B1)*10 #[dB]
            zero_gate_new = mid_le_B3 - (edge_width_B4/2)
            peak_gate_new = mid_le_B3 + (edge_width_B4/2)
                           
            "Check if there are any wavelets that have a zero gate sooner than range bin 5, i.e. if there are any possibly cut wavelets or bad wavelets:"
            #if zero_gate < 5: print('Zero gate at < 5, wavelet ' + str(i) + ' must be rejected. The peak gate is at: ' + str(peak_gate) + '. The minimum is at: ' + str(min_p_index[0]) + '.') #we also check where the wavelet minimum is and the leading edge ends to see the shape of the wavelet
            if zero_gate <5: 
                zero_gate_flag.append(1)
            else: 
                zero_gate_flag.append(0)
            
            "Check if the peak_gate (i.e. the maxima of leading edge is the highest peak in the waveform. If not, then flag the waveform.)"
            if peak_gate_height < max(peaks_heights):
                peak_gate_flag.append(1)
            else:
                peak_gate_flag.append(0)

            satelite_alt = sat_alt[i]
            
            "Get the values of corrections for this wavelet."
            corrections_wave = [interp_corr_df.dry_tropo_corr[i], interp_corr_df.wet_tropo_corr[i], interp_corr_df.ionosph_gim_corr[i], interp_corr_df.load_corr[i], interp_corr_df.solid_corr[i], interp_corr_df.pole_corr[i]]

            "Calculate range:"
            rb = 0.4684 #This is Range bin sample [m], it is different in SIN= 0.2342 [m]
            c = 299792458 #[m/s]
            w_delay = window_delay[i]
            ref_position = len(wave)/2 #it is still 64 as in baseline B
            corr_sum = sum(corrections_wave)
            
            le_width_m = edge_width_B4*rb #leading edge width in meters
 
            "Writing the wavelet data into a list per variable."
            max_power.append(np.float32(round(max_p,4)))
            mean_power.append(np.float32(round(mean_p, 4)))
            noises.append(np.float32(round(noise_B1, 4)))
            snrs.append(np.float32(round(snr_new, 4)))
            zero_gates.append(np.float32(round(zero_gate_new, 4)))
            peak_gates.append(np.float32(round(peak_gate_new, 4)))
            le_width.append(np.float32(round(le_width_m, 4)))
            te_slope.append(np.float32(round(tes_B5, 4)))
            backscatter.append(np.float32(round(b_coeff, 4)))
            peakiness.append(np.float32(round(pp, 4)))
            max_amplitude.append(np.float32(round(amplitude_B2, 4)))
            sum_corrections.append(np.float32(round(corr_sum, 4)))
            
            
            th_power = []
            th_index = []
            r_measurement = []
            elevs = []

            for ai, a in enumerate(alpha): #where alpha is a list of alphas [0.1, 0.2, 0.3, 0.4, 0.5] 
                thp_int_index = fit_output['x'][2] - ((fit_output['x'][3]*(1-alpha[ai]))/2)
                
                if thp_int_index > center_gravity:
                    interp_thp_indx = int(thp_int_index)
                    thold_power = wave[interp_thp_indx]
                elif thp_int_index < 0:
                    interp_thp_indx = 0
                    thold_power = wave[interp_thp_indx]
                elif thp_int_index < le_index:
                    interp_thp_indx = int(thp_int_index)
                    thold_power = wave[interp_thp_indx]
                else:
                    interp_thp_indx = int(round(thp_int_index - math.floor(le_index), 1) * 10)
                    thold_power = interp_wave[interp_thp_indx]

                    
                R = 0.5*c*w_delay + rb*(thp_int_index+1-ref_position) + corr_sum
                    
                "Substract from satelite altitude:"
                elevation = satelite_alt - R
                     
                th_power.append(np.float32(round(thold_power, 4)))
                th_index.append(np.float32(round(thp_int_index, 4)))
                r_measurement.append(np.float32(round(R, 4)))
                elevs.append(np.float32(round(elevation, 4)))
            
            
            threshold_power.append(th_power)
            threshold_index.append(th_index)
            range_measurement.append(r_measurement)
            elevations.append(elevs)

        out_data = [max_power, mean_power, noises, snrs, zero_gates, 
            zero_gate_flags, peak_gates, peak_gate_flags, le_width, 
            te_slope, backscatter, peakiness, max_amplitude, sum_corrections,
            threshold_power, threshold_index, range_measurement, elevations]

        
        output_retracked_file(path_to_outfile, out_data, alpha) #save to apropriate out directory
                
                        
        print('... retracking finished.')        
            
    except (IOError, RuntimeError, ValueError, IndexError) as e:
        print(e)
    
    

