# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

from astropy import log
from astropy.io import fits
from astropy.time import Time
import numba as nb
import numpy as np
from pandas import DataFrame
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import time
from scipy.interpolate import interp1d
from scipy import stats
from scipy.optimize import curve_fit
import csv

from sofia_redux.instruments.fifi_ls.make_header import make_header
from sofia_redux.instruments.fifi_ls.lambda_calibrate import wave
from sofia_redux.toolkit.interpolate \
    import interp_1d_point_with_error as interp
from sofia_redux.toolkit.utilities \
    import (hdinsert, gethdul, write_hdul)

from sofia_redux.instruments.fifi_ls.get_atran \
    import get_atran_interpolated
from sofia_redux.instruments.fifi_ls.get_resolution \
    import get_resolution
from sofia_redux.instruments.fifi_ls.apply_static_flat \
    import get_flat, calculate_flat




__all__ = ['classify_files', 'combine_extensions', 'combine_nods', 'telluric_scaling']


def _mjd(dateobs):
    """Get the MJD from a DATE-OBS."""
    try:
        mean_time = Time(dateobs).mjd
    except (ValueError, AttributeError):
        mean_time = 0
    return mean_time

def _unix(dateobs):
    """Get the Unix Time from a DATE-OBS."""
    try:
        mean_time = Time(dateobs).unix
    except (ValueError, AttributeError):
        mean_time = 0
    return mean_time


def _read_exthdrs(hdul, key, default=0):
    """Read FIFI-LS extension headers."""
    result = []
    if len(hdul) <= 1:
        return result
    ngrating = hdul[0].header.get('NGRATING', 1)
    for idx in range(ngrating):
        name = f'FLUX_G{idx}'
        header = hdul[name].header
        result.append(header.get(key, default))
    return np.array(result)


def _from_hdul(hdul, key):
    """Read a header keyword from the PHU."""
    return hdul[0].header[key.upper().strip()]

def _em_func(e, a, c):
    return a +c*e

def em_func_b(em_spax):
    def wow(lam ,a, b, c): 
        m = a + b*lam + c*em_spax
        return m
    return wow


def _apply_flat_for_telluric(hdul, flatdata, wave, skip_err=True):
# flat data from get_flat:
    flatfile, spatdata, specdata, specwave, specerr = flatdata

    # update the header for the output file;
    # add the flat file name to it
    primehead = hdul[0].header
    # hdul.info()

    data = np.asarray(hdul[1].data, dtype=float) # Flux
    var = np.asarray(hdul[2].data,               # Stddev   
                        dtype=float) ** 2
    if data.ndim < 3:
        data = data.reshape((1, *data.shape))
        var = var.reshape((1, *var.shape))
        do_reshape = True
    else:
        do_reshape = False

    hdu_result = calculate_flat(wave, data, var, spatdata, specdata,
                                specwave, specerr, skip_err)        
    # calculate_flat return has many values, only need actual data
    return hdu_result[0][0]

def _atransmission(hdul,row,hdr0, hdul0):
    # gets transmission from current A file

    # data.shape (numramp, 16, 25) -- OTF
    # data.shape (16, 25) -- non OTF, not covered here
    numramp, numspexel, numspaxel = hdul.data.shape
    dimspexel = 1
    # Values from main header for wavelength calibration
    dichroic = hdr0['DICHROIC']
    channel=hdr0['CHANNEL']    
    obsdate = hdr0['DATE-OBS']
    # Obsdate needs special format for lambda_calibrate.py 
    try:
        obsdate = [int(x) for x in obsdate[:10].split('-')]
    except (ValueError, TypeError, IndexError):
        log.error('Invalid DATE-OBS')
        return
    channel = hdr0['CHANNEL']
    b_order = int(hdr0['G_ORD_B', -1])
    if channel == 'BLUE':
        if b_order in [1, 2]:
            blue = 'B%i' % b_order
        else:
            log.error("Invalid Blue grating order")
            return
    else:
        blue = None

    # Values from extention header for wavelength calibration
    ind = row['indpos']

    # Get spectral resolution (RESFILE keyword is added to primehead)
    resolution = get_resolution(hdr0)
    if resolution is None:
        log.error("Unable to determine spectral resolution")
        return

    # Get ATRAN data from input file or default on disk, smoothed to
    # current resolution
    atran_data = get_atran_interpolated(hdr0, resolution=resolution,
                           atran_dir=None, use_wv=True,
                           get_unsmoothed=True)
    if atran_data is None or atran_data[0] is None:
        log.error("Unable to get ATRAN data")
        return
    atran, unsmoothed = atran_data
    # Split in wavelength (x) and transmission factor (y)
    w_atran=atran[0]
    t_atran=atran[1]
    

    # Get calibration from lambda_calibrate.py
    # calibration = {'wavelength': w, 'width': p, 'wavefile': wavefile}
    # wavecal is optional to hand over, DataFrame containing wave calibration data.  May be supplied in
    # order to remove the overhead of reading the file every iteration of this function.
    
    calibration = wave(ind, obsdate, dichroic, blue=blue, wavecal=None)
    w_cal = calibration['wavelength']

    # Initialize the arrays to return for each spaxel
    t =[np.empty(numspexel) * np.nan for _ in range(numspaxel)]

    for spaxel in range(numspaxel):
        
        w_cal_loop=w_cal[:,spaxel]
        lambda_min = np.min(w_cal_loop)
        lambda_max = np.max(w_cal_loop)

        # Limit the ATRAN wavelengths to lamnda_min and lambda_max
        mask = (w_atran >= lambda_min) & (w_atran <= lambda_max)
        w_atran_masked = w_atran[mask]
        t_atran_masked = t_atran[mask]        

        # Find spexel indices where data is not NaN 
        valid_spexel = np.where(~np.isnan(hdul.data[:, :, spaxel]))[dimspexel]  # Find non-NaN indices
        if len(valid_spexel) > 0:  # Proceed only if spaxel has valid data

            # Create a function for nearest neighbor interpolation
            interp_func = interp1d(w_atran_masked, t_atran_masked, kind='nearest', fill_value='extrapolate')

            # Interpolate the y values at low-resolution x values
            t[spaxel] = interp_func(w_cal_loop[valid_spexel])

            # Bring back emission and result to original size 16 and refill NaNs at originall positions. Can be done on 
            # original data, no change in NaNs 
            # Inizialize empty array with size 16
            t_full = np.empty(hdul.data.shape[dimspexel])
            t_full[:] = np.nan
            # Initialize array with ones at non-NaN indices
            nanarray = np.zeros(hdul.data.shape[dimspexel])
            nanarray[valid_spexel] = 1
            # Create two indices for the arrays of different lengths, then only increment the index of the longer
            # array (original data). 
            original_data_index = 0
            optimized_data_index = 0

            for value in nanarray:
                if value == 1:
                    if original_data_index < hdul.data.shape[dimspexel]:
                        t_full[original_data_index] = t[spaxel][optimized_data_index]
                        original_data_index += 1
                        optimized_data_index += 1                        
                elif value == 0:
                    original_data_index += 1  # Only increment index of original data 

            # Write back into return array
            t[spaxel] = t_full
    t = np.transpose(np.array(t))

    return t

def _telluric_scaling(hdul,brow,hdr0, hdul0, sig_rel):
    # print(hdul.data)
    # print(hdul0[1].data)
    # print(np.divide(hdul.data,hdul0[1].data))
    # print(hdr0)
    # print(hdul0[0].header)


    # print(hdul0[2].header)
    # print(hdul0[2].data)
    # data.shape (16, 25), extention
    numspexel, numspaxel = hdul.data.shape
    dimspexel = 0
    dimspaxel = 1
    # hdul0.info()
    stddev = hdul0[2].data
    # print('Go',hdul0['STDDEV_G0'].data)
    # Values from main header for wavelength calibration
    dichroic = hdr0['DICHROIC']
    channel=hdr0['CHANNEL']    
    obsdate = hdr0['DATE-OBS']
    # Obsdate needs special format for lambda_calibrate.py 
    try:
        obsdate = [int(x) for x in obsdate[:10].split('-')]
    except (ValueError, TypeError, IndexError):
        log.error('Invalid DATE-OBS')
        return
    channel = hdr0['CHANNEL']
    b_order = int(hdr0['G_ORD_B', -1])
    if channel == 'BLUE':
        if b_order in [1, 2]:
            blue = 'B%i' % b_order
        else:
            log.error("Invalid Blue grating order")
            return
    else:
        blue = None

    # Values from extention header for wavelength calibration
    ind = brow['indpos']

    # Get spectral resolution (RESFILE keyword is added to primehead)
    resolution = get_resolution(hdr0)
    if resolution is None:
        log.error("Unable to determine spectral resolution")
        return

    # Get ATRAN data from input file or default on disk, smoothed to
    # current resolution
    atran_data = get_atran_interpolated(hdr0, resolution=resolution,
                           atran_dir=None, use_wv=True,
                           get_unsmoothed=True)
    if atran_data is None or atran_data[0] is None:
        log.error("Unable to get ATRAN data")
        return
    atran, unsmoothed = atran_data
    # Split in wavelength (x) and transmission factor (y)
    w_atran=atran[0]
    t_atran=atran[1]
    

    ## Get calibration from lambda_calibrate.py
    # calibration = {'wavelength': w, 'width': p, 'wavefile': wavefile}
    # wavecal is optional to hand over, DataFrame containing wave calibration data.  May be supplied in
    # order to remove the overhead of reading the file every iteration of this function.
    
    calibration = wave(ind, obsdate, dichroic, blue=blue, wavecal=None)
    w_cal = calibration['wavelength']

    ## Perform flat fielding as in apply_static_flat.py to remove noise. There, this is performed as per grating position,
    # which is already the case here. This calculation is in a for loop per grating position. 
    # The flat fielded date is only used to get the a and c values from the curve fit, un-flatted data will be used
    # in the next steps of the pipeline as before! 
    flatdata = get_flat(hdr0)
    flatval = _apply_flat_for_telluric(hdul0, flatdata, w_cal, skip_err=True)  

    
    # Bounds in the physical range for all wavelengths, c must be positive as there are no negative
    # emissions
    param_bounds = ([0, 0], [np.inf, np.inf])
    param_bounds_b = ([-np.inf,-np.inf, 0], [ np.inf, np.inf ,np.inf])

    # Initialize the arrays to return for each spaxel
    em_opt =[np.empty(numspexel) * np.nan for _ in range(numspaxel)]
    em_opt_b =[np.empty(numspexel) * np.nan for _ in range(numspaxel)]
    em_opt_sig =[np.empty(numspexel) * np.nan for _ in range(numspaxel)]
    em_opt_b_sig =[np.empty(numspexel) * np.nan for _ in range(numspaxel)]
    em_opt_sig_rel =[np.empty(numspexel) * np.nan for _ in range(numspaxel)]
    em_opt_b_sig_rel =[np.empty(numspexel) * np.nan for _ in range(numspaxel)]
    t =[np.empty(numspexel) * np.nan for _ in range(numspaxel)]
    popt_b =[np.empty(3) * np.nan for _ in range(numspaxel)]  # a, b, c
    popt =[np.empty(2) * np.nan for _ in range(numspaxel)]  # a, c
    popt_b_sig =[np.empty(3) * np.nan for _ in range(numspaxel)]  # a, b, c
    popt_sig =[np.empty(2) * np.nan for _ in range(numspaxel)]  # a, c
    popt_b_sig_rel =[np.empty(3) * np.nan for _ in range(numspaxel)]  # a, b, c
    popt_sig_rel =[np.empty(2) * np.nan for _ in range(numspaxel)]  # a, c
    pr = 0
    for spaxel in range(numspaxel):
        pr +=1
        
        w_cal_loop=w_cal[:,spaxel]
        lambda_min = np.min(w_cal_loop)
        lambda_max = np.max(w_cal_loop)
        lambda_range = (lambda_max - lambda_min)



        # Limit the ATRAN wavelengths to lamnda_min and lambda_max
        # 30 because 15 spexel intervals and half pixel added
        mask_ext = (w_atran >= (lambda_min-(lambda_range/30))) & (w_atran <= (lambda_max+(lambda_range/30)))
        mask = (w_atran >= lambda_min) & (w_atran <= lambda_max)
        w_atran_masked = w_atran[mask]
        t_atran_masked = t_atran[mask]        

        # Find spexel indices where data is not NaN 
        valid_spexel = np.where(~np.isnan(hdul.data[:, spaxel]))[dimspexel]  # Find non-NaN indices
        if len(valid_spexel) > 0:  # Proceed only if spaxel has valid data

            # Create a function for nearest neighbor interpolation
            interp_func = interp1d(w_atran_masked, t_atran_masked, kind='nearest', fill_value='extrapolate')

            # Interpolate the y values at low-resolution x values
            t[spaxel] = interp_func(w_cal_loop[valid_spexel])

            # Emission = (1 − T ransmission)*1/(lambda^5*e^(h*c/k*lambda*T)-1)
            # Simplified: Emission = 1-Transmission
            # EmissionModel = a + b ∗ λ + c ∗ Emission(λ, T)
            # em = a + b*λ + c(1-T)
            e = 1-t[spaxel]
            #  # data.shape (16, 25)
            # numspexel, numspaxel = hdul.data.shape

            # Perform optimized linear fit  
            if sig_rel:
                e_mean = np.mean(e)

                e_kehr = 1/e
                e_mean_kehr= 1/e_mean
                e_kehr_rel = e_kehr/e_mean_kehr                
                sigma_mean = np.mean(stddev[valid_spexel,spaxel])            
                sigma_rel = stddev[valid_spexel,spaxel]/sigma_mean                
                sigma_rel_used = np.sqrt(np.square(sigma_rel)+np.square(e_kehr_rel))
                sigma_stddev = stddev[valid_spexel,spaxel] # for plotting purposes
                # print('============================================')
                # print('e_mean',e_mean)
                # print('e_kehr',e_kehr)
                # print('e_kehr_rel',e_kehr_rel)
                # print('sigma_original',stddev[valid_spexel,spaxel])
                # print('sigma_mean', sigma_mean)
                # print('sigma_rel', sigma_rel)
                # print('sigma_rel', sigma_rel)
                # print('sigma_used',sigma_used)
                # print('============================================')
            else: 
                sigma_stddev = stddev[valid_spexel,spaxel]
            popt_sig[spaxel], pcov = curve_fit(_em_func,e, np.array(flatval)[valid_spexel,spaxel],
                            sigma = sigma_stddev, bounds=param_bounds)
            popt_b_sig[spaxel], pcov = curve_fit(em_func_b(e), w_cal_loop[valid_spexel], np.array(flatval)[valid_spexel,spaxel],
                            sigma = sigma_stddev, bounds=param_bounds_b) 
            
            popt_sig_rel[spaxel], pcov = curve_fit(_em_func,e, np.array(flatval)[valid_spexel,spaxel],
                            sigma = sigma_rel_used, bounds=param_bounds)
            popt_b_sig_rel[spaxel], pcov = curve_fit(em_func_b(e), w_cal_loop[valid_spexel], np.array(flatval)[valid_spexel,spaxel],
                            sigma = sigma_rel_used, bounds=param_bounds_b)
            if sig_rel:
                popt_sig[spaxel] = popt_sig_rel[spaxel]
                popt_b_sig[spaxel] = popt_b_sig_rel[spaxel]

                
                
        
            popt[spaxel], pcov = curve_fit(_em_func,e, np.array(flatval)[valid_spexel,spaxel],
                        bounds=param_bounds)
            popt_b[spaxel], pcov = curve_fit(em_func_b(e), w_cal_loop[valid_spexel], np.array(flatval)[valid_spexel,spaxel],
                                        bounds=param_bounds_b)     

            a, c = popt[spaxel]
            em_opt[spaxel] = a + c*e
            a_sig, c_sig = popt_sig[spaxel]
            a_sig_rel, c_sig_rel = popt_sig_rel[spaxel]

            em_opt_sig[spaxel] = a_sig + c_sig*e
            em_opt_sig_rel[spaxel] = a_sig_rel + c_sig_rel*e # for plotting purposes

            a_b,b_b, c_b = popt_b[spaxel]
            em_opt_b[spaxel] = a_b + b_b*w_cal_loop[valid_spexel]+c_b*e
            a_b_sig,b_b_sig, c_b_sig = popt_b_sig[spaxel]
            em_opt_b_sig[spaxel] = a_b_sig + b_b_sig*w_cal_loop[valid_spexel]+c_b_sig*e

            a_b_sig_rel,b_b_sig_rel, c_b_sig_rel = popt_b_sig_rel[spaxel]
            em_opt_b_sig_rel[spaxel] = a_b_sig_rel + b_b_sig_rel*w_cal_loop[valid_spexel]+c_b_sig_rel*e # for plotting purposes


            # Bring back emission and result to original size 16 and refill NaNs at originall positions. Can be done on 
            # original data, no change in NaNs 
            # Inizialize empty array with size 16
            restored_em_opt = np.empty(hdul.data.shape[dimspexel])
            restored_em_opt_b = np.empty(hdul.data.shape[dimspexel])
            restored_em_opt_sig = np.empty(hdul.data.shape[dimspexel])
            restored_em_opt_b_sig = np.empty(hdul.data.shape[dimspexel])
            restored_em_opt_sig_rel = np.empty(hdul.data.shape[dimspexel]) # for plotting purposes
            restored_em_opt_b_sig_rel = np.empty(hdul.data.shape[dimspexel]) # for plotting purposes
            t_full = np.empty(hdul.data.shape[dimspexel])
            restored_em_opt[:] = np.nan
            restored_em_opt_b[:] = np.nan
            restored_em_opt_sig[:] = np.nan
            restored_em_opt_b_sig[:] = np.nan
            restored_em_opt_sig_rel[:] = np.nan # for plotting purposes
            restored_em_opt_b_sig_rel[:] = np.nan # for plotting purposes
            t_full[:] = np.nan
            # Initialize array with ones at non-NaN indices
            nanarray = np.zeros(hdul.data.shape[dimspexel])
            nanarray[valid_spexel] = 1
            # Create two indices for the arrays of different lengths, then only increment the index of the longer
            # array (original data). 
            original_data_index = 0
            optimized_data_index = 0

            for value in nanarray:
                if value == 1:
                    if original_data_index < hdul.data.shape[dimspexel]:
                        restored_em_opt[original_data_index] = em_opt[spaxel][optimized_data_index]
                        restored_em_opt_b[original_data_index] = em_opt_b[spaxel][optimized_data_index]
                        restored_em_opt_sig[original_data_index] = em_opt_sig[spaxel][optimized_data_index]
                        restored_em_opt_b_sig[original_data_index] = em_opt_b_sig[spaxel][optimized_data_index]
                        restored_em_opt_sig_rel[original_data_index] = em_opt_sig_rel[spaxel][optimized_data_index] # for plotting purposes
                        restored_em_opt_b_sig_rel[original_data_index] = em_opt_b_sig_rel[spaxel][optimized_data_index] # for plotting purposes
                        t_full[original_data_index] = t[spaxel][optimized_data_index]
                        original_data_index += 1
                        optimized_data_index += 1                        
                elif value == 0:
                    original_data_index += 1  # Only increment index of original data 

            # Write back into return array
            em_opt[spaxel] = restored_em_opt
            em_opt_b[spaxel] = restored_em_opt_b
            em_opt_sig[spaxel] = restored_em_opt_sig
            em_opt_b_sig[spaxel] = restored_em_opt_b_sig
            em_opt_sig_rel[spaxel] = restored_em_opt_sig_rel # for plotting purposes
            em_opt_b_sig_rel[spaxel] = restored_em_opt_b_sig_rel # for plotting purposes
            t[spaxel] = t_full  



            # Plot some stuff for testing 
            fig, ax1 = plt.subplots(figsize=(20, 15))  # Create a single figure with one subplot
            fontsize = 50
            markersize = 8
            font = {'size'   : fontsize}

            plt.rc('font', **font)
            plt.grid(True)
            plt.grid(which='minor')
            plt.rcParams["font.family"] = "Times New Roman"

            # plt.rc('grid', color='black')
            # ax1.plot(w[valid_spexel, i], hdul.data[valid_spexel, i], label=f' 63OI Spaxel {i + 1}, Länge {len(valid_spexel)}', marker='.')
            # plt.plot(w_cal_loop[valid_spexel], hdul.data[valid_spexel, spaxel], label=f' CII Spaxel {spaxel + 1}, Länge {len(valid_spexel)}', marker='.')
            # plt.plot(w_cal_loop[valid_spexel], hdul.data[valid_spexel,spaxel], label=f' Raw Spaxel {spaxel + 1}, Länge {len(valid_spexel)}', marker='.')
            lns1 = ax1.plot(w_cal_loop[valid_spexel], np.array(flatval)[valid_spexel,spaxel],marker='s',ms = markersize, 
                            color='burlywood', label=f' Schrift {fontsize}Flat Spaxel {spaxel + 1}, {len(valid_spexel)} datapoints')
            lns2 = ax1.plot(w_cal_loop[valid_spexel], em_opt[spaxel][valid_spexel],marker='o',ms = markersize ,
                            color='m',label='a  + c(1-T): a=%5.3f, c=%5.3f' % tuple(popt[spaxel]))
            lns3 = ax1.plot(w_cal_loop[valid_spexel], em_opt_b[spaxel][valid_spexel],marker='o',ms = markersize, 
                            color='r',label='a + b*\u03BB + c(1-T): a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt_b[spaxel]))
            
            lns4 = ax1.plot(w_cal_loop[valid_spexel], em_opt_sig[spaxel][valid_spexel],marker='*',ms = markersize, 
                            color='lime',label='\u03A3, a  + c(1-T): a=%5.3f, c=%5.3f' % tuple(popt_sig[spaxel]))
            lns5 = ax1.plot(w_cal_loop[valid_spexel], em_opt_b_sig[spaxel][valid_spexel],marker='*',ms = markersize, 
                            color='g',label='\u03A3, a + b*\u03BB + c(1-T): a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt_b_sig[spaxel]))
            lns6 = ax1.plot(w_cal_loop[valid_spexel], em_opt_sig_rel[spaxel][valid_spexel],marker='D',ms = markersize, 
                            color='deepskyblue',label='$\u03A3_{normed}$, a  + c(1-T): a=%5.3f, c=%5.3f' % tuple(popt_sig_rel[spaxel]))
            lns7 = ax1.plot(w_cal_loop[valid_spexel], em_opt_b_sig_rel[spaxel][valid_spexel],marker='D',ms = markersize, 
                            color='b',label='$\u03A3_{normed}$, a + b*\u03BB + c(1-T): a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt_b_sig_rel[spaxel]))


            

            # # Create a twin Axes sharing the same x-axis
            ax2 = ax1.twinx()
            # ax2.plot(w_cal_loop[valid_spexel], t_atran, marker='.', color='red', label='ATRAN factor')
            # ax2.plot(w[valid_spexel, i], tm_spax, marker='*', color='green', label='ATRAN factor nearest')
            lns8 = ax2.plot(w_cal_loop[valid_spexel], t[spaxel][valid_spexel], marker='P',ms = markersize, 
                            color='gray',linestyle='--' , label='ATRAN factor')


            # Find the min and max value of all ax1 data for y1 limits
            minnime = np.min(np.minimum.reduce([np.array(flatval)[valid_spexel, spaxel],
                                    em_opt[spaxel][valid_spexel],
                                    em_opt_b[spaxel][valid_spexel],
                                    em_opt_sig[spaxel][valid_spexel],
                                    em_opt_b_sig[spaxel][valid_spexel],
                                    em_opt_sig_rel[spaxel][valid_spexel],
                                    em_opt_b_sig_rel[spaxel][valid_spexel]]))
            maxime = np.max(np.maximum.reduce([np.array(flatval)[valid_spexel, spaxel],
                                    em_opt[spaxel][valid_spexel],
                                    em_opt_b[spaxel][valid_spexel],
                                    em_opt_sig[spaxel][valid_spexel],
                                    em_opt_b_sig[spaxel][valid_spexel],
                                    em_opt_sig_rel[spaxel][valid_spexel],
                                    em_opt_b_sig_rel[spaxel][valid_spexel]]))
            # Round to the nearest full 100
            min_value_rounded = np.floor(minnime / 100) * 100
            max_value_rounded = np.ceil(maxime / 100) * 100            

            # Set tick locations and labels for the first y-axis
            yticks1 = np.linspace(min_value_rounded, max_value_rounded, 5)  # Example tick locations
            ax1.set_yticks(yticks1)
            ax1.set_yticklabels([f'{int(val)}' for val in yticks1])
            ax1.set_ylim(top=max_value_rounded)
            ax1.set_ylim(bottom=min_value_rounded)
            ax1.set_xlabel('Wavelength [\u03BCm]')            
            ax1.set_ylabel('Flux [I.U.]')
            ax1.tick_params(axis='y', labelcolor='black')

            # Set tick locations and labels for the second y-axis
            yticks2 = np.linspace(0, 1, 5)  # Example tick locations
            ax2.set_yticks(yticks2)
            ax2.set_yticklabels([f'{val:}' for val in yticks2])
            ax2.set_ylim(bottom=0)
            ax2.set_ylim(top=1) 
            ax2.set_ylabel('Transmission Factor [-]')
            ax2.tick_params(axis='y', labelcolor='black')

            # Customize the grid
            ax1.grid(True, linestyle='--', alpha=0.8)
            ax2.grid(True, linestyle='--', alpha=0.8)


            ax1.set_xlim(lambda_min*0.9999, lambda_max*1.0001)
            lns = lns1 + lns2 + lns3 + lns4 + lns5 + lns6 + lns7 + lns8
            labs = [l.get_label() for l in lns]
            # ax1.legend(lns, labs,bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
            #     mode="expand", borderaxespad=0, ncol=1, prop={'size':fontsize})
                        
            plt.tight_layout()
            # plt.show()

            # filename = f"63OI_Spaxel_{spaxel+1}_sigma_normed_sigma_trans.png"
            filename = f"Test_63OI_{spaxel+1}_legend_node.png"
            plt.savefig(filename)
            
    em_opt = np.transpose(np.array(em_opt))
    em_opt_b = np.transpose(np.array(em_opt_b))
    popt = np.array(popt)
    popt_b = np.array(popt_b)
    popt_sig = np.array(popt_sig)
    popt_b_sig = np.array(popt_b_sig)
    t = np.transpose(np.array(t))
    return popt, t, flatval, popt_b, em_opt_b, popt_sig, popt_b_sig, em_opt, w_cal


def classify_files(filenames, offbeam=False):
    """
    Extract various properties of all files for subsequent combination.

    Parameters
    ----------
    filenames : array_like of str
        File paths to FITS files
    offbeam : bool, optional
        If True, swap 'A' nods with 'B' nods and the following
        associated keywords: DLAM_MAP <-> DLAM_OFF,
        DBET_MAP <->  DBET_OFF.

    Returns
    -------
    pandas.DataFrame
    """
    hduls = []
    fname_list = []
    for fname in filenames:
        hdul = gethdul(fname)
        if hdul is None:
            log.error("Invalid HDUList: %s" % fname)
            continue
        hduls.append(hdul)
        if not isinstance(fname, str):
            fname_list.append(_from_hdul(hdul, 'FILENAME'))
        else:
            fname_list.append(fname)
    filenames = fname_list
    n = len(filenames)
    if n == 0:
        log.error("No good files found.")
        return None

    keywords = ['nodstyle', 'detchan' ,'channel', 'nodbeam', 'dlam_map',
                'dbet_map', 'dlam_off', 'dbet_off', 'date-obs']

    init = dict((key, [_from_hdul(hdul, key) for hdul in hduls])
                for key in keywords)
    init['mjd'] = [_mjd(dateobs) for dateobs in init['date-obs']]

    init['indpos'] = [_read_exthdrs(hdul, 'indpos', default=0)
                      for hdul in hduls]
    init['bglevl'] = [_read_exthdrs(hdul, 'bglevl_a', default=0)
                      for hdul in hduls]
    init['asymmetric'] = [x in ['ASYMMETRIC', 'C2NC2']
                          for x in init['nodstyle']]
    init['tsort'] = [0.0] * n
    init['sky'] = [False] * n  # calculate later
    init['hdul'] = hduls
    init['chdul'] = [None] * n
    init['combined'] = [np.full(len(x), False) for x in init['indpos']]
    init['outfile'] = [''] * n
    init['mstddev'] = [_read_exthdrs(hdul, 'mstddev', default=0)
                      for hdul in hduls]


    df = DataFrame(init, index=filenames)


    # If any files are asymmetric, treat them all as asymmetric
    if df['asymmetric'].any() and not df['asymmetric'].all():
        log.warning("Mismatched NODSTYLE. Will attempt to combine anyway.")

    # Drop any bad dates
    bad_dates = df[df['mjd'] == 0]
    if len(bad_dates) > 0:
        for name, row in bad_dates.iterrows():
            log.error('DATE-OBS in header is %s for %s' %
                      (row['date-obs'], name))
        df = df.drop(bad_dates.index)

    # If there's a good detchan value, use it in place of channel,
    # then set channel to either 1 (BLUE) or 0 (RED)
    valid_detchan = (df['detchan'] != 0) & (df['detchan'] != '0')
    df['channel'] = np.where(valid_detchan, df['detchan'], df['channel'])
    df['channel'] = df['channel'].apply(lambda x: 1 if x == 'BLUE' else 0)


    # update headers if offbeam is True
    if offbeam:
        # Switch A and B beams
        df['nodbeam'] = np.where(df['nodbeam'] != 'A', 'A', 'B')

        df = df.rename(index=str, columns={
            'dlam_map': 'dlam_off',
            'dbet_map': 'dbet_off',
            'dlam_off': 'dlam_map',
            'dbet_off': 'dbet_map'})

        for key in ['dlam_map', 'dlam_off', 'dbet_map', 'dbet_off', 'nodbeam']:
            df.apply(lambda x: hdinsert(
                x.hdul[0].header, key.upper(), x[key]), axis=1)

    # set on-source exptime to 0 for asym B beams and track as 'sky' files
    df['sky'] = df['asymmetric'] & (df['nodbeam'] != 'A')
    for hdul in df[df['sky']]['hdul'].values:
        if isinstance(hdul, fits.HDUList):
            hdul[0].header['EXPTIME'] = 0.0
            hdul[0].header['NEXP'] = 0

    return df


@nb.njit(cache=True, nogil=False, parallel=False, fastmath=False)
def interp_b_nods(atime, btime, bdata, berr):   # pragma: no cover
    """
    Interpolate two B nods to the A time.

    Parameters
    ----------
    atime : array-like of float
        The UNIX time for each A nod sample.
    btime : array-like float
        Before and after time for the B nods.  Expected to have two
        elements; all `atime` values should fall between the first
        and second values.
    bdata : array-like of float
        2 x nw x ns B nod data to interpolate.
    berr : array-like of float
        2 x nw x ns B nod errors to interpolate.

    Returns
    -------
    bflux, bvar : array-like of float
        nw x ns interpolated B nod flux and variance.
    """
    nt = atime.size
    nn, nw, ns = bdata.shape

    bflux = np.empty((nt, nw, ns), dtype=nb.float64)
    bvar = np.empty((nt, nw, ns), dtype=nb.float64)

    for t in range(nt):
        for i in range(nw):
            for j in range(ns):
                # flux and error at this pixel
                bf = bdata[:, i, j]
                be = berr[:, i, j]

                # Interpolate B flux and error onto A time
                if np.any(np.isnan(bf)) or np.any(np.isnan(be)):
                    bflux[t, i, j] = np.nan
                    bvar[t, i, j] = np.nan
                else:
                    # f, e = interp(btime, bf, be, atime[t])
                    f, z = interp(btime, bf, be, atime[t])
                    e, arschlecken = interp(btime, be, bf, atime[t])        #interpolate the error like the flux without squared adding
                    bflux[t, i, j] = f
                    bvar[t, i, j] = e * e

    return bflux, bvar


def combine_extensions(df, b_nod_method='nearest', bg_scaling=False, telluric_scaling_on = False):
    """
    Find a B nod for each A nod.

    For asymmetric data, DLAM and DBET do not need to match,
    B data can be used more than once, and the B needs to be
    subtracted, rather than added (symmetric B nods are
    multiplied by -1 in chop_subtract)

    For the 'interpolate' option for B nod combination for most data, the
    time of interpolation is taken to be the middle of the observation,
    as determined by the FIFISTRT and EXPTIME keywords in the primary
    header.  For OTF data, the time is interpolated between RAMPSTRT
    and RAMPEND times in the extension header, for each ramp.

    Parameters
    ----------
    df : pandas.DataFrame
    b_nod_method : {'nearest', 'average', 'interpolate'}, optional
        Determines the method of combining the two nearest before
        and after B nods.

    Returns
    -------
    list of fits.HDUList
    """
    # check B method parameter
    if b_nod_method not in ['nearest', 'average', 'interpolate']:
        raise ValueError("Bad b_nod_method: should be 'nearest', "
                         "'average', or 'interpolate'.")
    get_two = b_nod_method != 'nearest'

    df.sort_values('mjd', inplace=True)
    blist = df[df['nodbeam'] == 'B']
    alist = df[df['nodbeam'] == 'A']

    # skip if no pairs available
    if len(blist) == 0:
        log.warning('No B nods found')
        return df
    elif len(alist) == 0:
        log.error('No A nods found')
        return df
    
    for afile, arow in alist.iterrows():

            asymmetric = arow['asymmetric']
            bselect = blist[(blist['channel'] == arow['channel'])
                            & (blist['asymmetric'] == asymmetric)]

            if not asymmetric:
                bselect = bselect[(bselect['dlam_map'] == arow['dlam_map'])
                                & (bselect['dbet_map'] == arow['dbet_map'])]
            # find closest matching B image in time
            if get_two and asymmetric:
                bselect['tsort'] = bselect['mjd'] - arow['mjd']
                after = (bselect[bselect['tsort'] > 0]).sort_values('tsort')
                bselect = (bselect[bselect['tsort'] <= 0]).sort_values(
                    'tsort', ascending=False)
            else:
                bselect['tsort'] = abs(bselect['mjd'] - arow['mjd'])
                bselect = bselect.sort_values('tsort')
                after = None

            primehead, combined_hdul = None, None
            for aidx, apos in enumerate(arow['indpos']):
                bidx, bidx2 = [], []
                bfile, bfile2 = None, None
                brow, brow2 = None, None
                for bfile, brow in bselect.iterrows():
                    bidx = brow['indpos'] == apos
                    if not asymmetric:
                        bidx &= ~brow['combined']
                    if np.any(bidx):
                        break

                if after is not None:
                    for bfile2, brow2 in after.iterrows():
                        # always asymmetric
                        bidx2 = brow2['indpos'] == apos
                        if np.any(bidx2):
                            break
                    if not np.any(bidx) and np.any(bidx2):
                        bidx = bidx2
                        brow = brow2
                        bfile = bfile2
                        bidx2 = []

                describe_a = f"A {os.path.basename(arow.name)} at ext{aidx + 1} " \
                            f"channel {arow['channel']} indpos {apos} " \
                            f"dlam {arow['dlam_map']} dbet {arow['dbet_map']}"
                if np.any(bidx):
                    arow['combined'][aidx] = True
                    a_fname = f'FLUX_G{aidx}'
                    a_sname = f'STDDEV_G{aidx}'
                    a_hdr = arow['hdul'][0].header

                    bgidx = np.nonzero(bidx)[0][0]
                    brow['combined'][bgidx] = True
                    b_background = brow['bglevl'][bgidx]
                    b_fname = f'FLUX_G{bgidx}'
                    b_sname = f'STDDEV_G{bgidx}'
                    b_flux = brow['hdul'][b_fname].data
                    b_var = brow['hdul'][b_sname].data ** 2
                    b_hdr = brow['hdul'][0].header

                    # check for offbeam with OTF mode: B nods
                    # can't have an extra dimension
                    if b_flux.ndim > 2:
                        msg = 'Offbeam option is not available for OTF mode'
                        log.error(msg)
                        raise ValueError(msg)

                    combine_headers = [a_hdr, b_hdr]

                    # check for a second B nod: if not found, will do
                    # 'nearest' for this A file
                    if np.any(bidx2):
                        # add in header for combination
                        b2_hdr = brow2['hdul'][0].header
                        combine_headers.append(b2_hdr)

                        # get A and B times
                        try:
                            # read number of chop cycles per grating position from 
                            # header of current file, might change over 
                            # over observations. Might be overkill, but still correct
                            if arow['detchan'] == 'BLUE': 
                                a_chpg = a_hdr['C_CYC_B']                            
                                b_chpg1 = b_hdr['C_CYC_B']
                                b_chpg2 = b2_hdr['C_CYC_B']
                            else:
                                a_chpg =  a_hdr['C_CYC_R']
                                b_chpg1 = b_hdr['C_CYC_R']
                                b_chpg2 = b2_hdr['C_CYC_R']
                            # unix time at middle of grating position, each time looking from A file 
                            # --> 3x adix as base as current A nod is reference 

                            atime = _unix(a_hdr['DATE-OBS']) \
                                + (aidx + 0.5)*((a_hdr['C_CHOPLN']*2/250)*a_chpg)
                            btime1 = _unix(b_hdr['DATE-OBS']) \
                                + (aidx + 0.5)*((b_hdr['C_CHOPLN']*2/250)*b_chpg1)
                            btime2 = _unix(b2_hdr['DATE-OBS']) \
                                + (aidx + 0.5)*((b2_hdr['C_CHOPLN']*2/250)*b_chpg2)
                            
                        except KeyError:
                            raise ValueError('Missing DATE-OBS, C_CHOPLN or C_CYC keys in headers.')

                        # get index for second B row
                        bgidx2 = np.nonzero(bidx2)[0][0]
                        brow2['combined'][bgidx2] = True
                        b_fname = f'FLUX_G{bgidx2}'
                        b_sname = f'STDDEV_G{bgidx2}' 
                        # Perform telluric scaling 
                        med = True        # True: Takes median of factors. False: Takes curve fit params for each spaxel
                        sig = False       # True: Sigma into curve fit. False: No sigma into curve fit
                        sig_rel = True    # True: Normed Sigma. False: STDDEV
                        ac = True       # True: only a and c fit, False: a, b and c fit
                        if telluric_scaling_on:  
                            popt1, t1 , b1_flat, popt1_b, b1_fitted_b, popt1_sig, popt1_b_sig, b1_fitted, lambda1   = _telluric_scaling(brow['hdul'][b_fname],brow, brow['hdul'][0].header, brow['hdul'], sig_rel) 
                            popt2, t2, b2_flat, popt2_b, b2_fitted_b, popt2_sig, popt2_b_sig, b2_fitted, lambda2 = _telluric_scaling(brow2['hdul'][b_fname],brow2, brow2['hdul'][0].header, brow2['hdul'], sig_rel)
                            ta = _atransmission(arow['hdul'][a_fname],arow, arow['hdul'][0].header, arow['hdul'])
                            # reshape into data.shape (16, 25)
                            numspexel, numspaxel = brow['hdul'][b_fname].data.shape 

                            if ac: 
                                if sig:
                                    # Only a and c curve fit parameters  
                                    a1 =popt1_sig[:, 0]  # Extracts the first column (a values)
                                    c1= popt1_sig[:, 1]  # Extracts the 2nd column (c values if no b values)
                                    a2 =popt2_sig[:, 0]  # Extracts the first column (a values)
                                    c2= popt2_sig[:, 1]  # Extracts the 2nd column (c values if no b values)
                                else: 
                                    a1 =popt1[:, 0]  # Extracts the first column (a values)
                                    c1= popt1[:, 1]  # Extracts the 2nd column (c values if no b values)
                                    a2 =popt2[:, 0]  # Extracts the first column (a values)
                                    c2= popt2[:, 1]  # Extracts the 2nd column (c values if no b values)
                                if med:
                                    a1 = np.nanmedian(a1)  
                                    a2 = np.nanmedian(a2)
                                    c1 = np.nanmedian(c1)
                                    c2 = np.nanmedian(c2)  
                                    # write array of size 16, 25 with one value
                                    a1_full= np.full((numspexel, numspaxel), a1)
                                    c1_full= np.full((numspexel, numspaxel), c1) 
                                    a2_full= np.full((numspexel, numspaxel), a2)
                                    c2_full= np.full((numspexel, numspaxel), c2)  
                                else:                                       
                                    # Reshape into a 2D array (16, 25)
                                    a1_full= np.tile(a1, (numspexel, 1))
                                    c1_full= np.tile(c1, (numspexel, 1)) 
                                    a2_full= np.tile(a2, (numspexel, 1))
                                    c2_full= np.tile(c2, (numspexel, 1))                         
                                b1_fitted = a1_full + np.multiply(c1_full,(1-t1))
                                b2_fitted = a2_full + np.multiply(c2_full,(1-t2))
                        
                                telfac1 = 1 +  np.divide(c1_full,b1_fitted)*(t1-ta)
                                telfac2 = 1 +  np.divide(c2_full,b2_fitted)*(t2-ta)
                                b_flux = np.multiply(b_flux,telfac1)
                                b_flux2 = np.multiply(brow2['hdul'][b_fname].data,telfac2)
                                bdata = np.array([b_flux, b_flux2])
                                berr = np.array([np.multiply(np.sqrt(b_var),telfac1),
                                        np.multiply(brow2['hdul'][b_sname].data,telfac2)])
                            else:   # a ,b and c curve fit parameters                                 
                                if sig: 
                                    a1_b = popt1_b_sig[:, 0]  
                                    b1_b = popt1_b_sig[:, 1]                            
                                    c1_b = popt1_b_sig[:, 2]  
                                    a2_b = popt2_b_sig[:, 0]  
                                    b2_b = popt2_b_sig[:, 1]
                                    c2_b = popt2_b_sig[:, 2] 
                                else:                                    
                                    a1_b = popt1_b[:, 0]  
                                    b1_b = popt1_b[:, 1]                            
                                    c1_b = popt1_b[:, 2]  
                                    a2_b = popt2_b[:, 0]  
                                    b2_b = popt2_b[:, 1]
                                    c2_b = popt2_b[:, 2] 

                                if med:
                                    a1_b = np.nanmedian(a1_b)
                                    b1_b = np.nanmedian(b1_b)  
                                    c1_b = np.nanmedian(c1_b)
                                    a2_b = np.nanmedian(a2_b)
                                    b2_b = np.nanmedian(b2_b)  
                                    c2_b = np.nanmedian(c2_b) 

                                    a1_b_full= np.full((numspexel, numspaxel), a1_b)
                                    b1_b_full= np.full((numspexel, numspaxel), b1_b)
                                    c1_b_full= np.full((numspexel, numspaxel), c1_b) 
                                    a2_b_full= np.full((numspexel, numspaxel), a2_b)
                                    b2_b_full= np.full((numspexel, numspaxel), b2_b)
                                    c2_b_full= np.full((numspexel, numspaxel), c2_b)  
                                else:                       
                                    
                                    # Reshape into a 2D array (16, 25)
                                    a1_b_full= np.tile(a1_b, (numspexel, 1))
                                    b1_b_full= np.tile(b1_b, (numspexel, 1))
                                    c1_b_full= np.tile(c1_b, (numspexel, 1)) 
                                    a2_b_full= np.tile(a2_b, (numspexel, 1))
                                    b2_b_full= np.tile(b2_b, (numspexel, 1))
                                    c2_b_full= np.tile(c2_b, (numspexel, 1))

                                b1_fitted_b = a1_b_full + np.multiply(b1_b_full, lambda1) + np.multiply(c1_b_full,(1-t1))
                                b2_fitted_b  = a2_b_full + np.multiply(b2_b_full, lambda2) + np.multiply(c2_b_full,(1-t2))  

                                telfac1_b = 1 +  np.divide(c1_b_full,b1_fitted_b)*(t1-ta)
                                telfac2_b = 1 +  np.divide(c2_b_full,b2_fitted_b)*(t2-ta)

                                b_flux = np.multiply(b_flux,telfac1_b)
                                b_flux2 = np.multiply(brow2['hdul'][b_fname].data,telfac2_b)
                                bdata = np.array([b_flux, b_flux2])
                                berr = np.array([np.multiply(np.sqrt(b_var),telfac1_b),
                                        np.multiply(brow2['hdul'][b_sname].data,telfac2_b)])

                        else:  
                            bdata = np.array([b_flux, brow2['hdul'][b_fname].data])
                            berr = np.array([np.sqrt(b_var),
                                            brow2['hdul'][b_sname].data])

                        if b_nod_method == 'interpolate':
                            # debug message
                            msg = f'Interpolating B {bfile} at {btime1} ' \
                                f'and {bfile2} at {btime2} ' \
                                f'to A time {atime} and subbing from '


                            
                            # UNIX time is a range of values for OTF data:
                            # retrieve from RAMPSTRT and RAMPEND keys
                            a_hdu_hdr = arow['hdul'][a_fname].header
                            a_shape = arow['hdul'][a_fname].data.shape

                            if len(a_shape) == 3 \
                                    and 'RAMPSTRT' in a_hdu_hdr \
                                    and 'RAMPEND' in a_hdu_hdr:
                                rampstart = a_hdu_hdr['RAMPSTRT']
                                rampend = a_hdu_hdr['RAMPEND']
                                nramp = a_shape[0]
                                ramp_incr = (rampend - rampstart) / (nramp - 1)
                                atime = np.full(nramp, rampstart)
                                atime += np.arange(nramp, dtype=float) * ramp_incr
                            else:
                                atime = np.array([atime])
                            btime = np.array([btime1, btime2])                    

                            b_flux, b_var = \
                                interp_b_nods(atime, btime, bdata, berr)

                            # reshape if there was only one atime
                            if atime.size == 1:
                                b_flux = b_flux[0]
                                b_var = b_var[0]

                            # # average over possibly 4 background files
                            # # grab the first max two entries of bselect, which is already sorted to tsort
                            # closest_before_files = bselect.iloc[:2]
                            # closest_after_files = after.iloc[:2]
                            # # Condition to filter rows
                            # condition_bef = abs(closest_before_files['tsort']) < 3 / (24 * 60)
                            # condition_aft = abs(closest_after_files['tsort']) < 3 / (24 * 60)
                            # before_filtered = closest_before_files[condition_bef]
                            # after_filtered = closest_after_files[condition_aft]
                            # # only take an equal number of files (one or two of each) before and after A Nod  
                            # # to avoid a biased calculation of the background 
                            # if len(before_filtered) < 2:
                            #     after_filtered_c = after_filtered.iloc[:1]
                            #     before_filtered_c = before_filtered                            
                            # elif len(after_filtered) < 2:
                            #     before_filtered_c = before_filtered.iloc[:1]
                            #     after_filtered_c = after_filtered
                            # else:
                            #     before_filtered_c = before_filtered
                            #     after_filtered_c = after_filtered

                            # concat_df = pd.concat([before_filtered_c, after_filtered_c])                      
                        
                            # # weighted average
                            # bglevl = np.array([item[0] for item in concat_df['bglevl']])
                            # mstddev = np.array([item[0] for item in concat_df['mstddev']])
                            # weights = 1/(mstddev*mstddev)
                            # b_background = np.average(bglevl, weights=weights) 

                            
                            # # average over two background files
                            b_background += brow2['bglevl'][bgidx2]
                            b_background /= 2.

                            # # interpolate background to header atime
                            # b_background = \
                            #     np.interp(atime, [btime1, btime2],
                            #               [b_background, brow2['bglevl'][bgidx2]])
                            # average background
                        else:
                            # debug message
                            msg = f'Averaging B {bfile} and {bfile2} ' \
                                f'and subbing from '

                            # average flux                           
                            b_flux += b_flux2
                            b_flux /= 2.

                            # propagate variance
                            b_var += brow2['hdul'][b_sname].data ** 2
                            b_var /= 4.

                            # average background
                            b_background += brow2['bglevl'][bgidx2]
                            b_background /= 2.

                    else:
                        if asymmetric:
                            msg = f'Subbing B {os.path.basename(brow.name)} from '
                        else:
                            msg = f'Adding B {os.path.basename(brow.name)} to '

                    log.debug(msg + describe_a)

                    # Note: in the OTF case, A data is a 3D cube with
                    # ramps x spexels x spaxels, and B data is a
                    # 2D array of spexels x spaxels.  The B data is
                    # subtracted at each ramp.
                    # For other modes, A and B are both spexels x spaxels.

                    flux = arow['hdul'][a_fname].data
                    stddev = arow['hdul'][a_sname].data ** 2 + b_var             
                                
                    if asymmetric:
                        # Optional background scaling for unchopped observations
                        if bg_scaling and a_hdr['C_AMP']==0:
                            a_background = arow['bglevl'][aidx]
                            flux -= b_flux*a_background/b_background
                        else:
                            flux -= b_flux
                    else:
                        # b_flux from source is negative for symmetric chops
                        # as result of subtract chops
                        flux += b_flux
                        # divide by two for doubled source
                        flux /= 2
                        stddev /= 4
                    stddev = np.sqrt(stddev)

                    if combined_hdul is None:
                        primehead = make_header(combine_headers)
                        primehead['HISTORY'] = 'Nods combined'
                        hdinsert(primehead, 'PRODTYPE', 'nod_combined')
                        outfile, _ = os.path.splitext(os.path.basename(afile))
                        outfile = '_'.join(outfile.split('_')[:-2])
                        outfile += '_NCM_%s.fits' % primehead.get('FILENUM')
                        df.loc[afile, 'outfile'] = outfile
                        hdinsert(primehead, 'FILENAME', outfile)
                        combined_hdul = fits.HDUList(
                            fits.PrimaryHDU(header=primehead))

                    exthead = arow['hdul'][a_fname].header
                    hdinsert(exthead, 'BGLEVL_B', b_background,
                            comment='BG level nod B (ADU/s)')
                    combined_hdul.append(fits.ImageHDU(flux, header=exthead,
                                                    name=a_fname))
                    combined_hdul.append(fits.ImageHDU(stddev, header=exthead,
                                                    name=a_sname))

                    # add in scanpos table from A nod if present
                    a_pname = f'SCANPOS_G{aidx}'
                    if a_pname in arow['hdul']:
                        combined_hdul.append(arow['hdul'][a_pname].copy())
                else:
                    msg = "No matching B found for "
                    log.debug(msg + describe_a)

            if combined_hdul is not None:
                df.at[afile, 'chdul'] = combined_hdul

    return df


def combine_nods(filenames, offbeam=False, b_nod_method='nearest',
                 outdir=None, write=False, bg_scaling = False, telluric_scaling_on = False):
    """
    Combine nods of ramp-fitted, chop-subtracted data.

    Writes a single FITS file to disk for each A nod found.  Each
    HDU list contains n_san binary table extensions, each containing
    DATA and STDDEV data cubes, each 5x5x18.  The output filename is
    created from the input filename, with the suffix 'CSB', 'RP0' or
    'RP1' replaced with 'NCM', and with input file numbers numbers
    concatenated.  Unless specified, the output directory is the same
    as the input files.

    Input files should have been generated by `subtract_chops`, or
    `fit_ramps` (for total power mode, which has no chops).

    The procedure is:

        1. Read header information from each extension in each of the
           input files, making lists of A data and B data, with relevant
           metadata (dither) position, date/time observed (DATE-OBS),
           inductosyn position, channel, nod style).

        2. Loop though all A data to find matching B data

            a. asymmetric nod style: find closest B nod in time with the
            same channel and inductosyn position.  Dither position does
            not have to match, B data can be used more than once, and
            data must be subtracted rather than added.

            b. symmetric nod style: find closest B nod in time with the
            same channel, inductosyn position, and dither position. Each
            B nod can only be used once, since it contains a source
            observation, and data must be added rather than subtracted.

        3. After addition or subtraction, create a FITS file and write
        results to disk.


    Parameters
    ----------
    filenames : array_like of str
        File paths to the data to be combined
    offbeam : bool, optional
        If True, swap 'A' nods with 'B' nods and the following
        associated keywords: DLAM_MAP <-> DLAM_OFF,
        DBET_MAP <->  DBET_OFF. This option cannot be used with
        OTF-mode A nods.
    b_nod_method : {'nearest', 'average', 'interpolate'}, optional
        For asymmetric, data this option controls how the nearest B nods
        are combined. The 'nearest' option takes only the nearest B nod
        in time.  The 'average' option averages the nearest before and
        after B nods.  The 'interpolate' option linearly interpolates the
        nearest before and after B nods to the time of the A data.
    outdir : str, optional
        Directory path to write output.  If None, output files
        will be written to the same directory as the input files.
    write : bool, optional
        If True, write to disk and return the path to the output
        file.  If False, return the HDUL.

    Returns
    -------
    pandas.DataFrame
        The output pandas dataframe contains a huge variety of
        information indexed by original filename.  The combined
        A-B FITS data are located in the 'chdul' column.  Note that
        only A nod files contain combined data in this 'chdul'
        column.  For example, in order to extract combined FITS
        data, one could issue the command::

            df = combine_nods(filenames)
            combined_hduls = df[df['nodbeam'] == 'A']['chdul']

        In order to extract rows from the dataframe that were not
        combined issue the command::

            not_combined = df[(df['nodbeam'] == 'A') & (df['chdul'] == None)]

        files are considered 'combined' if at least one A extension was
        combined for an A-nod file.  A true signifier of whether an
        extension was combined (both A and B nod files) can be found in the
        'combined' column as a list of bools, one for each extension.
    """
    if isinstance(filenames, str):
        filenames = [filenames]
    if not hasattr(filenames, '__len__'):
        log.error("Invalid input files type (%s)" % repr(filenames))
        return

    if isinstance(outdir, str):
        if not os.path.isdir(outdir):
            log.error("Output directory %s does not exist" % outdir)
            return
    df = classify_files(filenames, offbeam=offbeam)
    if df is None:
        log.error("Problem in file classification")
        return
    

    df = combine_extensions(df, b_nod_method=b_nod_method, bg_scaling=bg_scaling,  telluric_scaling_on = telluric_scaling_on)

    for filename, row in df[df['nodbeam'] == 'A'].iterrows():

        if outdir is not None:
            outdir = str(outdir)
        else:
            outdir = os.path.dirname(filename)

        if write and row['chdul'] is not None:
            write_hdul(row['chdul'], outdir=outdir, overwrite=True)
        if row['outfile'] is not None:
            df.at[filename, 'outfile'] = os.path.join(
                outdir, os.path.basename(row['outfile']))

    return df
