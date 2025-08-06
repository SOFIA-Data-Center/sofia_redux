# Licensed under a 3-clause BSD style license - see LICENSE.rst

import glob
import os
import re
import warnings

from astropy import log
from astropy.io import fits
import numpy as np

from sofia_redux.instruments import fifi_ls
from sofia_redux.instruments.fifi_ls.get_resolution import get_resolution
from sofia_redux.toolkit.utilities import goodfile, gethdul, hdinsert
from sofia_redux.spectroscopy.smoothres import smoothres


__all__ = ['clear_atran_cache', 'get_atran_from_cache',
           'store_atran_in_cache', 'get_atran', 'get_atran_interpolated']

__atran_cache = {}


def clear_atran_cache():
    """
    Clear all data from the atran cache.
    """
    global __atran_cache
    __atran_cache = {}


def get_atran_from_cache(atranfile, resolution):
    """
    Retrieves atmospheric transmission data from the atran cache.

    Checks to see if the file still exists, can be read, and has not
    been altered since last time.  If the file has changed, it will be
    deleted from the cache so that it can be re-read.

    Parameters
    ----------
    atranfile : str
        File path to the atran file
    resolution : float
        Spectral resolution used for smoothing.

    Returns
    -------
    tuple
        filename : str
            Used to update ATRNFILE in FITS headers.
        wave : numpy.ndarray
            (nwave,) array of wavelengths.
        unsmoothed : numpy.ndarray
            (nwave,) array containing the atran data from file
        smoothed : numpy.ndarray
            (nwave,) array containing the smoothed atran data
    """
    global __atran_cache

    key = atranfile, int(resolution)

    if key not in __atran_cache:
        return

    if not goodfile(atranfile):
        try:
            del __atran_cache[key]
        except KeyError:   # pragma: no cover
            # could happen in race conditions, working in parallel
            pass
        return

    modtime = str(os.path.getmtime(atranfile))
    if modtime not in __atran_cache.get(key, {}):
        return

    log.debug(f'Retrieving ATRAN data from cache '
              f'({key[0]}, resolution {key[1]})')
    return __atran_cache.get(key, {}).get(modtime)


def store_atran_in_cache(atranfile, resolution, filename, wave,
                         unsmoothed, smoothed):
    """
    Store atran data in the atran cache.

    Parameters
    ----------
    atranfile : str
        File path to the atran file
    resolution : float
        Spectral resolution used for smoothing.
    filename : str
        Used to update ATRNFILE in FITS headers.
    wave : numpy.ndarray
        (nwave,) array of wavelengths.
    unsmoothed : numpy.ndarray
        (nwave,) array containing the atran data from file
    smoothed : numpy.ndarray
        (nwave,) array containing the smoothed atran data

    """
    global __atran_cache
    key = atranfile, int(resolution)
    log.debug(f'Storing ATRAN data in cache '
              f'({key[0]}, resolution {key[1]})')
    __atran_cache[key] = {}
    modtime = str(os.path.getmtime(atranfile))
    __atran_cache[key][modtime] = (
        filename, wave, unsmoothed, smoothed)


def get_atran(header, resolution=None, filename=None,
              get_unsmoothed=False, use_wv=False,
              atran_dir=None):
    """
    Retrieve reference atmospheric transmission data.

    ATRAN files in the data/atran_files directory should be named
    according to the altitude, ZA, and wavelengths for which they
    were generated, as:

        atran_[alt]K_[za]deg_[wmin]-[wmax]mum.fits

    For example, the file generated for altitude of 41,000 feet,
    ZA of 45 degrees, and wavelengths between 40 and 300 microns
    should be named:

        atran_41K_45deg_40-300mum.fits

    If use_wv is set, files named for the precipitable water vapor
    values for which they were generated will be used instead, as:

        atran_[alt]K_[za]deg_[wv]pwv_[wmin]-[wmax]mum.fits

    The procedure is:

        1. Identify ATRAN file by ZA, Altitude, and WV , unless override is
           provided.
        2. Read ATRAN data from file and smooth to expected spectral
           resolution.
        3. Return transmission array.

    Parameters
    ----------
    header : fits.Header
        ATRNFILE keyword is written to the provided FITS header,
        containing the name of the ATRAN file used.
    resolution : float, optional
        Spectral resolution to which ATRAN data should be smoothed.
    filename : str, optional
        Atmospheric transmission file to be used.  If not provided,
        a default file will be retrieved from the data/atran_files
        directory.  The file with the closest matching ZA and
        Altitude to the input data will be used.  If override file
        is provided, it should be a FITS image file containing
        wavelength and transmission data without smoothing.
    get_unsmoothed : bool, optional
        If True, return the unsmoothed atran data in addition to the
        smoothed.
    atran_dir : str, optional
        Path to a directory containing ATRAN reference FITS files.
        If not provided, the default set of files packaged with the
        pipeline will be used.
    use_wv : bool, optional
        If set, water vapor values from the header will be used
        to select the correct ATRAN file.

    Returns
    -------
    numpy.ndarray
        Filename of the atmospheric transmission file used and a
        (2, nw) array containing wavelengths and (optionally
        smoothed) transmission data.
    """
    if filename is not None:
        if not goodfile(filename, verbose=True):
            log.warning(f'File {filename} not found; '
                        f'retrieving default')
            filename = None
    if not isinstance(header, fits.Header):
        log.error('Invalid header')
        return
    if resolution is None:
        log.warning('Getting default resolution from G_WAVE in header')
        resolution = get_resolution(header)

    if filename is None:
        za_start = float(header.get('ZA_START', 0))
        za_end = float(header.get('ZA_END', 0))
        if za_start > 0 >= za_end:
            za = za_start
        elif za_end > 0 >= za_start:
            za = za_end
        else:
            za = 0.5 * (za_start + za_end)

        alt_start = float(header.get('ALTI_STA', 0))
        alt_end = float(header.get('ALTI_END', 0))
        if alt_start > 0 >= alt_end:
            alt = alt_start
        elif alt_end > 0 >= alt_start:
            alt = alt_end
        else:
            alt = 0.5 * (alt_start + alt_end)
        alt /= 1000

        wv_obs = float(header.get('WVZ_OBS', 0))
        if wv_obs > 0:
            wv = wv_obs
        else:
            wv_start = float(header.get('WVZ_STA', 0))
            wv_end = float(header.get('WVZ_END', 0))
            if wv_start > 0 >= wv_end:
                wv = wv_start
            elif wv_end > 0 >= wv_start:
                wv = wv_end
            else:
                wv = 0.5 * (wv_start + wv_end)
        if use_wv and wv < 2:
            log.warning(f'Bad WV value: {wv}')
            log.warning('Using default ATRAN file.')
            use_wv = False

        log.debug(f'Alt, ZA, WV: {alt:.2f} {za:.2f} {wv:.2f}')

        true_value = [alt, za, wv]

        if atran_dir is not None:
            if not os.path.isdir(str(atran_dir)):
                log.warning(f'Cannot find ATRAN directory: {atran_dir}')
                log.warning('Using default ATRAN set.')
                atran_dir = None
        if atran_dir is None:
            atran_dir = os.path.join(os.path.dirname(fifi_ls.__file__),
                                     'data', 'atran_files')

        atran_files = glob.glob(os.path.join(atran_dir, 'atran*fits'))
        regex1 = re.compile(r'^atran_([0-9]+)K_([0-9]+)deg_40-300mum\.fits$')
        regex2 = re.compile(r'^atran_([0-9]+)K_([0-9]+)deg_'
                            r'([0-9]+)pwv_40-300mum\.fits$')

        # set up some values for tracking best atran match
        wv_overall_val = np.inf
        wv_best_file = None
        overall_val = np.inf
        best_file = None
        for f in atran_files:
            # check for WV match
            match = regex2.match(os.path.basename(f))
            if use_wv and wv > 0 and match is not None:
                match_val = 0
                for i in range(3):
                    # file alt, za, or wv
                    file_val = float(match.group(i + 1))
                    # check difference from true value
                    d_val = abs(file_val - true_value[i]) / true_value[i]
                    match_val += d_val
                if match_val < wv_overall_val:
                    wv_overall_val = match_val
                    wv_best_file = f
            else:
                # otherwise, check for non-WV match
                match = regex1.match(os.path.basename(f))
                if match is not None:
                    match_val = 0
                    for i in range(2):
                        # file alt or za
                        file_val = float(match.group(i + 1))
                        # check difference from true value
                        d_val = abs(file_val - true_value[i]) / true_value[i]
                        match_val += d_val
                    if match_val < overall_val:
                        overall_val = match_val
                        best_file = f

        if use_wv and wv_best_file is not None:
            log.debug('Using nearest Alt/ZA/WV')
            filename = wv_best_file
        else:
            log.debug('Using nearest Alt/ZA')
            filename = best_file

    if filename is None:
        log.error('No ATRAN file found')
        return

    # Read the atran data from cache if possible
    log.debug(f'Using ATRAN file: {filename}')
    atrandata = get_atran_from_cache(filename, resolution)
    if atrandata is not None:
        atranfile, wave, unsmoothed, smoothed = atrandata
    else:
        hdul = gethdul(filename, verbose=True)
        if hdul is None or hdul[0].data is None:
            log.error(f'Invalid data in ATRAN file {filename}')
            return
        data = hdul[0].data

        atranfile = os.path.basename(filename)
        wave = data[0]
        unsmoothed = data[1]
        smoothed = smoothres(data[0], data[1], resolution)

        store_atran_in_cache(filename, resolution, atranfile,
                             data[0], data[1], smoothed)

    hdinsert(header, 'ATRNFILE', atranfile)
    if not get_unsmoothed:
        return np.vstack((wave, smoothed))
    else:
        return (np.vstack((wave, smoothed)),
                np.vstack((wave, unsmoothed)))


def get_atran_interpolated(header, resolution=None,
                          get_unsmoothed=False, use_wv=False,
                          atran_dir=None):

    if not isinstance(header, fits.Header):
        log.error('Invalid header')
        return

    if resolution is None:
        log.warning('Getting default resolution from G_WAVE in header')
        resolution = get_resolution(header)

    # get ZA
    za_start = float(header.get('ZA_START', 0))
    za_end = float(header.get('ZA_END', 0))
    if za_start > 0 >= za_end:
        za = za_start
    elif za_end > 0 >= za_start:
        za = za_end
    else:
        za = 0.5 * (za_start + za_end)

    # get altitude in thousands feet
    alt_start = float(header.get('ALTI_STA', 0))
    alt_end = float(header.get('ALTI_END', 0))
    if alt_start > 0 >= alt_end:
        alt = alt_start
    elif alt_end > 0 >= alt_start:
        alt = alt_end
    else:
        alt = 0.5 * (alt_start + alt_end)
    alt /= 1000

    # get water vapor
    wv_obs = float(header.get('WVZ_OBS', 0))
    if wv_obs > 0:
        wv = wv_obs
    else:
        wv_start = float(header.get('WVZ_STA', 0))
        wv_end = float(header.get('WVZ_END', 0))
        if wv_start > 0 >= wv_end:
            wv = wv_start
        elif wv_end > 0 >= wv_start:
            wv = wv_end
        else:
            wv = 0.5 * (wv_start + wv_end)

    if use_wv and wv < 2:
        log.warning(f'Bad WV value: {wv}')
        log.warning('Using default ATRAN file.')
        use_wv = False

    log.debug(f'Alt, ZA, WV: {alt:.2f} {za:.2f} {wv:.2f}')

    if atran_dir is not None:
        if not os.path.isdir(str(atran_dir)):
            log.warning(f'Cannot find ATRAN directory: {atran_dir}')
            log.warning('Using default ATRAN set.')
            atran_dir = None
    if atran_dir is None:
        atran_dir = os.path.join(os.path.dirname(fifi_ls.__file__),
                                    'data', 'atran_files')

    za_values = list(range(30,75,5))
    wv_values = [1, 2, 3, 4, 5, 6, 7,
                8, 9, 10, 11, 12, 13,
                14, 15, 16, 17, 18, 20,
                22, 25, 27, 30, 32, 35,
                37, 40, 45, 50]
    
    # clip values to atran data range
    if not za_values[-1] >= za >= za_values[0]:
        log.warning('za={} outside of available ATRAN data.'.format(za))
        za = np.clip(za, a_min=za_values[0], a_max=za_values[-1])
        log.warning('Setting zenith angle to {} deg'.format(za))

    if use_wv and not wv_values[-1] >= wv >= wv_values[0]:
        log.warning('wv={} outside of available ATRAN data.'.format(wv))
        wv = np.clip(wv, a_min=wv_values[0], a_max=wv_values[-1])
        log.warning('Setting water vaper to {} um'.format(wv))

    if not 45 >= alt >= 38:
        log.warning('alt={} outside of available ATRAN data.'.format(alt))
        alt = np.clip(alt, a_min=38, a_max=45)
        log.warning('Setting altitude to {}K ft'.format(alt))

    za_high, za_low = np.inf, np.inf
    wv_high, wv_low = np.inf, np.inf

    # find the higher and lower boundaries
    for i, w in enumerate(wv_values):
        if w >= wv:
            wv_high = w
            wv_low = wv_values[i-1]
            break

    for i, z in enumerate(za_values):
        if z >= za:
            za_high = z
            za_low = za_values[i-1]
            break

    # build filenames
    current_atran_filenames = {}
    if use_wv:
        current_atran_filenames["za1_wv1"] = ("atran_{}K_{}deg_{}pwv_40-300mum.fits".format(round(alt), za_low, wv_low), za_low, wv_low)
        current_atran_filenames["za1_wv2"] = ("atran_{}K_{}deg_{}pwv_40-300mum.fits".format(round(alt), za_low, wv_high), za_low, wv_high)
        current_atran_filenames["za2_wv1"] = ("atran_{}K_{}deg_{}pwv_40-300mum.fits".format(round(alt), za_high, wv_low), za_high, wv_low)
        current_atran_filenames["za2_wv2"] = ("atran_{}K_{}deg_{}pwv_40-300mum.fits".format(round(alt), za_high, wv_high), za_high, wv_high)
    else:
        current_atran_filenames["za1"] = ("atran_{}K_{}deg_40-300mum.fits".format(round(alt), za_low,), za_low, None)
        current_atran_filenames["za2"] = ("atran_{}K_{}deg_40-300mum.fits".format(round(alt), za_high), za_high, None)

    # load all files
    atran_data = {}
    atrnfile_fits_keyword = ''
    for key, value in current_atran_filenames.items():
        filename = value[0]
        _za = value[1]  # the ZA value of this file
        _wv = value[2]  # the WV value of this file

        atrandata = get_atran_from_cache(filename, resolution)

        if atrandata is not None:
            atran_data[key] = atrandata, _za, _wv
        else:
            hdul = gethdul(os.path.join(atran_dir, filename), verbose=True)
            if hdul is None or hdul[0].data is None:
                log.error(f'Invalid data in ATRAN file {filename}')
                return
            data = hdul[0].data

            atranfile = os.path.basename(filename)
            wave = data[0]
            unsmoothed = data[1]
            smoothed = smoothres(data[0], data[1], resolution)

            atran_data[key] = (atranfile, wave, unsmoothed, smoothed), _za, _wv

            atrnfile_fits_keyword += filename + ', '
            store_atran_in_cache(os.path.join(atran_dir, filename), resolution, atranfile,
                                data[0], data[1], smoothed)

    if use_wv:
        # interpolate za for two pwv
        za1_za2_wv1 = interpolate_two_atran_files(atran_data["za1_wv1"],
                                                  atran_data["za2_wv1"],
                                                  za, 0)
        za1_za2_wv2 = interpolate_two_atran_files(atran_data["za1_wv2"],
                                                  atran_data["za2_wv2"],
                                                  za, 0)
        # interpolate WV and hold ZA
        interpolated = interpolate_two_atran_files(za1_za2_wv1,
                                                   za1_za2_wv2,
                                                   wv, 1)
    else:
        interpolated = interpolate_two_atran_files(atran_data["za1"],
                                                   atran_data["za2"],
                                                   za, 0)

    wave = interpolated[0][1]
    unsmoothed = interpolated[0][2]
    smoothed = smoothres(wave, unsmoothed, resolution)

    hdinsert(header, 'ATRNFILE', atrnfile_fits_keyword[:-2])

    if not get_unsmoothed:
        return np.vstack((wave, smoothed))
    else:
        return (np.vstack((wave, smoothed)),
                np.vstack((wave, unsmoothed)))


def interpolate_two_atran_files(atran_data1, atran_data2, des_val, itype=0):
    """
    Interpolate data from two atran files.

    Parameters
    ----------
    atran_data1, atran_data2 : tuple
        ((filename, wave, unsmoothed, smoothed), za, wv)
    des_val : float
        The value to be interpolated for. Can be zenith angle (ZA) or water
        vapor (WV)
    itype: int
        0 for ZA
        1 for WV

    Returns
    -------
    interpolated_atran_data: tuple
        same nested tuple format as input parameters
        (("", wave, unsmoothed, None), za, wv)
    """
    if not np.allclose(atran_data1[0][1], atran_data2[0][1]):
        raise ValueError("Incompatible atran wavelengths to interpolate")

    dy = atran_data2[0][2] - atran_data1[0][2]

    if itype == 0:
        # interpolating zenith angle
        dx = atran_data2[1] - atran_data1[1]
        t = dy/dx
        data_int = atran_data1[0][2] + t * (des_val - atran_data1[1])
        za = des_val
        wv = atran_data1[2]
    elif itype == 1:
        # interpolating water vapor
        dx = atran_data2[2] - atran_data1[2]
        t = dy/dx
        data_int = atran_data1[0][2] + t * (des_val - atran_data1[2])
        za = atran_data1[1]
        wv = des_val
    else:
        raise TypeError("itype has to be either 0 (ZA) or 1 (WV)")

    interpolated_atran_data0 = ("", atran_data1[0][1], data_int, None)
    return (interpolated_atran_data0, za, wv)
