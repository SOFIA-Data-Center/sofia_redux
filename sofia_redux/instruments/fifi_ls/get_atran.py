# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
from pathlib import Path

from astropy import log
from astropy.io import fits
from astropy.time import Time
import numpy as np

from sofia_redux.instruments import fifi_ls
from sofia_redux.instruments.fifi_ls.get_resolution import get_resolution
from sofia_redux.toolkit.utilities import goodfile, gethdul, hdinsert
from sofia_redux.toolkit.utilities.darus import get_file_from_darus
from sofia_redux.spectroscopy.smoothres import smoothres


__all__ = ['clear_atran_cache', 'get_atran_from_cache',
           'store_atran_in_cache', 'get_atran',
           'clear_ecmwf_cache']

__atran_cache = {}
__ecmwf_cache = {}


# doi per altitude
ATRAN_DATASETS = {
    35: "10.18419/DARUS-5705",
    36: "10.18419/DARUS-5707",
    37: "10.18419/DARUS-5708",
    38: "10.18419/DARUS-5709",
    39: "10.18419/DARUS-5710",
    40: "10.18419/DARUS-5711",
    41: "10.18419/DARUS-5712",
    42: "10.18419/DARUS-5713",
    43: "10.18419/DARUS-5714",
    44: "10.18419/DARUS-5715",
    45: "10.18419/DARUS-5716",
}

PWV_DOI = "10.18419/DARUS-5728"

ATRAN_ZA_VALUES = list(range(30,75,5))
ATRAN_WV_VALUES = [
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
    11, 12, 13, 14, 15, 16, 17, 18,
    20, 22, 25,     30, 32, 35, 37,
    40, 45, 50,
]

def clear_atran_cache():
    """
    Clear all data from the atran cache.
    """
    global __atran_cache
    __atran_cache = {}


def clear_ecmwf_cache():
    """
    Clear all data from the ECMWF cache.
    """
    global __ecmwf_cache
    __ecmwf_cache = {}


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


def get_atran_data(filename, resolution, atran_dir=None):
    """Get ATRAN data for given filename.

    Lookup Order:
        1. Memory cache if already read in
        2. ATRAN directory if provided
        3. DaRUS, SOFIA Astronomy Dataverse

    Parameters
    ----------
    filename : str
        Name of the ATRAN file to be read, also used for cache key
    resolution : float
        Spectral resolution to which ATRAN data should be smoothed.
    atran_dir : str
        Path to a directory containing ATRAN reference FITS files,
        comes from the recipe configuration.

    Returns
    -------
    atranfile : str
        The basename of the ATRAN file
    wave : numpy.ndarray
        (nwave,) array of wavelengths
    unsmoothed : numpy.ndarray
        (nwave,) array containing the transmission spectrum from file
    smoothed : numpy.ndarray
        (nwave,) array containing the smoothed transmission spectrum
    """
    atran_data = get_atran_from_cache(filename, resolution)

    if atran_data is not None:
        return atran_data

    if atran_dir is not None:
        if not os.path.isdir(str(atran_dir)):
            log.warning(f'Cannot find ATRAN directory: {atran_dir}')
            log.warning('Using default ATRAN set.')
            atran_dir = None
    if atran_dir is None:
        atran_dir = os.path.join(os.path.dirname(fifi_ls.__file__),
                                'data', 'atran_files')

    if not filename.startswith("atran_sdc_"):
        log.info('Non-standard ATRAN filename not starting with "atran_sdc_"'
                 f'Looking for {atran_dir}/{filename} only.')
        localpath = os.path.join(atran_dir, filename)
    else:
        parts = filename.split('_')
        if parts[2][-1] != 'K':
            raise ValueError(f'Invalid ATRAN filename: {filename}')
        alt = int(parts[2][:-1])
        localpath = os.path.join(atran_dir, f'{alt}K', filename)
        if not goodfile(localpath):
            log.debug(f'ATRAN file not found in ATRAN directory: {localpath}')
            localpath = get_atran_from_darus(alt, filename)

    atranfile = os.path.basename(filename)
    hdul = gethdul(localpath, verbose=True)
    if hdul is None:
        log.error(f'Invalid data in ATRAN file {localpath}')
        return

    if len(hdul) < 2 or hdul[1].header['CONTENT'] != "ATRAN_SDC Model":
        log.warning("Did not find ATRAN SDC model data in FITS BinTableHDU "
                    f"in {localpath}. Falling back to old format.")
        if hdul[0].data is None:
            log.error(f'Invalid data in ATRAN file {localpath}')
            return
        wave = hdul[0].data[0, :]
        unsmoothed = hdul[0].data[1, :]
    else:
        data_rec = hdul[1].data
        wave = data_rec['wavelength']
        unsmoothed = data_rec['transmission']

    smoothed = smoothres(wave, unsmoothed, resolution)
    store_atran_in_cache(localpath, resolution,
                         atranfile, wave, unsmoothed, smoothed)
    return (atranfile, wave, unsmoothed, smoothed)


def get_atran_from_darus(altitude, filename):
    """Download ATRAN file from DaRUS, SOFIA Astronomy Dataverse.

    Parameters
    ----------
    altitude : int
        Altitude in kft
    filename : str
        Name of the ATRAN file to be downloaded

    Returns
    -------
    str
        Local file path to the downloaded ATRAN file
    """
    dataset_doi = ATRAN_DATASETS.get(altitude)
    if dataset_doi is None:
        raise ValueError(f'No dataset DOI found for altitude {altitude}K')
    local_file = get_file_from_darus(dataset_doi, filename)
    log.info(f'ATRAN file in astropy cache: {local_file}')
    return local_file


def get_ecmwf_from_darus(filename):
    """Download ECMWF PWV file from DaRUS, SOFIA Astronomy Dataverse.

    Parameters
    ----------
    filename : str
        Name of the ECMWF PWV file to be downloaded,
        e.g. '855_FI_20220115_pwv.fits'.

    Returns
    -------
    str
        Local file path to the downloaded ECMWF file.
    """
    local_file = get_file_from_darus(PWV_DOI, filename)
    log.info(f'ECMWF file in astropy cache: {local_file}')
    return local_file


def get_ecmwf_from_cache(ecmwf_file):
    """Return cached ECMWF arrays for ecmwf_file, or None if not cached."""
    global __ecmwf_cache
    return __ecmwf_cache.get(str(ecmwf_file))


def store_ecmwf_in_cache(ecmwf_file, unixtime, pwv, posconst, posjump):
    """Store ECMWF arrays in cache keyed by file path."""
    global __ecmwf_cache
    __ecmwf_cache[str(ecmwf_file)] = (unixtime, pwv, posconst, posjump)


def get_wv_from_ecmwf(header, ecmwf_dir=None):
    """
    Get water vapor value from ECMWF reanalysis data.

    ECMWF files are stored as binary table FITS files named
    <flight>_<instrument>_<YYYYMMDD>_pwv.fits. This function
    finds the right file and extracts the water vapor value closest
    to the observation time.

    Lookup order:
        1. Local directory ecmwf_dir if set.
        2. DaRUS, SOFIA Astronomy Dataverse (PWV_DOI).

    Parameters
    ----------
    header : fits.Header
        FITS header containing MISSN-ID and DATE-OBS keywords.
    ecmwf_dir : str, optional
        Path to a local directory containing ECMWF FITS files.
        If not provided or the file is not found there, the file
        will be downloaded from DaRUS.

    Returns
    -------
    tuple or None
        Tuple of (WVZECMWF, wvz_fifi, wv_formula, filename) or None
        if ECMWF data could not be retrieved.
    """
    # Construct mission ID from header
    missn_id = header.get('MISSN-ID', '')
    if not missn_id:
        log.warning('MISSN-ID not found in header')
        return None

    # Format: last 3 chars + chars 10-14 + first 10 chars (no dashes)
    # Example: "2022-01-15_FI_F855" -> "855_FI_20220115"
    try:
        ecmwf_mission_id = (missn_id[-3:] + missn_id[10:14]
                           + missn_id[:10].replace('-', ''))
    except (IndexError, TypeError):
        log.warning(f'Could not get MISSN-ID: {missn_id}')
        return None

    # Get observation time
    date_obs = header.get('DATE-OBS', '')
    if not date_obs:
        log.warning('DATE-OBS not found in header')
        return None

    try:
        obs_time_unix = Time(date_obs).unix
    except Exception as e:
        log.warning(f'Could not get DATE-OBS: {date_obs} ({e})')
        return None

    # Find ECMWF file: try local directory first, then DaRUS
    filename = f'{ecmwf_mission_id}_pwv.fits'
    ecmwf_file = None
    if ecmwf_dir is not None:
        if not os.path.isdir(str(ecmwf_dir)):
            log.warning(f'ECMWF directory not found: {ecmwf_dir}')
        else:
            candidate = Path(ecmwf_dir, filename)
            if candidate.is_file():
                ecmwf_file = candidate
            else:
                log.debug(f'ECMWF file not found locally: {candidate}')
    if ecmwf_file is None:
        try:
            ecmwf_file = Path(get_ecmwf_from_darus(filename))
        except Exception as e:
            log.debug(f'Could not retrieve ECMWF file from DaRUS: {e}')
            return None

    WVZECMWF = None
    wvz_fifi = None
    try:
        cached = get_ecmwf_from_cache(ecmwf_file)
        if cached is not None:
            unixtime, pwv, posconst, posjump = cached
            log.debug(f'Retrieved ECMWF data from cache: {ecmwf_file}')
        else:
            with fits.open(ecmwf_file) as ecmwf_hdul:
                data = ecmwf_hdul['DATA'].data
                unixtime = data['unixtime'].copy()
                pwv = data['pwv'].copy()
                posconst = data['posconst'].copy()
                posjump = data['posjump'].copy()
            store_ecmwf_in_cache(ecmwf_file, unixtime, pwv, posconst, posjump)
            log.debug(f'Stored ECMWF data in cache: {ecmwf_file}')

        # Check if observation time is within the data range
        if unixtime[0] < obs_time_unix < unixtime[-1]:
            # Binary search for closest time index (unixtime is sorted)
            ecmwf_idx = np.searchsorted(unixtime, obs_time_unix)
            # Check quality flags (posconst and posjump should be 0)
            if posconst[ecmwf_idx] == 0 and posjump[ecmwf_idx] == 0:
                WVZECMWF = float(pwv[ecmwf_idx])
                # Convert to FIFI-LS scale
                wv_offset, wv_slope = 0.34, 0.55
                wvz_fifi = wv_offset + WVZECMWF * wv_slope
                wv_formula = (f'WVZ_USED = {wv_offset} '
                              f'+ WVZECMWF * {wv_slope}')
                log.debug(f'ECMWF WV: {WVZECMWF:.2f} -> '
                          f'FIFI-LS WV: {wvz_fifi:.2f} '
                          f'({wv_formula})')
    except Exception as e:
        log.debug(f'Error reading ECMWF file {ecmwf_file}: {e}')

    if wvz_fifi is None:
        log.warning('No valid ECMWF data found for this observation')
        return None

    return (WVZECMWF, wvz_fifi, wv_formula, filename)


def get_atran_parameters(header, use_ecmwf, ecmwf_dir):
    """Read FITS Header and get WV for ATRAN file selection.

    Parameters
    ----------
    header : fits.Header
        FITS header containing observation information.
    use_ecmwf : bool
        Whether to attempt to get water vapor from ECMWF data or
        from WVZ_OBS in the header.
    ecmwf_dir : str
        Directory containing ECMWF data files.
    """
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
    wv = None
    WVZECMWF = None
    wvz_fifi = None
    wv_formula = None
    wv_source = 'HEADER'

    if use_ecmwf:
        # Try to get water vapor from ECMWF reanalysis data
        ecmwf_result = get_wv_from_ecmwf(header, ecmwf_dir)
        if ecmwf_result is not None:
            WVZECMWF, wvz_fifi, wv_formula, wv_file = ecmwf_result
            wv = wvz_fifi
            wv_source = 'ECMWF'
            log.info(f'Using ECMWF water vapor: {wv:.2f}')

    if wv is None:
        # Try to get water vapor from the header WVZ_OBS keyword
        wvz_obs = float(header.get('WVZ_OBS', 0))
        if wvz_obs > 0:
            wv = wvz_obs
        else:
            log.warning('WVZ_OBS is missing or invalid in header.')
            if not use_ecmwf:
                log.warning('use_ecmwf=False but no valid WVZ_OBS available. '
                            'Automatically applying use_ecmwf=True')
                ecmwf_result = get_wv_from_ecmwf(header, ecmwf_dir)
                if ecmwf_result is not None:
                    WVZECMWF, wvz_fifi, wv_formula, wv_file = ecmwf_result
                    wv = wvz_fifi
                    wv_source = 'ECMWF'
                    log.info(f'Using ECMWF water vapor: {wv:.2f}')
            if wv is None:
                log.error('No valid water vapor value available. '
                          f'Using minimum WV ({ATRAN_WV_VALUES[0]} um).')
                wv = float(ATRAN_WV_VALUES[0])

    if wv < 1.:
        log.error(f'Invalid water vapor value: {wv}.')

    log.debug(f'Alt, ZA, WV: {alt:.2f} {za:.2f} {wv:.2f}')

    # Add water vapor source and values to header
    hdinsert(header, 'WVZ_SRC', wv_source,
             comment='Source of water vapor value (ECMWF or HEADER)')
    hdinsert(header, 'WVZ_USED', round(wv, 2),
             comment='[um] Water vapor used for ATRAN selection')
    if WVZECMWF is not None:
        hdinsert(header, 'WVZECMWF', round(WVZECMWF, 2),
                 comment='[um] Raw ECMWF WV (before FIFI-LS conversion)')
        hdinsert(header, 'WVZ_FORM', wv_formula,
                 comment='Formula used to convert ECMWF WV to FIFI-LS scale')
        hdinsert(header, 'WVZ_FILE', wv_file,
                 comment='ECMWF PWV file used')

    return alt, za, wv


def get_atran(header, resolution=None, atran_file=None,
              get_unsmoothed=False, atran_dir=None,
              use_ecmwf=True, ecmwf_dir=None, ozone_model=39,
              interpolated=True):
    """
    Retrieve reference atmospheric transmission data.

    ATRAN files in the data/atran_files directory should be named
    according to description in https://doi.org/10.18419/DARUS-5705:

        atran_sdc_[alt]K_[za]deg_[wv]pwv_[ozone]deg_2nlayer_[wmin]-[wmax]mum_bt.fits

    "[ozone]deg" denotes the Ozone column model (9, 30, 39, 43, or 59);
    "2nlayer" (currently fixed) denotes the number of atmospheric
    layers used in ATRAN.
    The FIFI-LS files have a wavelength range of 40 to 300 microns.

    Two selection modes are available via the `interpolated` parameter:

    - ``interpolated=True`` (default): Loads the four ATRAN files that
      bracket the observed ZA and WV values and linearly interpolates
      between them.

    - ``interpolated=False``: Selects the single ATRAN file with the
      nearest matching ZA, altitude, and WV values.

    In both modes, if `atran_file` is provided it is used directly and
    the selection logic is skipped entirely.

    Parameters
    ----------
    header : fits.Header
        ATRNFILE keyword is written to the provided FITS header,
        containing the name of the ATRAN file used.
    resolution : float, optional
        Spectral resolution to which ATRAN data should be smoothed.
    atran_file : str, optional
        Exact path to an ATRAN file to use directly. If provided, the
        normal file selection and interpolation are skipped. If the file
        is not found, a warning is issued and automatic selection is
        used as fallback.
    get_unsmoothed : bool, optional
        If True, return the unsmoothed atran data in addition to the
        smoothed.
    atran_dir : str, optional
        Path to a directory containing ATRAN reference FITS files.
        If not provided, files will be downloaded from the SOFIA
        Astronomy Dataverse
        (https://darus.uni-stuttgart.de/dataverse/irs-sofia-ad).
    use_ecmwf : bool, optional
        Whether to get water vapor from ECMWF data. If False, WVZ_OBS
        from the header is used. If WVZ_OBS is missing, ECMWF is applied
        automatically with a warning.
    ecmwf_dir : str, optional
        Directory containing ECMWF data files, required if use_ecmwf is True.
    ozone_model : int, optional
        Ozone column model in degrees. Currently fixed at 39, which is
        the only model available via DaRUS. Default is 39.
    interpolated : bool, optional
        If True (default), linearly interpolate between the four ATRAN
        files bracketing the observed ZA and WV. If False, use the single
        nearest-matching file.

    Returns
    -------
    numpy.ndarray or tuple of numpy.ndarray
        A (2, nw) array of [wavelengths, smoothed transmission]. If
        `get_unsmoothed` is True, returns a tuple
        ``(smoothed_array, unsmoothed_array)``.
    """
    if not isinstance(header, fits.Header):
        log.error('Invalid header')
        return

    if resolution is None:
        log.warning('Getting default resolution from G_WAVE in header')
        resolution = get_resolution(header)

    # If an exact file is provided, load it directly and skip all selection
    if atran_file is not None:
        if not goodfile(atran_file, verbose=True):
            log.warning(f'File {atran_file} not found; retrieving default')
            atran_file = None
        else:
            log.debug(f'Using exact ATRAN file: {atran_file}')
            atran_data = get_atran_data(atran_file, resolution)
            if atran_data is None:
                return None
            atranfile, wave, unsmoothed, smoothed = atran_data
            hdinsert(header, 'ATRNFILE', atranfile)
            if not get_unsmoothed:
                return np.vstack((wave, smoothed))
            else:
                return (np.vstack((wave, smoothed)),
                        np.vstack((wave, unsmoothed)))

    # Get parameters from header (ZA, altitude, water vapor)
    alt, za, wv = get_atran_parameters(header, use_ecmwf, ecmwf_dir)

    # atran_layers is currently not variable
    O3_model = f'{ozone_model}deg'
    atran_layers = '2nlayer'

    if not interpolated:
        # Nearest-file selection
        alt = int(round(alt))
        za = ATRAN_ZA_VALUES[np.argmin(np.abs(np.array(ATRAN_ZA_VALUES) - za))]
        wv = ATRAN_WV_VALUES[np.argmin(np.abs(np.array(ATRAN_WV_VALUES) - wv))]

        log.debug(f'Using nearest Alt {alt}K, ZA {za}deg, WV {wv}um')

        # see https://doi.org/10.18419/DARUS-5705 for naming convention.
        atran_file = (
            "atran_sdc"
            f"_{alt}K_{za}deg_{wv}pwv"
            f"_{O3_model}_{atran_layers}_40-300mum_bt.fits")

        log.debug(f'Using ATRAN file: {atran_file}')
        atran_data = get_atran_data(atran_file, resolution, atran_dir)
        if atran_data is None:
            return None
        atranfile, wave, unsmoothed, smoothed = atran_data
        hdinsert(header, 'ATRNFILE', atranfile)

    else:
        # interpolation over ZA and WV

        # Clip values to ATRAN data range
        if not ATRAN_ZA_VALUES[-1] >= za >= ATRAN_ZA_VALUES[0]:
            log.warning('za={} outside of available ATRAN data.'.format(za))
            za = np.clip(za, a_min=ATRAN_ZA_VALUES[0],
                         a_max=ATRAN_ZA_VALUES[-1])
            log.warning('Setting zenith angle to {} deg'.format(za))
            za_high, za_low = ATRAN_ZA_VALUES[-1], ATRAN_ZA_VALUES[-2]
        else:
            za_high, za_low = np.inf, np.inf

        if not ATRAN_WV_VALUES[-1] >= wv >= ATRAN_WV_VALUES[0]:
            log.warning('wv={} outside of available ATRAN data.'.format(wv))
            wv = np.clip(wv, a_min=ATRAN_WV_VALUES[0],
                         a_max=ATRAN_WV_VALUES[-1])
            log.warning('Setting water vapor to {} um'.format(wv))
            wv_high, wv_low = ATRAN_WV_VALUES[-1], ATRAN_WV_VALUES[-2]
        else:
            wv_high, wv_low = np.inf, np.inf

        if not 45 >= alt >= 35:
            log.warning('alt={} outside of available ATRAN data.'.format(alt))
            alt = np.clip(alt, a_min=35, a_max=45)
            log.warning('Setting altitude to {}K ft'.format(round(alt)))

        log.debug(f'Alt, ZA, WV: {alt:.2f} {za:.2f} {wv:.2f}')

        # find the higher and lower boundaries
        for i, z in enumerate(ATRAN_ZA_VALUES):
            if z > za:
                za_high = z
                za_low = ATRAN_ZA_VALUES[i-1]
                break

        for i, w in enumerate(ATRAN_WV_VALUES):
            if w > wv:
                wv_high = w
                wv_low = ATRAN_WV_VALUES[i-1]
                break

        # load all files
        grid_data = {}
        atrnfile_fits_keyword = ''
        for key, _za, _wv in (("za1_wv1", za_low,  wv_low),
                               ("za1_wv2", za_low,  wv_high),
                               ("za2_wv1", za_high, wv_low),
                               ("za2_wv2", za_high, wv_high)):
            # see https://doi.org/10.18419/DARUS-5705 for naming convention.
            grid_file = (
                "atran_sdc"
                f"_{round(alt)}K_{_za}deg_{_wv}pwv"
                f"_{O3_model}_{atran_layers}_40-300mum_bt.fits")

            single_atran_data = get_atran_data(grid_file, resolution, atran_dir)
            if single_atran_data is None:
                return None
            grid_data[key] = (single_atran_data, _za, _wv)
            atrnfile_fits_keyword += grid_file + ', '

        # interpolate za for two pwv
        log.debug(f'Interpolating between ZA: {za_low:.2f}, {za_high:.2f}'
                  f' & WV: {wv_low:.2f}, {wv_high:.2f}')
        za1_za2_wv1 = interpolate_two_atran_files(
            grid_data["za1_wv1"], grid_data["za2_wv1"], za, 0)
        za1_za2_wv2 = interpolate_two_atran_files(
            grid_data["za1_wv2"], grid_data["za2_wv2"], za, 0)
        result = interpolate_two_atran_files(za1_za2_wv1, za1_za2_wv2, wv, 1)

        wave = result[0][1]
        unsmoothed = result[0][2]
        smoothed = smoothres(wave, unsmoothed, resolution)
        hdinsert(header, 'ATRNFILE', atrnfile_fits_keyword[:-2])
        log.debug(f'Interpolated ATRAN: alt={round(alt)}K, '
                  f'za={za:.2f}deg, wv={wv:.2f}um')

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
