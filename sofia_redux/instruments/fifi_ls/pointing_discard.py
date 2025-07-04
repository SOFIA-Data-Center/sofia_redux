import os

from astropy.time import Time
import numpy as np
import pandas as pd


def get_timing(primehead):
    '''
    Set up the pointing-discard-feature by calculating the duration of one ramp
    and grating position. The duration of the ramp and grating position is
    calculated by extracting the corresponding keywords from the FITS header of
    the currently loaded FITS file.

    Parameters
    ----------
    primehead : astropy.io.fits.Header
        Primary header of the loaded FITS file.

    Returns
    -------
    dt_ramp, dt_grating : float, float
        Duration in seconds of one ramp and one grating position.
    '''
    # Calculate the number of ramps per grating position
    ramps_per_grating = 2 * primehead['C_CHOPLN'] \
        / primehead['RAMPLN_{}'.format(primehead['CHANNEL'][0])] \
        * primehead['C_CYC_{}'.format(primehead['CHANNEL'][0])]

    # Calculate the time interval in seconds of one ramp
    dt_ramp = primehead['RAMPLN_{}'.format(primehead['CHANNEL'][0])] / 250

    # Calculate the time interval in seconds  of one grating position
    dt_grating = ramps_per_grating * dt_ramp

    return dt_ramp, dt_grating

def load_pointing_data(primehead, pointing_directory=None):
    '''
    Load the file listing the pointing errors. The file must be in the “.tsv”
    format and named as the mission ID, which is extracted from the FITS header.
    The default path to the pointing error file is “./data/pointing_error/” but
    this can be adjusted with the “pointing_directory” parameter.

    Parameters
    ----------
    primehead : astropy.io.fits.Header
        Primary header of the loaded FITS file.
    pointing_directory : str, optional
        Path to the directory containing the pointing error files.

    Returns
    -------
    df_pointing_error : pandas.DataFrame
        Dataframe listing the pointing errors at a given time.
    '''
    # Load pointing error information
    if pointing_directory is None:
        pointing_directory = os.path.join(os.path.dirname(__file__), 'data',
                                          'pointing_error')

    pointing_err_file = os.path.join(pointing_directory,
                                     '{}.tsv'.format(primehead['missn-id']))
    df_pointing_error = pd.read_csv(pointing_err_file,
                                    sep='\t', index_col='timestamp')

    return df_pointing_error

def get_pointing_mask(primehead, pointing_threshold, grating_idx, dt_ramp,
                      dt_grating, df_pointing_error):
    '''
    Create a boolean mask to set ramp values to NaN if the pointing at the ramp
    time is below the defined threshold.

    Parameters
    ----------
    primehead : astropy.io.fits.Header
        Primary header of the loaded FITS file.
    pointing_threshold : astropy.io.fits.Header
        Absolute threshold of the pointing error in arcsec. The slopes of ramps
        with a pointing error greater than the define threshold are set to NaN.
    grating_idx : int
        Index of the grating position.
    dt_ramp : float
        Duration of one ramp in seconds.
    dt_grating : float
        Duration of one grating position in seconds.
    df_pointing_error : pandas.DataFrame
        Dataframe listing the pointing errors at a given time.

    Returns
    -------
    pointing_mask : pandas.Series
        Boolean mask used to discard ramps with bad pointing.
    '''
    # Calculate timestamp of each ramp in the current grating position
    ramp_times = np.arange(round(dt_grating/dt_ramp))*dt_ramp + dt_ramp/2 \
        + grating_idx * dt_grating + Time(primehead['DATE-OBS']).unix

    # Calculate the combined pointing error of RA and Dec
    df_pointing_error['error_combined'] = pd.Series(
        np.sqrt(df_pointing_error.error_RA**2
                + df_pointing_error.error_Dec**2))

    # Calculate the mean ramp time for each chop cycle
    n_ramps_chop_cyc = int(2 * primehead['C_CHOPLN'] \
        / primehead['RAMPLN_{}'.format(primehead['CHANNEL'][0])])
    chop_cycle_times = np.mean(ramp_times.reshape(-1,n_ramps_chop_cyc), axis=1)

    # Create an empty dataframe with the calculated chop cycle times.
    # The pointing error columns are filled with NaNs.
    df_chop_cycle_times = pd.DataFrame(index=chop_cycle_times,
                                       columns=df_pointing_error.columns)

    # Merge the new empty dataframe with the loaded pointing error dataframe
    df_merged = pd.concat([df_pointing_error, df_chop_cycle_times])

    # Interpolate the pointing error at the mean chop cycle times
    df_merged = df_merged.sort_index().interpolate('index', limit_area='inside')

    # Remove all rows but the rows listing the chop cycle times
    df_chop_cycle_pointing_errors = df_merged.loc[chop_cycle_times]

    # Create boolean pointing mask
    pointing_mask = df_chop_cycle_pointing_errors.error_combined \
        < pointing_threshold

    # Repeat the pointing mask to adjust shape
    pointing_mask = np.repeat(pointing_mask, 2)

    return pointing_mask
