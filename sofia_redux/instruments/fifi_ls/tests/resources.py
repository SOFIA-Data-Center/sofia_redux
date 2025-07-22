# Licensed under a 3-clause BSD style license - see LICENSE.rst

import glob
import os
from pathlib import Path

from astropy import log
from astropy.io import fits
import numpy as np

from sofia_redux.instruments import fifi_ls
from sofia_redux.instruments.fifi_ls.make_header import make_header

seed = 42

# classes used to mock a FIFI data HDU
class MockCols(object):
    def __init__(self):
        self.names = ['DATA', 'HEADER']


class MockData(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.columns = MockCols()


class MockHDU(object):
    def __init__(self):
        self.data = MockData()
        self.header = {}
        self.data['DATA'] = np.zeros((10, 10))
        self.data['HEADER'] = np.chararray((10, 10))


def raw_testdata(nod='A', obsid=None):
    """
    Make raw FIFI-LS data for testing purposes

    Returns
    -------
    HDUList
    """
    # seed the random module for consistent tests
    global seed
    rand = np.random.RandomState(seed)
    # modify for next call to ensure some difference between data sets
    seed += 1

    n = 3840
    data = np.recarray(n, dtype=[('header', '>i2', (8,)),
                                 ('data', '>i2', (18, 26))])

    # header array:
    # channel, sample, ramp (index 3, 4, 5) matter
    # for chop split
    data['header'][:, 3] = 1
    if nod == 'A':
        floor = 500
    else:
        floor = 0
    for i in range(n // 32):
        i32 = i * 32
        data['header'][i32:i32 + 32, 4] = np.arange(32)
        data['header'][i32:i32 + 32, 5] = i
        # put a noisy line in each ramp
        data['data'][i32:i32 + 32, :, :] = \
            (250 * rand.rand()
             * np.arange(32) - 20000 + floor)[:, None, None] \
            + rand.randint(10, size=(32, 18, 26))

    # set empty spexel to zero
    data['data'][:, 0, :] = -2**15

    header = fits.Header()
    header['SIMPLE'] = True
    header['EXTEND'] = True

    # different values for A and B
    header['NODBEAM'] = nod
    if nod == 'A':
        if obsid is None:
            header['OBS_ID'] = 'R001'
        else:
            header['OBS_ID'] = obsid
        header['DATE-OBS'] = '2016-03-01T10:38:39'
        header['UTCSTART'] = '10:38:39'
        header['UTCEND'] = '10:38:44'
    else:
        if obsid is None:
            header['OBS_ID'] = 'R002'
        else:
            header['OBS_ID'] = obsid
        header['DATE-OBS'] = '2016-03-01T10:38:51'
        header['UTCSTART'] = '10:38:51'
        header['UTCEND'] = '10:38:56'

    # standardize some defaults
    header = make_header([header])

    # reduction required keys
    header['CHANNEL'] = 'RED'
    header['C_SCHEME'] = '2POINT'
    header['C_AMP'] = 60.0
    header['C_CHOPLN'] = 64
    header['RAMPLN_B'] = 32
    header['RAMPLN_R'] = 32
    header['G_PSUP_B'] = 4
    header['G_PSUP_R'] = 4
    header['G_SZUP_B'] = 200
    header['G_SZUP_R'] = 510
    header['G_STRT_B'] = 713285
    header['G_STRT_R'] = 463923
    header['G_PSDN_B'] = 0
    header['G_PSDN_R'] = 0
    header['G_SZDN_B'] = 0
    header['G_SZDN_R'] = 0
    header['G_CYC_B'] = 1
    header['G_CYC_R'] = 1
    header['C_CYC_B'] = 10
    header['C_CYC_R'] = 10
    header['G_ORD_B'] = 2
    header['PRIMARAY'] = 'BLUE'
    header['DICHROIC'] = 105
    header['NODSTYLE'] = 'NMC'
    header['DLAM_MAP'] = -14.1
    header['DBET_MAP'] = -5.1
    header['DLAM_OFF'] = 0.
    header['DBET_OFF'] = 0.
    header['OBSLAM'] = 0.
    header['OBSBET'] = 0.
    header['OBSRA'] = 0.
    header['OBSDEC'] = 0.
    header['OBJ_NAME'] = 'Mars'
    header['OBJECT'] = 'Mars'
    header['G_WAVE_B'] = 51.807
    header['G_WAVE_R'] = 162.763
    header['TELRA'] = 15.7356
    header['TELDEC'] = -18.4388
    header['PLATSCAL'] = 4.2331334
    header['DET_ANGL'] = 70.0

    # other important keys
    header['AOR_ID'] = '90_0001_01'
    header['MISSN-ID'] = '2016-03-01_FI_F282'
    header['ALTI_STA'] = 41000.
    header['ALTI_END'] = 41000.
    header['ZA_START'] = 45.
    header['ZA_END'] = 45.
    header['LAT_STA'] = 40.
    header['LON_STA'] = -120.
    header['CHPFREQ'] = 1.
    header['DATASRC'] = 'ASTRO'
    header['DETCHAN'] = 'RED'
    header['EXPTIME'] = 2.56
    header['OBSTYPE'] = 'OBJECT'
    header['SPECTEL1'] = 'FIF_BLUE'
    header['SPECTEL2'] = 'FIF_RED'
    header['ALPHA'] = header['EXPTIME'] / n
    header['START'] = 1456857531.0
    header['FIFISTRT'] = 0

    hdulist = fits.HDUList([fits.PrimaryHDU(header=header),
                            fits.BinTableHDU(data)])
    hdulist[1].header['EXTNAME'] = 'FIFILS_rawdata'
    return hdulist

FIFI_LS_TESTPATH = Path(__file__).parent / 'data'
FIFI_LS_TESTFILES = {
    'raw': [
        '00001_123456_00001_TEST_A_lw.fits',
        '00002_123456_00001_TEST_B_lw.fits',
        '00003_123456_00001_TEST_A_lw.fits',
        '00004_123456_00001_TEST_B_lw.fits',
    ],
    'split': [
        'F0282_FI_IFS_90000101_RED_CP0_001.fits',
        'F0282_FI_IFS_90000101_RED_CP0_002.fits',
        'F0282_FI_IFS_90000101_RED_CP0_003.fits',
        'F0282_FI_IFS_90000101_RED_CP0_004.fits',
        'F0282_FI_IFS_90000101_RED_CP1_001.fits',
        'F0282_FI_IFS_90000101_RED_CP1_002.fits',
        'F0282_FI_IFS_90000101_RED_CP1_003.fits',
        'F0282_FI_IFS_90000101_RED_CP1_004.fits',
    ],
    'chop': [
        'F0282_FI_IFS_90000101_RED_RP0_001.fits',
        'F0282_FI_IFS_90000101_RED_RP0_002.fits',
        'F0282_FI_IFS_90000101_RED_RP0_003.fits',
        'F0282_FI_IFS_90000101_RED_RP0_004.fits',
    ],
    'chop1': [  # not queried but used in subtraction
        'F0282_FI_IFS_90000101_RED_RP1_001.fits',
        'F0282_FI_IFS_90000101_RED_RP1_002.fits',
        'F0282_FI_IFS_90000101_RED_RP1_003.fits',
        'F0282_FI_IFS_90000101_RED_RP1_004.fits',
    ],
    'csb': [
        'F0282_FI_IFS_90000101_RED_CSB_001.fits',
        'F0282_FI_IFS_90000101_RED_CSB_002.fits',
        'F0282_FI_IFS_90000101_RED_CSB_003.fits',
        'F0282_FI_IFS_90000101_RED_CSB_004.fits',
    ],
    'ncm': [
        'F0282_FI_IFS_90000101_RED_NCM_001-002.fits',
        'F0282_FI_IFS_90000101_RED_NCM_003-004.fits',
    ],
    'wav': [
        'F0282_FI_IFS_90000101_RED_WAV_001-002.fits',
        'F0282_FI_IFS_90000101_RED_WAV_003-004.fits',
    ],
    'xyc': [
        'F0282_FI_IFS_90000101_RED_XYC_001-002.fits',
        'F0282_FI_IFS_90000101_RED_XYC_003-004.fits',
    ],
    'flf': [
        'F0282_FI_IFS_90000101_RED_FLF_001-002.fits',
        'F0282_FI_IFS_90000101_RED_FLF_003-004.fits',
    ],
    'scm': [
        'F0282_FI_IFS_90000101_RED_SCM_001-002.fits',
        'F0282_FI_IFS_90000101_RED_SCM_003-004.fits',
    ],
    'tel': [
        'F0282_FI_IFS_90000101_RED_TEL_001-002.fits',
        'F0282_FI_IFS_90000101_RED_TEL_003-004.fits',
    ],
    'cal': [
        'F0282_FI_IFS_90000101_RED_CAL_001-002.fits',
        'F0282_FI_IFS_90000101_RED_CAL_003-004.fits',
    ],
    'wsh': [
        'F0282_FI_IFS_90000101_RED_WSH_001-002.fits',
        'F0282_FI_IFS_90000101_RED_WSH_003-004.fits',
    ],
    'wxy': [
        'F0282_FI_IFS_90000101_RED_WXY_001-004.fits',
    ],
}
FIFI_LS_TESTFILES = {
    k: [str(FIFI_LS_TESTPATH / f) for f in v]
    for k, v in FIFI_LS_TESTFILES.items()
}
FIFI_LS_FILEPATTERN = {
    'raw': '0*',
    'split': '*_CP*_*',
    'chop': '*_RP0_*',
    'chop1': '*_RP1_*',
    'csb': '*_CSB_*',
    'ncm': '*_NCM_*',
    'wav': '*_WAV_*',
    'xyc': '*_XYC_*',
    'flf': '*_FLF_*',
    'scm': '*_SCM_*',
    'tel': '*_TEL_*',
    'cal': '*_CAL_*',
    'wsh': '*_WSH_*',
    'wxy': '*_WXY_*',
}


def create_files():
    """Create raw data files and products for testing."""
    for i in range(4):
        obsid = 'R{:03d}'.format(i + 1)
        if i % 2 == 0:
            nod = 'A'
        else:
            nod = 'B'
        fn = '0000{}_123456_00001_TEST_{}_lw.fits'.format(i + 1, nod)

        hdul = raw_testdata(nod=nod, obsid=obsid)
        hdul[0].header['FILENAME'] = fn

        if i > 1:
            hdul[0].header['DLAM_MAP'] *= -1
            hdul[0].header['DBET_MAP'] *= -1

        hdul.writeto(FIFI_LS_TESTPATH / fn, overwrite=True)

    # run default steps in order
    current_dir = os.getcwd()
    FIFI_LS_TESTPATH.mkdir(exist_ok=True)
    os.chdir(FIFI_LS_TESTPATH)

    input_files = FIFI_LS_TESTFILES['raw']
    from sofia_redux.instruments.fifi_ls.split_grating_and_chop \
        import wrap_split_grating_and_chop
    wrap_split_grating_and_chop(input_files, write=True)

    split_files = FIFI_LS_TESTFILES['split']
    from sofia_redux.instruments.fifi_ls.fit_ramps \
        import wrap_fit_ramps
    wrap_fit_ramps(split_files, write=True)

    ramp0_files = FIFI_LS_TESTFILES['chop']
    from sofia_redux.instruments.fifi_ls.subtract_chops \
        import wrap_subtract_chops
    wrap_subtract_chops(ramp0_files, write=True)

    csb_files = FIFI_LS_TESTFILES['csb']
    from sofia_redux.instruments.fifi_ls.combine_nods \
        import combine_nods
    combine_nods(csb_files, write=True)

    ncm_files = FIFI_LS_TESTFILES['ncm']
    from sofia_redux.instruments.fifi_ls.lambda_calibrate \
        import wrap_lambda_calibrate
    wrap_lambda_calibrate(ncm_files, write=True)

    wav_files = FIFI_LS_TESTFILES['wav']
    from sofia_redux.instruments.fifi_ls.spatial_calibrate \
        import wrap_spatial_calibrate
    wrap_spatial_calibrate(wav_files, rotate=True, write=True)

    xyc_files = FIFI_LS_TESTFILES['xyc']
    from sofia_redux.instruments.fifi_ls.apply_static_flat \
        import wrap_apply_static_flat
    wrap_apply_static_flat(xyc_files, write=True)

    flf_files = FIFI_LS_TESTFILES['flf']
    from sofia_redux.instruments.fifi_ls.combine_grating_scans \
        import wrap_combine_grating_scans
    wrap_combine_grating_scans(flf_files, write=True)

    scm_files = FIFI_LS_TESTFILES['scm']
    from sofia_redux.instruments.fifi_ls.telluric_correct \
        import wrap_telluric_correct
    wrap_telluric_correct(scm_files, write=True)

    tel_files = FIFI_LS_TESTFILES['tel']
    from sofia_redux.instruments.fifi_ls.flux_calibrate \
        import wrap_flux_calibrate
    wrap_flux_calibrate(tel_files, write=True)

    cal_files = FIFI_LS_TESTFILES['cal']
    from sofia_redux.instruments.fifi_ls.correct_wave_shift \
        import wrap_correct_wave_shift
    wrap_correct_wave_shift(cal_files, write=True)

    wsh_files = FIFI_LS_TESTFILES['wsh']
    from sofia_redux.instruments.fifi_ls.resample import resample
    resample(wsh_files, write=True)

    # cd back to original directory
    os.chdir(current_dir)


def check_test_files():
    """Check whether test data files exist in the current pipeline version."""
    for prodtype, expected_files in FIFI_LS_TESTFILES.items():
        existing_files = glob.glob(f'{FIFI_LS_TESTPATH}'
                                   f'/{FIFI_LS_FILEPATTERN[prodtype]}.fits')
        if sorted(existing_files) != expected_files:
            log.debug(f'Expected test files: {expected_files}')
            log.debug(f'Found test files: {existing_files}')
            return False
    expected_pipevers =  fifi_ls.__version__.replace('.', '_')
    pipevers = fits.getval(FIFI_LS_TESTFILES['raw'][0], 'PIPEVERS')
    if pipevers != expected_pipevers:
        log.debug(f'Expected pipeline version: {expected_pipevers}')
        log.debug(f'Found pipeline version: {pipevers}')
        return False
    return True
