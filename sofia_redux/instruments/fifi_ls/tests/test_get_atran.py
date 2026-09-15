# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import shutil
import time
from unittest.mock import patch

import numpy as np
import pytest
from astropy.io import fits

from sofia_redux.instruments.fifi_ls.get_atran import (
    clear_atran_cache,
    get_atran,
    get_atran_from_cache,
    get_atran_from_darus,
    get_wv_from_ecmwf,
    store_atran_in_cache,
)


def test_atran_cache(tmpdir):

    tempdir = str(tmpdir.mkdir('test_get_atran'))
    atranfile = os.path.join(tempdir, 'test01')
    res = 100.0

    with open(atranfile, 'w') as f:
        print('this is the atran file', file=f)

    wave = unsmoothed = smoothed = np.arange(10)

    filename = 'ATRNFILE_header_value'
    store = atranfile, res, filename, wave, unsmoothed, smoothed

    clear_atran_cache()
    assert get_atran_from_cache(atranfile, res) is None
    store_atran_in_cache(*store)

    # It should be in there now
    result = get_atran_from_cache(atranfile, res)
    assert result[0] == filename

    for r in result[1:]:
        assert np.allclose(r, wave)

    # Check it's still in there
    assert get_atran_from_cache(atranfile, res) is not None

    # Check that a different resolution does not retrieve this file
    assert get_atran_from_cache(atranfile, res * 2) is None

    # Modify the file - the result should be None,
    # indicating it was removed from the file and
    # should be processed and stored again.
    time.sleep(0.5)
    with open(atranfile, 'w') as f:
        print('a modification', file=f)

    assert get_atran_from_cache(atranfile, res) is None

    # Store the data again
    store_atran_in_cache(*store)

    # Make sure it's there
    assert get_atran_from_cache(atranfile, res) is not None

    # Check clear works
    clear_atran_cache()
    assert get_atran_from_cache(atranfile, res) is None

    # Store then delete the atran file -- check that bad file
    # can't be retrieved
    store_atran_in_cache(*store)
    assert get_atran_from_cache(atranfile, res) is not None
    os.remove(atranfile)
    assert get_atran_from_cache(atranfile, res) is None


def test_filename(tmpdir, capsys, test_files):
    atranfile = tmpdir.join('test_file.fits')

    filename = test_files('scm')[0]
    header = fits.open(filename)[0].header

    # default: gets alt/za/resolution from header, no unsmoothed data
    default = get_atran(header)
    assert default is not None
    assert isinstance(default, np.ndarray)

    # provide a bad filename -- warns and gets default
    result = get_atran(header, atran_file=str(atranfile))
    assert np.allclose(result, default)
    capt = capsys.readouterr()
    assert 'not found; retrieving default' in capt.err

    # provide a good filename, bad file
    atranfile.write('Test data\n')
    result = get_atran(header, atran_file=str(atranfile))
    assert result is None
    capt = capsys.readouterr()
    assert 'Invalid data' in capt.err

    # good filename, minimum required data
    data = np.arange(100).reshape(2, 50).astype(float)
    hdul = fits.HDUList(fits.PrimaryHDU(data=data))
    hdul.writeto(str(atranfile), overwrite=True)
    result = get_atran(header, atran_file=str(atranfile))
    # wavelengths will match
    assert np.allclose(result[0], data[0])
    # data will be all nan
    assert np.all(np.isnan(result[1]))

ATRAN_TEST_NAME = 'atran_sdc_41K_45deg_1pwv_39deg_2nlayer_40-300mum_bt.fits'
ATRAN_TEST_WAVE = np.linspace(40., 300., 50)


@pytest.fixture
def atran_table_file(tmpdir):
    """Make a minimal ATRAN binary table file in a user directory.
    """
    def _write(content='ATRAN_SDC Model'):
        path = tmpdir.join(ATRAN_TEST_NAME)
        columns = fits.ColDefs([
            fits.Column(name='wavelength', format='E', unit='um',
                        array=ATRAN_TEST_WAVE),
            fits.Column(name='transmission', format='E',
                        array=np.ones(ATRAN_TEST_WAVE.size))])
        table = fits.BinTableHDU.from_columns(columns)
        if content is not None:
            table.header['CONTENT'] = content
        fits.HDUList([fits.PrimaryHDU(), table]).writeto(str(path))
        return str(path)

    return _write


def test_atran_file_override_path(capsys, test_files, atran_table_file):
    """An override given as a full path is read from that path directly."""
    filename = test_files('scm')[0]
    header = fits.open(filename)[0].header

    atranfile = atran_table_file()

    result = get_atran(header, atran_file=atranfile)
    assert np.allclose(result[0], ATRAN_TEST_WAVE)

    # the file is used as given and no directory lookup
    capt = capsys.readouterr()
    assert 'Non-standard ATRAN filename' not in capt.out
    assert 'Non-standard ATRAN filename' not in capt.err


def test_atran_file_missing_content(capsys, test_files, atran_table_file):
    """A table without CONTENT warns and falls back, rather than raising."""
    filename = test_files('scm')[0]
    header = fits.open(filename)[0].header

    atranfile = atran_table_file(content=None)

    # no KeyError: falls back to the old format, finds no image data
    assert get_atran(header, atran_file=atranfile) is None
    capt = capsys.readouterr()
    assert 'Falling back to old format' in capt.err
    assert 'Invalid data' in capt.err


def test_header(capsys, test_files):
    filename = test_files('scm')[0]
    header = fits.open(filename)[0].header

    # default
    default = get_atran(header)
    assert default is not None
    assert isinstance(default, np.ndarray)

    # bad header
    result = get_atran(['test'])
    assert result is None
    capt = capsys.readouterr()
    assert 'Invalid header' in capt.err

    # bad alti_sta -- should get end
    hdr = header.copy()
    hdr['ALTI_STA'] = -9999
    result = get_atran(hdr)
    assert np.allclose(result, default)

    # bad alti_end -- should get start
    hdr = header.copy()
    hdr['ALTI_END'] = -9999
    result = get_atran(hdr)
    assert np.allclose(result, default)

    # bad za_start -- should get end
    hdr = header.copy()
    hdr['ZA_START'] = -9999
    result = get_atran(hdr)
    assert np.allclose(result, default)

    # bad za_end -- should get start
    hdr = header.copy()
    hdr['ZA_END'] = -9999
    result = get_atran(hdr)
    assert np.allclose(result, default)


def test_header_wv(capsys, tmp_path, test_files):
    filename = test_files('scm')[0]
    header = fits.open(filename)[0].header

    # default
    default = get_atran(header)
    assert default is not None
    assert isinstance(default, np.ndarray)

    hdr = header.copy()

    # wvz_obs is used if not ecmwf
    hdr['WVZ_OBS'] = 5.0
    get_atran(hdr, use_ecmwf=False)
    capt = capsys.readouterr()
    assert "Alt, ZA, WV: 41.00 45.00 5.00" in capt.out

    # WVZ_OBS invalid with use_ecmwf=False: warns and falls back to ECMWF
    hdr['WVZ_OBS'] = -9999.
    get_atran(hdr, use_ecmwf=False)
    capt = capsys.readouterr()
    assert "WVZ_OBS is missing or invalid in header" in capt.err
    assert "use_ecmwf=False but no valid WVZ_OBS available" in capt.err
    assert "Automatically applying use_ecmwf=True" in capt.err


def test_atran_dir(tmp_path, capsys, test_files):
    filename = test_files('scm')[0]
    header = fits.open(filename)[0].header

    # get from cache or download from DaRUS
    default = get_atran(header)

    # specify bad directory: should get default
    result = get_atran(header, atran_dir='badval')
    assert np.allclose(result, default)
    capt = capsys.readouterr()
    assert 'Cannot find ATRAN directory' in capt.err

    # specify empty directory: should get the same as cached before
    result = get_atran(header, atran_dir=str(tmp_path))
    capt = capsys.readouterr()
    assert 'ATRAN file not found in ATRAN directory' in capt.out
    assert np.allclose(result, default)

    # put 2 files into the atran_directory
    atran_alt_dir = tmp_path / "41K"
    atran_alt_dir.mkdir()
    ref_files = [
        'atran_sdc_41K_45deg_1pwv_39deg_2nlayer_40-300mum_bt.fits',
        'atran_sdc_41K_45deg_5pwv_39deg_2nlayer_40-300mum_bt.fits',
    ]
    for pwvf in ref_files:
        cachefile = get_atran_from_darus(41, pwvf)
        (atran_alt_dir / pwvf).symlink_to(cachefile)

    for i, wv in [(0, 1.), (1, 5.)]:
        hdr = header.copy()
        # wvz_obs is used if not ecmwf
        hdr['WVZ_OBS'] = wv
        result = get_atran(hdr, atran_dir=str(tmp_path), use_ecmwf=False,
                           interpolated=False)
        capt = capsys.readouterr()
        assert f"Using ATRAN file: {ref_files[i]}" in capt.out
        assert 'ATRAN file not found in ATRAN directory' not in capt.out


def test_get_wv_from_ecmwf(tmp_path):
    """Test WV retrieval using a synthetic ECMWF file."""
    from astropy.time import Time

    header = fits.Header()
    header['MISSN-ID'] = '2022-08-25_FI_F906'
    header['DATE-OBS'] = '2022-08-25T09:40:20'

    obs_unix = Time(header['DATE-OBS']).unix
    pwv_value = 6.0

    # Build a ECMWF table with observation time
    unixtime = np.array([obs_unix - 3600, obs_unix, obs_unix + 3600])
    pwv = np.array([5.0, pwv_value, 7.0])
    posconst = np.zeros(3, dtype=np.int32)
    posjump = np.zeros(3, dtype=np.int32)

    data_hdu = fits.BinTableHDU.from_columns([
        fits.Column(name='unixtime', format='D', array=unixtime),
        fits.Column(name='pwv', format='D', array=pwv),
        fits.Column(name='posconst', format='J', array=posconst),
        fits.Column(name='posjump', format='J', array=posjump),
    ], name='DATA')
    ecmwf_file = tmp_path / '906_FI_20220825_pwv.fits'
    fits.HDUList([fits.PrimaryHDU(), data_hdu]).writeto(str(ecmwf_file))

    wv_ecmwf, wv_fifi, wv_formula, filename = get_wv_from_ecmwf(
        header, str(tmp_path))

    assert wv_ecmwf == pwv_value
    assert abs(wv_fifi - (0.34 + pwv_value * 0.55)) < 0.001


def test_missing_directory(caplog):
    """
    Test that function returns None when ECMWF data is unavailable

    When the local directory is missing, the function falls back to DaRUS.
    None is only returned if both the local directory and DaRUS fail.
    """
    header = fits.Header()
    header['MISSN-ID'] = '2022-08-25_FI_F906'
    header['DATE-OBS'] = '2022-08-25T09:40:20'

    darus_fail = patch(
        'sofia_redux.instruments.fifi_ls.get_atran.get_ecmwf_from_darus',
        side_effect=FileNotFoundError("DaRUS unavailable"))

    with darus_fail:
        # Non-existent directory: logs warning - falls back to DaRUS
        # DaRUS fails too then None
        result = get_wv_from_ecmwf(header, '/nonexistentpath')
        assert result is None
        assert "ECMWF directory not found" in caplog.text

        # None directory: goes straight to DaRUS, DaRUS fails then None
        result = get_wv_from_ecmwf(header, None)
        assert result is None


@pytest.fixture(scope='function')
def header_for_atran():
    """Provide a fake header sufficient for atran file selection."""
    header = fits.Header()
    header['ZA_START'] = "45.5"
    header['ZA_END'] = "45.7"
    header['ALTI_STA'] = "41014.0"
    header['ALTI_END'] = "41008.0"
    header['WVZ_OBS'] = "47.9"
    header['G_WAVE_B'] = "51.819"
    header['G_WAVE_R'] = "157.741"
    header['CHANNEL'] = "RED"
    header['G_ORD_B'] = "2"
    return header


@pytest.fixture
def fake_atran_dir(tmp_path):
    """Create a temporary fake ATRAN directory for the whole data range.

    Get ATRAN spectra from the cache or downloads from DaRUS and place
    them into a temporary directory for the atran_dir argument.
    """
    atran_dir = tmp_path / 'atran'
    atran_dir.mkdir()

    for alt in [38, 41, 45]:
        atran_alt_dir = atran_dir / f'{alt}K'
        atran_alt_dir.mkdir()
        for za in [30, 35, 45, 50, 65, 70]:
            for wv in [1, 2, 45, 50]:
                ref_atran = (
                    f'atran_sdc_{alt}K_{za}deg_{wv}pwv'
                    '_39deg_2nlayer_40-300mum_bt.fits'
                )
                cachefile = get_atran_from_darus(alt, ref_atran)
                shutil.copy(cachefile, atran_alt_dir / ref_atran)
    return atran_dir


def test_get_atran(header_for_atran):
    # default: gets alt/za/resolution from header
    atran_smoothed, atran_unsmoothed = get_atran(
        header_for_atran, get_unsmoothed=True
    )
    assert atran_smoothed.ndim == 2
    assert atran_unsmoothed.ndim == 2
    # default: no unsmoothed data
    atran_smoothed = get_atran(
        header_for_atran
    )
    assert atran_smoothed.ndim == 2


@pytest.mark.parametrize(
    "hdrval, expected_atran_string, expected_warnings",
    [
        (
            {},  # no change, standard header, no warning
            "Alt, ZA, WV: 41.01 45.60 47.90",
            None
        ),
        (
            {"ZA_START": "29.0", "ZA_END": "29.2"},
            "Alt, ZA, WV: 41.01 30.00 47.90",
            [
                "za=29.1 outside of available ATRAN data",
                "Setting zenith angle to 30.0 deg"
            ]
        ),
        (
            {"ZA_START": "70.0", "ZA_END": "71.0"},
            "Alt, ZA, WV: 41.01 70.00 47.90",
            [
                "za=70.5 outside of available ATRAN data",
                "Setting zenith angle to 70.0 deg"
            ]
        ),
        (
            {"WVZ_OBS": "55.0"},
            "Alt, ZA, WV: 41.01 45.60 50.00",
            [
                "wv=55.0 outside of available ATRAN data",
                "Setting water vapor to 50.0 um"
            ]
        ),
        (
            {"WVZ_OBS": "0.5"},
            "Alt, ZA, WV: 41.01 45.60 1.0",
            [
                "wv=0.5 outside of available ATRAN data",
                "Setting water vapor to 1.0 um"
            ]
        ),
        (
            {"ALTI_STA": "33900.", "ALTI_END": "35000."},
            "Alt, ZA, WV: 35.00 45.60 47.90",
            [
                "alt=34.45 outside of available ATRAN data",
                "Setting altitude to 35K ft"
            ]
        ),
    ]
)
def test_get_atran_clipped(
        header_for_atran, caplog, fake_atran_dir,
        hdrval, expected_atran_string, expected_warnings):
    """Test that get_atran handles the desired clipping."""
    for k, v in hdrval.items():
        header_for_atran[k] = v

    atran_spectrum = get_atran(
        header_for_atran, atran_dir=fake_atran_dir)

    assert atran_spectrum.ndim == 2

    assert expected_atran_string in caplog.text

    assert "Invalid data in ATRAN file" not in caplog.text

    if not expected_warnings:
        assert "outside of availalbe ATRAN data " not in caplog.text
    else:
        for expected_warning in expected_warnings:
            assert expected_warning in caplog.text
