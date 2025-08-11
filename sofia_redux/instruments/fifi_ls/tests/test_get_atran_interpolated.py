# Licensed under a 3-clause BSD style license - see LICENSE.rst

from pathlib import Path
import shutil

from astropy.io import fits
import matplotlib.pyplot as plt
import pytest

from sofia_redux.instruments.fifi_ls import get_atran


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

    These are not the real spectra at those altitudes, zenith angles,
    and water vapor values!
    """
    atran_dir = tmp_path / 'atran'
    atran_dir.mkdir()

    real_atran_dir = Path(get_atran.__file__).parent / 'data' / 'atran_files'
    real_atran = [
        real_atran_dir / 'atran_41K_45deg_45pwv_40-300mum.fits',
        real_atran_dir / 'atran_41K_45deg_50pwv_40-300mum.fits'
    ]

    shutil.copy(real_atran[0], atran_dir / real_atran[0].name)
    shutil.copy(real_atran[1], atran_dir / real_atran[1].name)

    for alt in [38, 41, 45]:
        for za in [30, 35, 45, 50, 65, 70]:
            for wv in [1, 2, 45, 50]:
                fake_atran = f"atran_{alt}K_{za}deg_{wv}pwv_40-300mum.fits"
                shutil.copy(real_atran[0], atran_dir / fake_atran)
    return atran_dir


def test_get_atran_interpolated(header_for_atran):
    # default: gets alt/za/resolution from header
    atran_smoothed, atran_unsmoothed = get_atran.get_atran_interpolated(
        header_for_atran, use_wv=True, get_unsmoothed=True
    )
    assert atran_smoothed.ndim == 2
    assert atran_unsmoothed.ndim == 2
    # default: no unsmoothed data
    atran_smoothed = get_atran.get_atran_interpolated(
        header_for_atran, use_wv=True
    )
    assert atran_smoothed.ndim == 2


@pytest.mark.parametrize(
    "hdrval, expected_atran_string, expected_warnings",
    [
        (
            {},  ## no change, standard header, no warning
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
            {"ALTI_STA": "36900.", "ALTI_END": "38000."},
            "Alt, ZA, WV: 38.00 45.60 47.90",
            [
                "alt=37.45 outside of available ATRAN data",
                "Setting altitude to 38K ft"
            ]
        ),
    ]
)
def test_get_atran_interpolated_clipped(
    header_for_atran, capsys, caplog, fake_atran_dir,
    hdrval, expected_atran_string, expected_warnings):
    """Test that get_atran_interpolated handles the desired clipping."""
    for k, v in hdrval.items():
        header_for_atran[k] = v

    atran_spectrum = get_atran.get_atran_interpolated(
        header_for_atran, use_wv=True, atran_dir=fake_atran_dir)

    assert atran_spectrum.ndim == 2

    capt = capsys.readouterr()
    assert expected_atran_string in capt.out

    assert "Invalid data in ATRAN file" not in caplog.text

    if not expected_warnings:
        assert "outside of availalbe ATRAN data " not in caplog.text
    else:
        for expected_warning in expected_warnings:
            assert expected_warning in caplog.text


def plot_single_atran_file(filename):
    with fits.open(filename) as hdul:
        # hdul.info()
        wavelength = hdul[0].data[0,:]
        transmission = hdul[0].data[1,:]

        print("plotting file ", filename)
        plt.plot(wavelength, transmission, label=filename)


def plot_all_atran_files(filenames):
    [plot_single_atran_file(fn) for fn in filenames]

    plt.title("Transmission over Wavelength")
    plt.xlabel(r"Wavelength [$\mu$m]")
    plt.ylabel("Transmission")
    plt.grid(True)
    plt.legend(loc='lower left')
    plt.xlim(63.,63.5)

    plt.show()


if __name__ == "__main__":
    test_get_atran_interpolated()
