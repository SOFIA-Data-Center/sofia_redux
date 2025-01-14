# Licensed under a 3-clause BSD style license - see LICENSE.rst

import matplotlib.pyplot as plt
from astropy.io import fits

from sofia_redux.instruments.fifi_ls.get_atran import get_atran_interpolated

def test_get_atran_interpolated():

    # create fake header
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

    # default: gets alt/za/resolution from header
    atran_smoothed, atran_unsmoothed = \
        get_atran_interpolated(header, use_wv=True, get_unsmoothed=True)
    assert atran_smoothed.ndim == 2
    assert atran_unsmoothed.ndim == 2
    # default: no unsmoothed data
    atran_smoothed = get_atran_interpolated(header, use_wv=True)
    assert atran_smoothed.ndim == 2

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
