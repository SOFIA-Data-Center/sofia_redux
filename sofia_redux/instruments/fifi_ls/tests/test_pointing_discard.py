# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

from astropy.io import fits
import numpy as np

from sofia_redux.instruments.fifi_ls.pointing_discard \
    import get_timing, load_pointing_data, get_pointing_mask
from sofia_redux.instruments.fifi_ls.tests.resources \
    import FIFITestCase, get_split_files


class TestFitRamps(FIFITestCase):

    def test_get_timing(self):
        filename = get_split_files()[0]
        hdul = fits.open(filename)
        primehead = hdul[0].header

        dt_ramp, dt_grating = get_timing(primehead)

        np.testing.assert_almost_equal(dt_ramp, 0.128)
        np.testing.assert_almost_equal(dt_grating, 5.12)

    def test_load_pointing_data(self):
        filename = get_split_files()[0]
        hdul = fits.open(filename)
        primehead = hdul[0].header

        # Path to simulated test pointing data
        pointing_directory = os.path.join(os.path.dirname(__file__), 'data',
                                          'pointing_error')

        df_pointing_error = load_pointing_data(primehead, pointing_directory)

        assert df_pointing_error.index.name == 'timestamp'
        assert all(df_pointing_error.columns == ['error_RA', 'error_Dec',
                                                 'error_RoF'])
        assert len(df_pointing_error) == 1001

    def test_get_pointing_mask(self):
        filename = get_split_files()[0]
        hdul = fits.open(filename)
        primehead = hdul[0].header

        dt_ramp, dt_grating = get_timing(primehead)

        # Path to simulated test pointing data
        pointing_directory = os.path.join(os.path.dirname(__file__), 'data',
                                          'pointing_error')
        df_pointing_error = load_pointing_data(primehead, pointing_directory)

        pointing_threshold = 15

        for grating_idx in range(2):
            pointing_mask = get_pointing_mask(primehead, pointing_threshold,
                                              grating_idx, dt_ramp, dt_grating,
                                              df_pointing_error)

            if grating_idx == 0:
                mask_grating_pos0 = [True, True, True, True, True, True, True,
                                     True, True, True, False, False, False,
                                     False, False, False, True, True, True,
                                     True]
                assert all(pointing_mask == mask_grating_pos0)
            else:
                mask_grating_pos1 = [True, True, False, False, True, True,
                                     True, True, True, True, True, True, True,
                                     True, True, True, True, True, True, True]
                assert all(pointing_mask == mask_grating_pos1)
