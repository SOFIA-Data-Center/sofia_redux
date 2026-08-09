# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy import log
import numpy as np
from sofia_redux.toolkit.utilities.fits import gethdul
from sofia_redux.toolkit.image.adjust import rotate90 as idlrot

__all__ = ['readwavecal']


def readwavecal(filename, rotate=None, info=None):
    """
    Read a Spextool wavecal file

    Parameters
    ----------
    filename : str
        The path to a 2-d Spextool wavecal file.
    rotate : int, optional
        Rotation direction (not angle) to pass to rotate wavecal
        and spatcal images when passed to
        `sofia_redux.toolkit.image.adjust.rotate90`.
    info : dict, optional
        If supplied will be updated with wctype, wavefmt, spatfmt,
        wdisp, flatname, orders, and norders.

    Returns
    -------
    wavecal, spatcal : numpy.ndarray, numpy.ndarray
        - wavecal (nrow, ncol) array where each pixel is set to its
          wavelength (column in this case).
        - spatcal (nrow, ncol) array where each pixel is set to its
          angular position on the sky.
    """
    hdul = gethdul(filename, verbose=True)
    if hdul is None:
        return

    header = hdul[0].header

    wavecal, spatcal = hdul[0].data[:2].copy()
    if rotate is not None:
        wavecal = idlrot(wavecal, rotate)
        spatcal = idlrot(spatcal, rotate)

    orders = np.array(header.get('ORDERS', '0').split(',')).astype(int)
    norders = int(header.get('NORDERS', 0))
    if 'ORDERS' not in header:
        log.warning("ORDERS missing from header; using '0' - orders will "
                    "contain a single placeholder order (0) instead of the "
                    "wavecal file's real orders.")
    if 'NORDERS' not in header:
        log.warning('NORDERS missing from header; using 0 - no per-order '
                    'wavelength/spatial coefficients will be extracted, '
                    'even if the header defines them.')

    if isinstance(info, dict):
        info['wctype'] = str(header.get('WCTYPE')).strip()
        info['orders'] = orders
        info['norders'] = norders
        keys = ['wdeg', 'odeg', 'homeordr', 'wxdeg',
                'wydeg', 'sxdeg', 'sydeg']
        for key in keys:
            info[key] = header.get(key.upper(), 0)
        if (info['wdeg'] + info['odeg']) > 0:
            info['xo2w'] = np.array(list(header['1DWC*'].values()))
        else:
            info['xo2w'] = np.zeros(0)

        if np.sum([info[x] for x in
                   ['wxdeg', 'wydeg', 'sxdeg', 'sydeg']]) > 0:
            nwxy = (info['wxdeg'] + 1) * (info['wydeg'] + 1)
            nsxy = (info['sxdeg'] + 1) * (info['sydeg'] + 1)
            info['xy2w'] = np.zeros((norders, nwxy))
            info['xy2s'] = np.zeros((norders, nsxy))
            for i in range(norders):
                ordi = str(orders[i]).zfill(2)
                info['xy2w'][i] = np.array(
                    list(header['OR%sWC*' % ordi].values()))
                info['xy2s'][i] = np.array(
                    list(header['OR%sSC*' % ordi].values()))
        else:
            info['xy2w'] = np.zeros((norders, 0))
            info['xy2s'] = np.zeros((norders, 0))

    return wavecal, spatcal
