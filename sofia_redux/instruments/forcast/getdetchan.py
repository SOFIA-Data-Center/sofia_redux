# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy import log
from astropy.io.fits.header import Header

__all__ = ['getdetchan']


def getdetchan(header):
    """
    Retrieve DETCHAN keyword value from header as either SW or LW

    A helper function to support the change to the DETCHAN keyword from
    value 0|1 to value SW|LW.  Returns SW if any other value is found.

    Parameters
    ----------
    header : astropy.io.fits.header.Header
        FITS header

    Returns
    -------
    str
        SW or LW
    """
    if not isinstance(header, Header):
        log.warning("not a valid header: returning SW")
        return 'SW'

    if 'DETCHAN' not in header:
        log.warning("DETCHAN missing; using SW — wrong "
                    "channel-specific bias level or non-linearity "
                    "correction may be applied.")
        return 'SW'

    value = header.get('DETCHAN', 'SW')
    value = str(value).strip().upper()
    if value == '1':
        return 'LW'
    elif value == 'LW':
        return 'LW'
    elif value in ('0', 'SW'):
        return 'SW'
    else:
        log.warning("DETCHAN value '%s' not recognized; using SW — "
                    "wrong channel-specific bias level or "
                    "non-linearity correction may be applied." % value)
        return 'SW'
