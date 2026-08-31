# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Download files from DaRUS (darus.uni-stuttgart.de).

DaRUS is the data repository at the University of Stuttgart where the
SOFIA pipeline's reference files (ATRAN, ECMWF PWV, EXES and FLITECAM
calibration files, ...) are hosted. Each instrument's files have their
own dataset, identified by a unique DOI.
"""
from pathlib import Path

from astropy import log
from astropy.utils.data import download_file
import requests

__all__ = ['DARUS_URL_BASE', 'DarusError', 'get_file_from_darus']

DARUS_URL_BASE = "https://darus.uni-stuttgart.de"


class DarusError(OSError):
    """Raised when a file could not be retrieved from DaRUS.

    Carries a single message describing what was being retrieved and why
    it failed. Callers raise it in place of the underlying transport
    exception once that exception has been logged, so that the reduction
    stops with a readable error instead of a chained retry traceback.
    """

# cache for list of files in each DaRUS dataset in the DOI
__darus_files_in_ds = {}


def get_file_from_darus(doi, filename):
    """Download a file from a DaRUS dataset by DOI and filename.

    The downloaded file is stored in astropy's on-disk download cache
    for the sofia_redux package. See
    https://docs.astropy.org/en/stable/utils/data.html
    for information on how to find and configure the cache location.

    Parameters
    ----------
    doi : str
        DOI of the DaRUS dataset, e.g. ``"10.18419/DARUS-5981"``.
    filename : str
        Name of the file to be downloaded

    Returns
    -------
    str
        Local path to the downloaded file.

    Raises
    ------
    FileNotFoundError
        If the dataset does not contain a file with that name.
    requests.HTTPError
        If the dataset listing request fails.
    """
    global __darus_files_in_ds
    if __darus_files_in_ds.get(doi) is None:
        r = requests.get(f'{DARUS_URL_BASE}/api/datasets/:persistentId',
                         params={'persistentId': f'doi:{doi}'})
        r.raise_for_status()
        __darus_files_in_ds[doi] = r.json()['data']['latestVersion']['files']

    for file in __darus_files_in_ds[doi]:
        if file['label'] != filename:
            continue
        fid = file['dataFile']['id']
        download_url = f'{DARUS_URL_BASE}/api/access/datafile/{fid}'
        cached_file = download_file(
            download_url, cache=True, pkgname="sofia_redux")
        log.debug(f"DaRUS file {filename} from {doi} cached in "
                  f"{Path(cached_file).parent}")
        return cached_file

    raise FileNotFoundError(
        f'{filename} not found in DaRUS dataset {doi}')
