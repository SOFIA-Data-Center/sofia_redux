# Licensed under a 3-clause BSD style license - see LICENSE.rst

from pathlib import Path
import shutil

import pytest
from astropy import log
from filelock import FileLock

from sofia_redux.instruments.fifi_ls.tests.resources import (FIFI_LS_TESTFILES,
                                                             FIFI_LS_TESTPATH,
                                                             check_test_files,
                                                             create_files)


@pytest.fixture(autouse=True, scope='function')
def set_debug_level():
    """Sets the log level to DEBUG on test entry."""
    orig_level = log.level
    log.setLevel('DEBUG')
    # let tests run
    yield
    # reset log level
    log.setLevel(orig_level)


@pytest.fixture(scope="function")
def test_files(tmp_path):
    """Return a function providing test data files, (re-)create if necessary."""
    with FileLock(FIFI_LS_TESTPATH / 'generating.lock'):
        # Use a lock to prevent multiple processes from checking and re-creating
        # at the same time
        if not check_test_files():
            for files in FIFI_LS_TESTFILES.values():
                for file in files:
                    Path(file).unlink(missing_ok=True)
            log.info('Creating test files')
            create_files()
        if not check_test_files():
            raise RuntimeError('Test data files could not be created')

    def get_test_files(prodtype):
        test_files = []
        for file in FIFI_LS_TESTFILES[prodtype]:
            # copy files to temporary directory
            dst = tmp_path / Path(file).name
            shutil.copy(file, dst)
            if prodtype == 'chop':
                # also copy the chop pair
                chop_pair = file.replace('RP0', 'RP1')
                dst_pair = tmp_path / Path(chop_pair).name
                shutil.copy(chop_pair, dst_pair)
            test_files.append(str(dst))
        return test_files

    return get_test_files
