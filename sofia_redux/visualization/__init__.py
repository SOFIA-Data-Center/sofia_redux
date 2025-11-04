# Licensed under a 3-clause BSD style license - see LICENSE.rst

from .version import version as __version__

from sofia_redux.visualization.utils.logger import _init_log
log = _init_log()

__all__ = ['log']
