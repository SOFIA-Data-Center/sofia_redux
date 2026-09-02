============
Installation
============

Stable release
--------------

The `sofia-redux` package is available via PyPI::

   pip install sofia-redux


From source
-----------

Obtain the source code for this package from the `SOFIA Redux GitHub project
<https://github.com/SOFIA-Data-Center/sofia_redux>`__, then install it from the
source directory::

  pip install .

or for development the package can be installed *editable*::

  pip install -e .

Optional Requirements
---------------------

PyQt6
^^^^^

If sofia_redux is installed via pip, the PyQt6 package, required for
the pipeline GUI interface, is not automatically installed as a dependency.
To use the GUI tools, install PyQt6 via pip::

  pip install PyQt6

or use the [display] extra::

  pip install 'sofia-redux[display]'

Please note that there may be some incompatibilities between some versions
of PyQt6, some versions of the package dependencies, and some host OS versions.

DS9
^^^

Some optional visualization tools in the SOFIA Redux interface also
use the `ds9samp` and `regions` packages to interface with the external
SAOImage DS9 tool. To use these tools, install
`DS9 <https://sites.google.com/cfa.harvard.edu/saoimageds9>`__, then
install ds9samp and regions directly via pip::

  pip install ds9samp regions

or using the [display] extra::

  pip install 'sofia-redux[display]'

On MacOS, you will need to make a `ds9`
executable available in your PATH environment variable; see the
`DS9 FAQs <http://ds9.si.edu/doc/faq.html#MacOSX>`__ for more information.

Reference data
^^^^^^^^^^^^^^

Some pipeline modes require additional reference data for optimal data
quality.  These files are too large to distribute with the pipeline code,
so they are provided separately.

Atmospheric models
~~~~~~~~~~~~~~~~~~
For optimal telluric correction, FORCAST, FLITECAM, and FIFI-LS spectroscopic
reductions require a library of FITS files, containing model atmospheric
transmission spectra, derived either from Steve Lord's ATRAN model
`(bibcode: 1992nstc.rept.....L) <https://ntrs.nasa.gov/citations/19930010877>`__.
(older models), or from ATRAN_SDC, a SOFIA-specific version of the former
that was created with permission by the author.
`ATRAN_SDC is provided as a web service <https://atran.sofiadatacenter.de>`
by the SOFIA Data Center where you can generate and download atmospheric
transmission spectra with parameters specific to your reduction analysis.

The EXES pipeline does not use ATRAN models for telluric correction, but it
does attach a reference atmospheric model at a matching altitude and zenith
angle to output spectral products, if available. The models used are derived
from the `Planetary Spectrum Generator (PSG) <https://psg.gsfc.nasa.gov/>`__.

- FORCAST:

  For FORCAST, two versions of the ATRAN libraries are available, with and without water vapor parameterization.

  - Approximate models, not accounting for water vapor variation

    - Download: `atran_forcast_standard.tgz <https://irsa.ipac.caltech.edu/data/SOFIA/ATRAN_FITS/atran_forcast_standard.tgz>`__
    - Size: 531.3 MB
    - MD5 checksum: 61141843b245eea1fdfd45167a1a750b

  - More accurate models, enabling programmatic optimization of
    the telluric correction

    - Download: `atran_forcast_wv.tgz <https://irsa.ipac.caltech.edu/data/SOFIA/ATRAN_FITS/atran_forcast_wv.tgz>`__
    - Size: 37.5 GB
    - MD5 checksum: 49264dfd6c3288af2553f73f8082ec97

- FIFI-LS:

  For FIFI-LS, pre-computed ATRAN_SDC models, parameterized by water vapor, are stored in the DaRUS data repository.
  The Pipeline will attempt to retrieve these automatically, however they can also be downloaded and stored locally
  in a directory referenced by the `atran_dir` or `atran_file` configuration parameters of the `telluric_correct` step.
  ECMWF water vapor data are also stored in DaRUS, and can also either be automatically retrieved
  by the pipeline, or downloaded and locally stored. See the user manual for details.

  - ATRAN_SDC models, organized by flight altitude

    - Link: `DaRUS dataverse <https://darus.uni-stuttgart.de/dataverse/irs-sofia-ad?q=atran&types=datasets>`__
    - Size: 3.6 GB
    - MD5 checksum: see individual files on DaRUS

  - ECMWF water vapor data for all FIFI-LS flights

    - Link: `DaRUS dataverse <https://darus.uni-stuttgart.de/dataset.xhtml?persistentId=doi:10.18419/DARUS-5728>`__
    - Size: 340.6 MB
    - MD5 checksum: see individual files on DaRUS

- FLITECAM:

  For FLITECAM, only approximate ATRAN models, with water vapor parameterized by flight altitude, are available.

  - Approximate models, not accounting for water vapor variation

    - Download: `atran_flitecam_standard.tgz <https://irsa.ipac.caltech.edu/data/SOFIA/ATRAN_FITS/atran_flitecam_standard.tgz>`__
    - Size: 875 MB
    - MD5 checksum: 6576883144bcc381eacdfe16688ad4d2

- EXES:

  - Approximate models, not accounting for water vapor variation

    - Download: `psg_exes_standard.tgz <https://irsa.ipac.caltech.edu/data/SOFIA/ATRAN_FITS/psg_exes_standard.tgz>`__
    - Size: 5.4 GB
    - MD5 checksum: 147cf56cf15f2626b75a600e1ede5410


After downloading and unpacking the library, its location can be provided
to the pipeline as an optional parameter in the *Calibrate Flux* step for
FORCAST or FLITECAM, the *Telluric Correct* step for FIFI-LS, or the
*Extract Spectra* step for EXES.

Standard flux models
~~~~~~~~~~~~~~~~~~~~
In addition to the ATRAN_SDC models, a library of standard flux models is
required to reduce FORCAST or FLITECAM standard spectra to instrumental
response curves. This should be rarely needed for standard scientific reductions,
since reference response curves are provided for most data.  If needed for
re-deriving spectral flux calibrations, the standard model spectra are
provided in the
`source distribution <https://github.com/SOFIA-Data-Center/sofia_redux>`__ of
this package, at sofia_redux/instruments/forcast/data/grism/standard_models
or sofia_redux/instruments/flitecam/data/grism/standard_models.

FLITECAM and EXES auxiliary data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The FLITECAM and EXES pipelines rely upon a number of auxiliary data files.
For FLITECAM, these include nonlinearity correction coefficients,
spectroscopic order masks, and wavelength calibration files. For EXES,
these are bad pixel masks, reset dark files, and nonlinearity correction coefficients.

These datasets are too large to be included in the software packages
provided via PyPI; therefore, they are stored in the DaRUS data repository, and
will be automatically retrieved by the pipeline as needed. Downloaded files
are stored in astropy's on-disk download cache for the sofia_redux package. See
https://docs.astropy.org/en/stable/utils/data.html for how to find and configure the cache location.

Troubleshooting
---------------

The last working set of installed versions of all dependencies is recorded in the
`requirements-min.txt`
file in this package. If errors are encountered in the other listed installation
methods, it may be useful to install the frozen versions directly. For example, to install
from source create a new Python environment from the sofia_redux package
directory::

   python -v venv .redux_troubleshoot_venv
   source .redux_troubleshoot_venv/bin/activate
   python -m pip install -r requirements-min.txt
   python -m pip install -e .


