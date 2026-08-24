1.5.2 (2026-XX-YY)
==================

- FLITECAM pipeline v2.2.1:

  - Numpy 2.5 compatibility: removed np.row_stack usage

- FIFI-LS pipeline v2.11.1:

  - Numpy 2.5 compatibility: removed test mockup for np.char.chararray usage
    in FITS headers.

- toolkit module v2.6.0:

  - Numpy 2.5 compatibility: removed unused utilities for np.char.chararray
    to FITS header conversion and vice versa.


1.5.1 (2026-07-22)
==================

- Reformat and release PDF versions of Pipeline User's Manuals:

  - SDC-MAN-0005_Iss01 EXES Pipeline User's Manual
  - SDC-MAN-0006_Iss01 FLITECAM Pipeline User's Manual
  - SDC-MAN-0007_Iss01 FORCAST Pipeline User's Manual
  - SDC-MAN-0008_Iss01 EXES Pipeline User's Manual
  - Hyperlinks to PDF and HTML versions of the manual in the
    `Redux online documentation <https://redux.sofiadatacenter.de/stable/sofia_redux/index.html#manuals>`__.

- calibration module v1.7.0

  - Fix for changed behavior in Astroquery JPL Horizons minor body queries

- pipeline module v2.16.0

  - Clarify FIFI-LS pipeline parameter descriptions

- spectroscopy module v1.5.0
  - Remove custom error handling around Astropy IERS data queries


1.5.0 (2026-06-29)
==================

- FIFI-LS pipeline v2.11.0:

  - Document changes to spectral/spatial resolution data file,
    and clarify functionality of resample.
  - Enhanced telluric correction to automatically retrieve ATRAN_SDC
    transmission spectra and ECMWF water vapor data from the DaRUS Data Repository.
  - The pipeline now supports ATRAN_SDC and ECMWF binary table data formats.
  - Implemented SDC-derived flux response curves, based on an updated
    flux calibration procedure.
  - As part of the new flux calibration, systematic errors are now stored
    in their own dedicated FITS extension, separate from the statistical errors.
  - Fixed a major issue in the B nod interpolation option of the combine nods step.
  - Bumped FIFI-LS User Manual from Iss. 03 to Iss. 04.

- FORCAST pipeline v2.7.3:

  - Maintain compatibility with Photutils 3.0 in the FORCAST peak finding code,
    specifically the
    `rename of column names <https://photutils.readthedocs.io/en/stable/whats_new/3.0.html#deprecated-table-column-names>`__.
    This is a breaking change, thus the minimum required version of Photutils
    has been updated to 3.0.

- EXES pipeline v3.1.0:

  - Download and cache bad pixel masks and dark files from DaRUS if not
    provided locally. Some of these files have been hosted in the Git
    repository but were not included in the published PyPI packages. While the
    README files directed the user to download directly from Github, the
    pipeline had functionality to download missing files from IRSA. However,
    the
    `IRSA Dataset <https://irsa.ipac.caltech.edu/data/SOFIA/PIPELINE_REFERENCE/EXES/>`__
    is not complete. The new
    `SOFIA Redux EXES pipeline reference files <https://doi.org/10.18419/DARUS-5981>`__
    dataset at DaRUS contains the missing files, so we removed them from the
    code repository.

- FLITECAM pipeline v2.1.0:

  - Download and cache auxiliary calibration files from the
    `SOFIA Redux FLITECAM pipeline reference files <https://doi.org/10.18419/DARUS-6126>`__
    dataset at DaRUS instead of
    `IRSA <https://irsa.ipac.caltech.edu/data/SOFIA/PIPELINE_REFERENCE/FLITECAM/>`__.

- General Upkeep

  - Work on the automatic documentation generation and deployment for the
    documentation site.
  - As a consequence of the Photutils 3.0 bump, minimum versions of other
    libraries also have been updated, including Astropy 6.1.4 and NumPy 2.1

1.4.3 (2026-03-25)
==================

- Fix compatibility with NumPy 2.4

  - Numba 0.64 has been released and unlocked NumPy 2.4 for us, revealing some
    issues with the handling of scalars and 0-dimensional arrays.

- New home for the Redux documentation: https://redux.sofiadatacenter.de

  - The documentation has been migrated to a new server hosted by the SOFIA Data Center.
  - The new documentation uses the
    `PyData Sphinx Theme <https://pydata-sphinx-theme.readthedocs.io/en/stable/>`__
    following most of the Scientific Python stack projects.
  - A version selector dropdown now lets you easily switch between documentation
    for different versions of the software.

- FIFI-LS User Manual Iss. 03:

  - Add a description for the Narrow Line Mode (NLM) and redshift
    and rest wavelength override parameters to the documentation.

1.4.2 (2026-01-30)
==================

- Fix compatibility with pandas 3.0

1.4.1 (2026-01-29)
==================
- FIFI-LS User Manual Iss. 02:

  - Document DSI response curves for flux calibration.
  - Additional documentation for new pipeline features, including
    telluric scaling, background scaling & ATRAN interpolation features.

- Support for Python 3.14
- Update Zenodo DOI for citations.

1.4.0 (2025-12-04)
==================
- This is the first full release of SOFIA Redux from the SOFIA Data Center
- Update dependency versions for current support windows

  - Support Python 3.11 to 3.13
  - Update API usage for Astropy, NumPy, SciPy, Matplotlib, Pandas, ...
  - Switch from pyds9 to ds9samp for interaction with DS9
  - Use PyQt6 for the GUI

- Update CI/CD setup for Uni Stuttgart internal GitHub Enterprise
  with local runners
- FIFI-LS pipeline v2.10.0, updates developed at DSI and SDC:

  - Telluric scaling
  - ATRAN spectra interpolation
  - Updated calibration files (to be documented in v1.4.1)
  - New discard logic for pointing errors (not used yet)

1.3.3  (2023-08-06)
===================
- Updated FIFI-LS pipeline to v2.8.0, including separate spatial flats for
  each dichroic.
- Updated calibration files for HAWC, FORCAST, and EXES for prior cycles with
  new analysis from data reprocessing efforts.
- Minor updates for compatibility with third-party dependency changes.
- Update packaging for PEP compliance.

1.3.2 (2023-03-29)
==================
- Minor bug fixes for EXES edge cases (raster flats with non-standard
  overlaps; occasional spurious zero-length orders).
- Updated calibration files for EXES.
- Environment update for astropy management.


1.3.1 (2023-02-11)
==================
- Minor fixes for the spectral viewer, relating to plotting multiple
  spectra with potentially mismatched units.


1.3.0 (2022-12-21)
==================

- Added EXES support (EXES Redux v3.0.0).
- Updated FORCAST pipeline to v2.7.0 to allow the use of measured
  water vapor values in calibration routines.
- Updated FIFI-LS pipeline to v2.7.1, including a bug fix to spectral
  maps rotated to align detector orientation.
- Updated HAWC pipeline to v3.2.0, including improved identification and
  cleaning for discrepant artifacts in scan maps, caused by detector
  level jumps.
- Updated the spectral viewer interface to better control spectra in
  multiple files, apertures, and orders.


1.2.8 (2022-10-25)
==================

- Calibration updates for FORCAST flight series OC9V and HAWC flight
  series OC9W.
- Expanded acceptable formats for reference line lists in the spectral
  viewer.


1.2.7 (2022-09-15)
==================

- Updated FORCAST pipeline to v2.6.0, including additional parameters to
  allow reduction in detector coordinates and a bug fix for water vapor
  optimization in spectroscopy calibration.
- Updated FIFI-LS pipeline to v2.7.0, including improved spatial calibration
  and astrometry fixes for the final spectral cube and an experimental
  scan reduction improvement for OTF data.
- Updated HAWC pipeline to v3.1.0, including a new parameter for the
  scan mapping steps, allowing reduction at alternate grid scales.
- Calibration data updates for FIFI-LS flight series OC9T.
- Added FIFI-LS support to the scan module.
- Added support for more general spectral formats, e.g. text files and simple
  FITS files, to the spectral viewer.


1.2.6 (2022-07-25)
==================

- Calibration data updates for HAWC flight series OC9Q-R.


1.2.5 (2022-06-10)
==================

- Calibration data updates for FORCAST flight series OC9P.
- Add pixel unit display to spectral viewer tool.


1.2.4 (2022-05-25)
==================

- Updated FORCAST pipeline to v2.5.0, including support for NXCAC slit scans.
- Updated FIFI-LS pipeline to v2.6.1, including performance improvements
  and additional memory management for large data sets.
- Calibration data updates for FORCAST and FIFI-LS.


1.2.3 (2022-04-06)
==================

- Updated FORCAST pipeline to v2.4.0, removing scikit-image as a dependency.
- Added a line list overplot feature to the spectral viewer tool.
- Improved log messages for the QAD image viewer tools.
- Completed test coverage and associated minor fixes for the scan module.


1.2.2 (2022-03-15)
==================

- Updated FORCAST and FIFI-LS calibration files.
- Increased test coverage and minor bug fixes for the scan module.


1.2.1 (2022-02-22)
==================

- Various build and test error fixes, including missing scan
  package configuration files.


1.2.0 (2022-02-18)
==================

- Added HAWC+ support (HAWC DRP v3.0.0).
- Updated FIFI-LS bad pixel and flat files.


1.1.1 (2021-12-13)
==================

- Updated FIFI-LS pipeline to v2.6.0, including simplified
  wavelength calibration format.
- Refactored toolkit algorithms for more generalized resampling
  and improved multiprocessing support.


1.1.0 (2021-09-27)
==================

- Added FLITECAM support (FLITECAM Redux v2.0.0) for all observation
  modes.


1.0.0 (2021-07-15)
==================

- Initial release, including FORCAST Redux v2.3.0 and FIFI-LS Redux v2.5.1.
