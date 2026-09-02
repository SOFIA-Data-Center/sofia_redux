Appendix C: Calibration Data Generation
=======================================

The FLITECAM Redux pipeline requires several kinds of auxiliary reference
calibration files, listed in :numref:`flitecam_auxiliary`.  Some of these
are produced by tools packaged with the pipeline.  This section describes the
procedures used to produce these auxiliary files.

.. |ref_wavecal_plots| replace:: :numref:`flitecam_wavecal_plots`

.. |ref_spatcal_plots| replace:: :numref:`flitecam_spatcal_plots`

.. |ref_wavecal_residuals| replace:: :numref:`flitecam_wavecal_residuals`

.. |ref_spatcal_residuals| replace:: :numref:`flitecam_spatcal_residuals`

.. include::  ../../forcast/users/spectral_calibration.rst

.. figure:: images/wavecal_plots.png
   :name: flitecam_wavecal_plots
   :alt: The Redux GUI window, several spectral plot displays with lines marked,
         and a DS9 window showing a spectral image.

   Wavecal mode reduction and diagnostic plots.

.. figure:: images/spatcal_plots.png
   :name: flitecam_spatcal_plots
   :alt: GUI window, spatial profile plot display, and a DS9 window with a spectral image.

   Spatcal mode reduction and diagnostic plots.

.. figure:: images/wavecal_residuals.png
   :name: flitecam_wavecal_residuals
   :alt: An image marked with positions and vertical fit lines and a plot window
         showing fit residuals in X and Y.

   Wavecal mode fit surface and residuals.

.. figure:: images/spatcal_residuals.png
   :name: flitecam_spatcal_residuals
   :alt: An image marked with positions and horizontal fit lines and a plot window
         showing fit residuals in X and Y.

   Spatcal mode fit surface and residuals.
