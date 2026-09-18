.. _fifi_appendix_fluxcal_history:

Appendix C: Flux Calibration History
====================================

As described in the section :ref:`fifi_subsection_reduction_algorithms`, the flux calibration
response curves are generated from comparisons of FIFI-LS Mars Observations
with a time-resolved Mars flux model. Over the years, this procedure has been
refined and reprocessed using newer pipeline versions, resulting in a number
of different versions of the response curves. Here, the three major iterations
of flux calibration generation are described.

"USRA" Response Curves
----------------------
Fadda *et al*. 2023 [#Fadda2023]_ describes in detail the first major update to
the flux calibration procedure since the release of Redux (details of prior calibrations
are spurious and will not be documented here). This update was developed and implemented
by the Universities Space Research Association (USRA) in 2023, and was used for FIFI-LS
Redux versions 2.0.0 through 2.9.0.

This approach used mostly Mars as a calibrator, supplemented in the blue channel by a small
number of Callisto and Uranus measurements. The full dataset spanned 77 different flight-channel
configurations, from October 2015 to August 2022, and thus a considerable range of atmospheric
conditions and calibrator separations. The datasets were reduced to WXY products, using mostly standard 
Redux parameters. From the resulting data cube(s), a wavelength-dependent flux
was computed from a summation of the flux in all spaxels within each wavelength slice, excluding
the partially illuminated spaxel stack. These fluxes were then compared to the modelled
calibrator fluxes, including the model of Lellouch & Amri for Mars [#fn_fifi_mars2]_, and the
ESA2 models for Uranus and Callisto [#fn_fifi_uranus]_, [#fn_fifi_callisto]_. Resulting
plots of response datapoints and curves are shown in :numref:`fifi_response_fadda`.

There exist three noteworthy aspects to this approach. First, as outlined in Table 13 of the
FIFI-LS Handbook for Archive Users [#fn_fifi_archiveman]_, the telluric correction performed
on the reduced Mars data include flight-wise scaling factors used to convert derived ECMWF
water vapor values to FIFI-LS values. Such factors, and the flightwise variation between them,
introduce significant systematic variation and potential errors into the derived flux response,
and can be circumvented with the use of the new on-the-fly ECMWF-based telluric correction
developed at the SDC.

Second, as outlined in Tables 5 & 6 of the Archive Handbook, additional flightwise factors
were introduced to scale the flightwise response curves to a common median. As displayed in
:numref:`fifi_fadda_fluxcal_factors`, this scaling is quite successful in minimising the
significant flightwise scatter in the derived response datapoints. However, similar to the
telluric correction factors, these factors introduce further systematic variation to the
flux calibration process. At this time, the origin of these residual variations between flights
can possibly be understood as responsivity changes of the detector system that have not yet been
characterized, or related to any other systematic effect. Therefore, it is preferable to
avoid using such factors, until the cause of the scatter has been identified.

Third, the data reduction process utilises the standard Redux resample algorithm, which includes
a 3-dimensional interpolation onto an even spatial-spectral grid. Spatial interpolation of the
highly undersampled Mars disk may introduce significant errors in the derived fluxes when resampling
onto the output grid, particularly at low Mars separations, where a significant fraction of the Mars
flux may fall outside of the FIFI-LS FOV. While unlikely to be the cause of the significant scatter
seen in the response datapoints, this systematic effect may nevertheless be circumvented with the use
of a Drizzle-based resampling algorithm, which preserves the total flux between the input and
output grids.

.. [#Fadda2023] \D. Fadda *et al*. 2023 ApJ 166 237, https://doi.org/10.3847/1538-3881/acffb4
.. [#fn_fifi_mars2] See https://lira.obspm.fr/perso/emmanuel-lellouch/mars/index.php
.. [#fn_fifi_uranus] \G. Orton *et al*. 2014 Icarus 243 471-493, https://doi.org/10.1016/j.icarus.2014.07.012
.. [#fn_fifi_callisto] \T. G. Müller *et al*. 2016 A&A 588 A109, https://doi.org/10.1051/0004-6361/201527371
.. [#fn_fifi_archiveman] See https://irsa.ipac.caltech.edu/data/SOFIA/docs/instruments/handbooks/FIFI-LS_Handbook_for_Archive_Users_Ver1.0.pdf

.. figure:: images/fifi_response_fadda.png
   :alt: Fadda response curves
   :name: fifi_response_fadda

   Response curve fits derived from USRA flux calibration, overplotted on telluric-corrected
   data divided by modelled calibrator fluxes, for all channel/dichoric/order configurations.
   Each datapoint color corresponds to a different FIFI-LS flight. Upper 2 panels are for
   post-filter change data, lower 2 panels for pre-filter change data (also indicated by the dates on the plots).
   Fitted response curves indicated by colored lines. From Fadda *et al*. 2023, Fig. 11.

.. figure:: images/fifi_fadda_fluxcal_factors.png
   :alt: Fadda flux calibration factors
   :name: fifi_fadda_fluxcal_factors

   Response datapoints derived from USRA flux calibration, for the post filter-change RD105
   configuration. Each color corresponds to a different FIFI-LS flight. Top panel displays
   original response datapoints. Middle panel displays datapoints after scaling with flightwise
   factors, to a median response curve, displayed in black. Bottom panel displays residual in
   scaled responses from the median after scaling. From Fadda *et al*. 2023, Fig. 11.


"DSI" Response Curves
---------------------
From v2.10.0 onwards, the pipeline used response curves generated at the
Deutsches SOFIA Institut (DSI) in 2024. This method utilises only Mars datasets,
reducing them similarly to the USRA method. Then, for each wavelength slice, it identifies the centre of
Mars via the brightest pixel, and sums all flux within a static aperture around that pixel. This approach
better captures the true measured flux, albeit without correcting for flux lost outside of the
FIFI-LS field of view. Additionally, flightwise atmospheric calibration factors
that were applied to the prior calibration have been removed in post-filter change
response curves. The pre-filter change response curves remain the same as before,
except those for the B2 filters, which have been scaled from their post-filter change
counterparts. This corrects for artefacts in the response curves, and also extends
response coverage for the [OIII] 52 micron line in the old filter configuration.
:numref:`fifi_fadda_dsi_response_comparison` displays the difference between the USRA
and DSI response versions for all FIFI-LS filter configurations.

.. figure:: images/fifi_fadda_dsi_response_comparison.png
   :alt: Fadda DSI response comparison
   :name: fifi_fadda_dsi_response_comparison
   :height: 800

   Comparison of Fadda et al. 2023 (orange) and DSI 2024 (blue) response curves
   for all FIFI-LS filter configurations.


"SDC" Response Curves
---------------------
Based on issues identified in both the USRA and DSI flux calibration approaches, an
updated flux calibration was developed at the SOFIA Data Center (SDC) in 2026.
This method utilises the aforementioned Drizzle-based resampling algorithm, along with
a fully time-resolved Mars model, and input Mars datasets reduced using newly developed
Redux features, such as the on-the-fly telluric correction and updated spectral resolution
tables.

During development, much effort was made to attempt to identify and mitigate
sources of scatter in the response datapoints, both between flights, and within individual
flights, without the use of flightwise scaling factors. As shown in :numref:`fifi_response_scatter`,
the scatter largely remains, however the resulting response curves indeed show some deviation
from the previous two approaches. Nevertheless, thanks to the use of flux-preserving spatial-spectral
resampling and the updated Mars model, these SDC response curves are expected to be
more robust than prior versions. Remaining response scatter is propagated into the Redux
pipeline in the form of a wavelength-varying systematic error (extension ``CAL_ERROR``), stored separately to the existing
statistical error. Prior versions stored a far more generic systematic error as a
single value in the ``CALERR`` header keyword.

.. figure:: images/fifi_response_scatter.png
   :alt: SDC response scatter
   :name: fifi_response_scatter

   SDC-derived response datapoints for RD105 new filter configuration, with fitted response curve.
   Considerable scatter remains in the datapoints both between flights and within individual flights.

Response values for the B2 D130 old filter configuration below 56.5 µm are scaled from the previously
used DSI response curve, which itself was scaled from their post-filter change counterparts, due to 
the lack of calibration data in this region as described in the "DSI" response curve section. The
response errors in this region are linearly extrapolated from the existing errors derived
from SDC-computed response datapoints. Therefore, calibration error maps covering wide wavelength
ranges within this region will exhibit a linear increase in error towards shorter wavelengths.

.. figure:: images/fifi_response_scaling.png
   :alt: SDC response scaling
   :name: fifi_response_scaling

   SDC-derived response curve for the B2 D130 old filter configuration, consisting of a fitted curve 
   for data above 56.5 µm and a scaled response curve for data below 56.5 µm.

Order Filter Keywords
---------------------
As of FIFI-LS Redux v2.0.0., the pipeline accomodates an additional blue order keyword,
``G_FLT_B``, alongside the existing ``G_ORD_B`` keyword. For the vast majority of raw datasets
produced after this update, the two keywords are identical. The remaining datasets are all of
flux calibrators, and the additional keyword allows the pipeline to apply spectral flatfields
and response curves from one filter to the data taken with the other filter, necessary at
some specific wavelengths.
