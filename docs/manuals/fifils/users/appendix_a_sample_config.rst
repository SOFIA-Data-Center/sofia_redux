.. _fifi_appendix_sample_config_file:

Appendix A: Sample parameter config file
========================================

The below is a sample FIFI-LS Redux parameter config file in INI format.
If a value for any parameter is present, it will override the
corresponding default defined by the FIFI-LS reduction object. If not present, the
default value will be used - the sole exception is ``xy_pixel_size``, which is
computed at runtime based on the detector channel of the input data. The parameters
displayed here are the current default values.

.. include:: include/redux_param.cfg
   :literal:

.. raw:: latex

    \clearpage

