.. _config_appendix:

Appendix B: Sample Configuration Files
======================================

Full Pipeline Configuration File
--------------------------------

Below is a copy of the full configuration file used by the pipeline in
the Redux environment (*pipeconf.cfg*). It is in INI format, and is readable
by the configobj Python module.

.. include:: include/pipeconf.cfg
   :literal:

Pipeline Override Configuration File
------------------------------------

Below is a sample override configuration file that demonstrates how to set
override parameters to provide to the HAWC+ pipeline. The parameters
listed here are those most likely to change from one flight series to
another.

.. include:: include/pipeconf_oc8e.cfg
   :literal:

Full Scan Map Configuration File
--------------------------------

Below is a copy of the default global configuration file for the scan
map algorithm. Other configuration files specifying values for specific
instruments or modes may override values in this file.

.. include:: include/default.cfg
   :literal:

HAWC+ Scan Map Configuration File
---------------------------------

Below is the HAWC+ configuration file for the scan map algorithm. Values in
this file override those in the global configuration file for HAWC+ reductions.

.. include:: include/hawc_default.cfg
   :literal:

