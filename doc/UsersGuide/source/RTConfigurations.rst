.. role:: raw-html(raw)
    :format: html



Atmospheric Model Configurations
====================================

The atmospheric model configurations all use the UFS WM atmospheric component 
and may couple it with other models (e.g., a wave or aerosol model).
The standalone atmospheric model (:term:`ATM`) is an :term:`FV3`-based prognostic 
atmospheric model that can be used for short- and medium-range research and operational 
forecasts. In standalone mode, ``ATM`` is not coupled to any other model. Current ATM regression tests cover a wide variety of functionality and involve several 
physics tests. 

Input files required for ATM configurations can be viewed in :numref:`Section %s <atm-in>`
or in the `UFS WM RT Data Bucket <https://registry.opendata.aws/noaa-ufs-regtests/>`_. 
Information on ``ufs.configure`` files is available in :numref:`Section %s <ufs-conf>`,
and a sample ATM ``ufs.configure`` file (``ufs.configure.atm.IN``) is available 
`here <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/parm/ufs.configure.atm.IN>`__.

Additional configuration options are listed in :numref:`Table %s <UFS-configurations>`.

Rapid Refresh Forecast System (RRFS)
--------------------------------------

The RRFS configurations use an :term:`ATM`-only configuration on a high-resolution regional grid. 
These tests use the default values set in the ``export_fv3``, ``export_rap_common``, ``export_rrfs_v1``, and/or ``export_hrrr_conus13km`` functions of ``default_vars.sh`` unless other values are explicitly set in a given test file. In all tests, the values in ``export_fv3`` are set first. Depending on the test, some of these values may be overriden by ``export_rrfs_v1`` (which includes values from ``export_rap_common``) or ``export_hrrr_conus13km``. 

.. note:: 

   ``export_rrfs_v1`` calls ``export_rap_common``, which calls ``export_fv3``. Values from ``export_fv3`` are set first, followed by values in ``export_rap_common`` and then values in ``export_rrfs_v1``. Values in italics indicate that the value is inherited from a previously-called function. 

Input files required for RRFS ATM configurations can be downloaded from the `UFS WM RT Data Bucket <https://registry.opendata.aws/noaa-ufs-regtests/>`_. Users who wish to run additional (unsupported) cases may also find useful data in the `NOAA RRFS data bucket <https://registry.opendata.aws/noaa-rrfs/>`_. 

Information on ``ufs.configure`` files is available in :numref:`Section %s <ufs-conf>`. The supported RRFS WM RTs use the same ``ufs.configure`` file that ATM-only tests do (``ufs.configure.atm.IN``). This file can be viewed in the :wm-repo:`ufs-weather-model/tests/parm directory <tree/develop/tests/parm>`. Additionally, users can find examples of various RRFS configuration files in the ``ufs-weather-model/tests/parm`` directory. These files include ``model_configure_*``, ``*_run.IN`` (input run), ``*.nml.IN`` (input namelist), ``field_table_*``, and ``diag_table_*`` files.
