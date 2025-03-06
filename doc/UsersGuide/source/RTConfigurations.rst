.. |nbsp| unicode:: 0xA0 
   :trim:

.. role:: raw-html(raw)
    :format: html

.. _Configurations:

*******************************
Regresson Test Configuration
*******************************

The UFS Weather Model (WM) can be run in any of several configurations, from a single-component atmospheric 
model to a fully coupled model with multiple earth system components (e.g., atmosphere, ocean, sea-ice, land, and 
mediator). This chapter explains how forecasts are configured using the regression test (RT) framework. For a full list of 
supported RT configurations, view the `rt.conf <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/rt.conf>`_ file
or visit the :wm-repo:`tests/tests <tree/develop/tests/tests>` directory.

====================
Test Configuration
====================

Each RT test configuration file (located in the ``ufs-weather-model/tests/tests`` 
:wm-repo:`directory <tree/develop/tests/tests>`) 
sets default variables by calling functions from ``ufs-weather-model/tests/default_vars.sh`` 
(view :wm-repo:`default_vars.sh here <blob/develop/tests/default_vars.sh>`). 
Then, the test configuration file sets test-specific variables. These values will override 
the defaults. 

``default_vars.sh`` 
=====================

``default_vars.sh`` sets a series of machine-specific variables. It also contains several functions that set defaults for different types of tests. :numref:`Table %s <def-funcs>` describes what each function does. 

.. _def-funcs:

.. list-table:: ``default_vars.sh`` functions
   :widths: 10 70
   :header-rows: 1
   
   * - Function Name
     - Description
   * - export_fv3_v16
     - Set variables to the FV3 default values for GFS v16 cases. This section will be removed once support for GFSv16 is officially depricated.
   * - export_fv3
     - Set variables to the FV3 default values.
   * - export_tiled
     - Set default values for tiled grid namelist.
   * - export_ugwpv1() 
     - Set default values for the Unified Gravity Wave Drag Physics v1. 
   * - export_cice6()
     - Set default values for the CICE6 model namelist and ``mx100``. 
   * - export_mom6() 
     - Set default values for the MOM6 model namelist and ``mx100``. 
   * - export_ww3()
     - Set default values for the WW3 global model. 
   * - export_fire_behavior()
     - Set default values for the Fire Behavior model. 
   * - export_cmeps()
     - Set default values for the coupled 5-component tests using CMEPS.
   * - export_cpl
     - Set variables to the default values for *coupled* / S2S configurations. 
   * - export_35d_run
     - Set variables to the default values for EMC's weekly coupled benchmark 35d tests (see `rt_35d.conf <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/rt_35d.conf>`__). 
   * - export_datm_cdeps
     - Set variables to the default values for configurations that use the data atmosphere (:term:`DATM`) component. 
   * - export_hafs_datm_cdeps
     - Set variables to the default values for HAFS configurations that use the data atmosphere (DATM) component. 
   * - export_hafs_docn_cdeps
     - Set variables to the default values for HAFS configurations that use the data ocean (:term:`DOCN`) component. 
   * - export_hafs_regional
     - Set variables to the default values for regional HAFS configurations. 
   * - export_hafs
     - Set variables to the default values for HAFS configurations. 
   * - export_hrrr() 
     - Set default values for HRRR test configurations. 
   * - export_hrrr_conus13km()
     - Set default values for hrrr_conus13km test configurations. 
   * - export_rap_common()
     - Set default values that are common to RAP and RRFS v1 test configurations. 
   * - export_rap()
     - Set default values for RAP test configurations. 
   * - export_rrfs_v1()
     - Set default values for RRFS v1 test configurations.
   
Multiple ``default_vars.sh`` functions may be called in a given test, usually starting with the most general function and ending with the most specific. Values set in one function will be overridden when the same values are set in a subsequent function. 

Test Configuration Files
=========================

Individual test configuration files typically start with an ``export TEST_DESCR`` statement describing the test, followed by an ``export CNTL_DIR`` statement indicating the name of the directory that contains the baselines for the experiment. Next, an ``export LIST_FILES`` statement indicates which files the test expects to output from the model run. This list often includes RESTART files. After the LIST_FILES statement, the tests typically call functions from ``default_vars.sh`` to set default values. 

For example, the ``hafs_regional_atm_ocn_wav`` test file lists the files that it will output and then calls three ``export_*`` functions from ``default_vars.sh``, starting in order from the most general to the most specific:

.. code-block:: console

   export LIST_FILES="atmf006.nc \
                   sfcf006.nc \
                   archv.2019_241_06.a \
                   archs.2019_241_06.a \
                   20190829.060000.out_grd.ww3 \
                   20190829.060000.out_pnt.ww3 \
                   ufs.hafs.ww3.r.2019-08-29-21600.nc \
                   ufs.hafs.cpl.r.2019-08-29-21600.nc"

   export_fv3
   export_hafs
   export_hafs_regional

Lastly, the :wm-repo:`test configuration file <blob/develop/tests/tests/hafs_regional_atm_ocn_wav>` sets any test-specific variables for the experiment. These variables will override the default values from ``default_vars.sh``. In the excerpt below, ``...`` indicates omitted lines: 

.. code-block:: console

   export HAFS=true
   export FHMAX=6
   export RESTART_N=${FHMAX}
   export DT_ATMOS=180
   export IDEFLATE=1
   export OUTPUT_FH='3 -1'
   export OUTPUT_FILE="'netcdf' 'netcdf'"
   export SDAY=29
   export SHOUR=00
   export SMONTH=08
   export SYEAR=2019

   ...

   export CDEPS_DOCN=false
   export OCEAN_START_DTG=43340.00000

   export atm_model=fv3
   export ocn_model=hycom
   export wav_model=ww3
   OCN_tasks=60
   WAV_tasks=60
   export coupling_interval_sec=360
   export MESH_ATM=unset

   export FIELD_TABLE=field_table_hafs
   export DIAG_TABLE=diag_table_hafs_template
   export INPUT_NML=input_regional_hafs.nml.IN
   export MODEL_CONFIGURE=model_configure_hafs.IN
   export UFS_CONFIGURE=ufs.configure.hafs_atm_ocn_wav.IN
   export FV3_RUN="hafs_fv3_run.IN hycom_hat10_run.IN hafs_ww3_run.IN"

   if [[ $MACHINE_ID = orion ]]; then
   WLCLK=40
   fi
   ...


Creating New Test Configurations
=================================

Users are welcome to modify current tests for their own use or create new tests to facilitate their own research. 
When creating a test, users will need to add a row for the test in ``rt.conf`` or in their own custom file. 
See :numref:`Section %s <rt.conf>` for more information. 

====================================
Supported Configurations
====================================


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
