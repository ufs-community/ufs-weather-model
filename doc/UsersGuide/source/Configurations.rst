.. |nbsp| unicode:: 0xA0 
   :trim:

.. role:: raw-html(raw)
    :format: html

.. _Configurations:

*************************
Configurations
*************************

The UFS Weather Model (WM) can be run in any of several configurations, from a single-component atmospheric 
model to a fully coupled model with multiple earth system components (e.g., atmosphere, ocean, sea-ice, land, and 
mediator). This chapter documents a few of the currently supported configurations. For a full list of 
supported configurations, view the `rt.conf <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/rt.conf>`__ file. 

.. attention::

   This chapter is a work in progress. There are a multitude of options for configuring the UFS WM, 
   and this chapter merely details a few supported configurations. It will be expanded over time
   to include the full set of configurations supported for WM regression tests (RTs). 

.. _UFS-configurations-documented:

.. list-table:: *Documented UFS Weather Model Configuration Categories*
   :widths: 10 70
   :header-rows: 1
   
   * - Configuration Category
     - Description
   * - :ref:`ATM <atm-documented>`
     - Standalone Atmospheric Model (:term:`ATM`)
   * - :ref:`ATML <atml-documented>`
     - Coupled :term:`ATM` and :term:`LND`
   * - :ref:`LND <lnd-documented>`
     - Coupled :term:`CDEPS` - :term:`DATM` - :term:`LND` -:term:`CMEPS`
   * - :ref:`RRFS <rrfs-documented>`
     - :term:`ATM` with :term:`data assimilation`

This chapter details the supported build/run options for each supported configuration. 
Click on the configuration category in :numref:`Table %s <UFS-configurations-documented>` 
to go to that section. Each configuration category includes sample code for setting ``CMAKE_FLAGS`` and ``CCPP_SUITES``. 
Additionally, there is a list of preferred physics suites, examples of ``ufs.configure`` files, 
and links to information on other input files required to run the model. 

============
Background
============

Each RT configuration file (located in the ``ufs-weather-model/tests/tests`` 
`directory <https://github.com/ufs-community/ufs-weather-model/tree/develop/tests/tests>`__) 
sets default variables by calling setup functions from ``ufs-weather-model/tests/default_vars.sh`` 
(see defaults `here <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/default_vars.sh>`__). 
Then, the RT configuration file sets test-specific variablesthese values will override 
the defaults. For example, the ``control_c48`` test file sets a list of files that 
it will use, calls the ``export_fv3`` function from ``default_vars.sh``, and then exports 
test-specific variables. An excerpt is included below (``...`` indicates omitted lines): 

.. code-block:: console

   export LIST_FILES="sfcf000.nc \
                   sfcf024.nc \
                   atmf000.nc \
                   atmf024.nc \
                   RESTART/20210323.060000.coupler.res \
                   RESTART/20210323.060000.fv_core.res.nc \
                   ...
                   RESTART/20210323.060000.sfc_data.tile5.nc \
                   RESTART/20210323.060000.sfc_data.tile6.nc"

   export_fv3

   export INPES=1
   export JNPES=1
   export WRTTASK_PER_GROUP=2
   export NPZ=127
   export NPZP=128
   export NPX=49
   export NPY=49
   export DT_ATMOS=1200
   ...

``default_vars.sh`` contains eight functions that set defaults for different types of tests. :numref:`Table %s <def-funcs>` describes what each function does. 

.. _def-funcs:

.. list-table:: *default_vars.sh functions*
   :widths: 10 70
   :header-rows: 1
   
   * - Function Name
     - Description
   * - export_fv3
     - Set variables to the FV3 default values (first common variables, then model-specific ones). Different machines may have different defaults for some variables. 
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

Multiple ``default_vars.sh`` functions may be called in a given test. Values set in one
function will be overridden when the same values are set in a subsequent function. 

The most up-to-date list of ``develop`` branch data required for each test is available in 
the `UFS WM RT Data Bucket <https://registry.opendata.aws/noaa-ufs-regtests/>`__.
Users should click on "Browse Bucket" and navigate to the most recent date (in ``develop-YYYY-MM-DD`` format).
Then, users should select *Intel* or *GNU* based on the compiler used in the test they 
want to run and then select the test name to see the required data. 

====================================
Atmospheric Model Configurations
====================================

The atmospheric model configurations all use the UFS WM atmospheric component 
and may couple it with other models (e.g., a wave or aerosol model).

.. _atm-documented:

ATM - Standalone Atmospheric Model
=====================================

The standalone atmospheric model (:term:`ATM`) is an :term:`FV3`-based prognostic 
atmospheric model that can be used for short- and medium-range research and operational 
forecasts. In standalone mode, ``ATM`` is not coupled to any other model. 

Current ATM regression tests cover a wide variety of functionality and involve several 
physics tests. :numref:`Table %s <atm-rts>` contains a small selection of ATM-only RTs; 
it will be expanded to cover the full range of ATM-only supported configurations in time: 

.. _atm-rts:

.. list-table:: *ATM regression test descriptions*
   :widths: 10 40 10 10 15 5
   :header-rows: 1

   * - Test Name
     - Description
     - Physics Suite (see `namelist options <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_c_c_p_psuite_nml_desp.html>`__)
     - DT_ATMOS
     - Start Date
     - Forecast Length (hours)
   * - `control_c48 <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/control_c48>`__
     - Compare global control C48L127 results with previous trunk version
     - FV3_GFS_v16
     - 1200
     - 2021-03-22 06:00:00
     - 24
   * - `control_p8 <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/control_p8>`__
     - Compare global control results with previous trunk version
     - FV3_GFS_v17_p8
     - 720
     - 2021-03-22 06:00:00
     - 24
   * - `regional_control <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/regional_control>`__
     - FV3 regional control (hi-res 3km, small domain) test
     - FV3_GFS_v15_thompson_mynn_lam3km
     - 1800
     - 2016-10-03 00:00:00
     - 6

**Sample** ``CMAKE_FLAGS`` **Setting**

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATM -DCCPP_SUITES=FV3_GFS_v16,FV3_GFS_v17_p8,FV3_GFS_v15_thompson_mynn_lam3km -D32BIT=ON"

**Supported Physics Suites**

.. list-table:: *Physics suites used in the ATM configurations above*
   :widths: 10 50
   :header-rows: 1

   * - Physics Suite
     - Description
   * - FV3_GFS_v16
     - The :term:`CCPP` GFS_v16 physics suite is described in the CCPP documentation `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_g_f_s_v16_page.html>`__.
   * - FV3_GFS_v17_p8
     - The CCPP GFS_v17_p8 physics suite is described in the CCPP documentation `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_g_f_s_v17_p8_page.html>`__. 
   * - FV3_GFS_v15_thompson_mynn_lam3km
     - The CCPP GFS_v15 physics suite with the Thompson Aerosol-Aware Cloud Microphysics Scheme 
       (see `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_t_h_o_m_p_s_o_n.html>`__) and 
       Mynn Surface Layer Module (see `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/group__mynn__sfc.html>`__) 
       tailored for a limited area model (LAM) 3-km resolution grid.

**Additional Information**

Input files required for ATM configurations can be viewed in :numref:`Section %s <atm-in>`
or in the `UFS WM RT Data Bucket <https://registry.opendata.aws/noaa-ufs-regtests/>`__. 
Information on ``ufs.configure`` files is available in :numref:`Section %s <ufs-conf>`,
and a sample ATM ``ufs.configure`` file (``ufs.configure.atm.IN``) is available 
`here <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/parm/ufs.configure.atm.IN>`__.


ATMW
=======

**COMING SOON!**

ATMAERO
=========

**COMING SOON!**

ATMAQ
=======

**COMING SOON!**

.. _atml-documented:

ATML
======

The ATML configuration couples :term:`ATM` with :term:`LND`. 
These tests use default values set in the ``export_fv3`` function of ``default_vars.sh``. 

.. attention::
   There is an issue with ``-D32BIT=ON`` in the ATM-LND tests, and NoahMP requires r8 libraries.

.. COMMENT: Should "r8" be "p8"?

.. _atml-rts:

.. list-table:: *ATML regression test descriptions*
   :widths: 10 40 10 10 15 5
   :header-rows: 1

   * - Test Name
     - Description
     - Physics Suite (see `namelist options <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_c_c_p_psuite_nml_desp.html>`__)
     - DT_ATMOS
     - Start Date
     - Forecast Length (hours)
   * - control_p8_atmlnd_sbs
     - Compare global control results with previous trunk version
     - FV3_GFS_v17_p8
     - 720
     - 2021-03-22 06:00:00
     - 24

**Sample** ``CMAKE_FLAGS`` **Setting**

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATML -DCCPP_SUITES=FV3_GFS_v17_p8"


**Supported Physics Suites**

.. list-table:: *Physics suites used in the ATM configurations above*
   :widths: 10 50
   :header-rows: 1

   * - Physics Suite
     - Description
   * - FV3_GFS_v17_p8
     - The :term:`CCPP` GFS_v17_p8 physics suite is described in the CCPP documentation `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_g_f_s_v17_p8_page.html>`__. 

**Additional Information**

Input files required for ATML configurations can be viewed in :numref:`Section %s (ATM) <atm-in>` 
and :numref:`Section %s (LND) <lnd-in>` or in the `UFS WM RT Data Bucket <https://registry.opendata.aws/noaa-ufs-regtests/>`__. 
Information on ``ufs.configure`` files is available in :numref:`Section %s <ufs-conf>`,
and a sample ATML ``ufs.configure`` file (``ufs.configure.atm_lnd.IN``) is available 
`here <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/parm/ufs.configure.atm_lnd.IN>`__.


.. _rrfs-documented:

=======================================
Rapid Refresh Forecast System (RRFS)
=======================================

The RRFS configurations use an :term:`ATM`-only configuration on a high-resolution 
regional grid with data assimilation capabilities. 
These tests use the default values set in the ``export_fv3``, ``export_rap_common``, ``export_rrfs_v1``, and/or ``export_hrrr_conus13km`` functions of ``default_vars.sh`` unless other values are explicitly set in a given test file. In all tests, the values in ``export_fv3`` are set first. Depending on the test, some of these values may be overriden by ``export_rrfs_v1`` (which includes values from ``export_rap_common``) or ``export_hrrr_conus13km``. :numref:`Table %s <rrfs-default-vars-comparison>` compares the values set in ``export_fv3`` to the values set in the other functions. 

.. note:: 

   ``export_rrfs_v1`` calls ``export_rap_common``, which calls ``export_fv3``. Values from ``export_fv3`` are set first, followed by values in ``export_rap_common`` and then values in ``export_rrfs_v1``. Values in italics indicate that the value is inherited from a previously-called function. 

.. _rrfs-default-vars-comparison:

.. csv-table:: *RRFS Default Variables*
   :file: tables/RRFSDefaultVariables.csv
   :widths: 50 10 10 10 10
   :header-rows: 1
   :stub-columns: 1

Current RRFS regression tests cover a wide variety of functionality and involve several 
physics tests. :numref:`Table %s <rrfs-rts>` (below) contains a selection of RTs for RRFS functionality. Blanks indicate that the value comes from the default setting file. These default values are listed in :numref:`Table %s <rrfs-default-vars-comparison>` above. 

.. _rrfs-rts:

.. csv-table:: *RRFS regression test descriptions*
   :file: tables/rrfs-rts.csv
   :widths: 20 20 30 50 10 10 10
   :header-rows: 1

**Sample** ``CMAKE_FLAGS`` **Setting**

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATM -DCCPP_SUITES=FV3_RAP,FV3_HRRR,FV3_RRFS_v1beta,FV3_RRFS_v1nssl -D32BIT=ON"

**Supported Physics Suites**

.. list-table:: *Physics suites used in the RRFS configurations above*
   :widths: 10 50
   :header-rows: 1

   * - Physics Suite
     - Description
   * - FV3_HRRR
     - The FV3_HRRR physics suite is described in the :term:`CCPP` documentation `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_h_r_r_r_suite_page.html>`__.
   * - FV3_RRFS_v1beta 
     - The FV3_RRFS_v1beta physics suite is described in the CCPP documentation `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_r_r_f_s_v1beta_page.html>`__.
   * - FV3_RRFS_v1nssl
     - The FV3_RRFS_v1nssl physics suite is similar to the *FV3_RRFS_v1beta* suite; however, it uses the NSSL 2-moment microphysics scheme instead of the Thompson microphysics scheme.


**Additional Information**

Each test file lists the input files required for a given test. Input files required for RRFS ATM configurations can be downloaded from the `UFS WM RT Data Bucket <https://registry.opendata.aws/noaa-ufs-regtests/>`__. Users who wish to run additional (unsupported) cases may also find useful data in the `NOAA RRFS data bucket <https://registry.opendata.aws/noaa-rrfs/>`__. 

Information on ``ufs.configure`` files is available in :numref:`Section %s <ufs-conf>`. The supported RRFS WM RTs use the same ``ufs.configure`` file that ATM-only tests do (``ufs.configure.atm.IN``). This file can be viewed in the ``ufs-weather-model/tests/parm`` `directory <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/parm/ufs.configure.atm.IN>`__. 

Additionally, users can find examples of various RRFS configuration files in the ``ufs-weather-model/tests/parm`` `directory <https://github.com/ufs-community/ufs-weather-model/tree/develop/tests/parm>`__. These files include ``model_configure_*``, ``*_run.IN`` (input run), ``*.nml.IN`` (input namelist), ``field_table_*``, and ``diag_table_*`` files.  

.. _lnd-documented:

=======
LND
=======

The LND configuration couples :term:`DATM`, :term:`CDEPS`, and :term:`CMEPS` with :term:`LND`. These tests use default values set in the ``export_datm_cdeps`` function of ``default_vars.sh``. 

.. _lnd-rts:

.. list-table:: *LND regression test descriptions*
   :widths: 10 40 10 10 15 5
   :header-rows: 1

   * - Test Name
     - Description
     - Physics Suite
     - DT_ATMOS
     - Start Date
     - Forecast Length (hours)
   * - datm_cdeps_lnd_gswp3
     - DATM_CDEPS_NOAHMP_GSWP3 - control
     - N/A
     - N/A
     - 2000-01-01 00:00:00
     - 24
   * - datm_cdeps_lnd_gswp3_rst
     - DATM_CDEPS_NOAHMP_GSWP3_RST - control restart
     - N/A
     - N/A
     - 2000-01-01 12:00:00
     - 12

**Sample** ``CMAKE_FLAGS`` **Setting**

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=LND"

**Additional Information**

Input files required for LND configurations can be viewed in :numref:`Section %s (LND) <lnd-in>` 
or in the `UFS WM RT Data Bucket <https://registry.opendata.aws/noaa-ufs-regtests/>`__. 
Information on ``ufs.configure`` files is available in :numref:`Section %s <ufs-conf>`,
and a sample ATML ``ufs.configure`` file (``ufs.configure.atm_lnd.IN``) is available 
`here <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/parm/ufs.configure.atm_lnd.IN>`__.


=============================================
Seasonal to Subseasonal (S2S) Configurations
=============================================

**COMING SOON!**

==============
NG-GODAS
==============

**COMING SOON!**

========================================================
Hurricane Analysis and Reforecast System Configurations
========================================================

**COMING SOON!**





