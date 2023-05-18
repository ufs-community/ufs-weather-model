.. _Configurations:

*************************
Configurations
*************************

The UFS Weather Model (WM) can be run in any of several configurations, from a single-component atmospheric 
model to a fully coupled model with multiple earth system components (e.g., atmosphere, ocean, sea-ice and 
mediator). This chapter documents a few of the curently supported configurations. For a full list of 
supported configurations, view the `rt.conf <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/rt.conf>`__ 
and `rt.gnu.conf <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/rt_gnu.conf>`__ files. 

.. attention::

   This chapter is a work in progress. There are a multitude of options for configuring the UFS WM, 
   and this chapter merely details a few supported configurations. It will be expanded over time
   to include the full set of configurations supported for WM regression tests (RTs). 

.. _UFS-configurations-documented:

.. list-table:: *Documented UFS Weather Model Configurations*
   :widths: 10 70
   :header-rows: 1
   
   * - Configuration Name
     - Description
   * - :ref:`ATM <atm_documented>`
     - Standalone Atmospheric Model (:term:`ATM`)
   * - :ref:`LND <lnd_documented>`
     - Coupled :term:`CDEPS` - :term:`DATM` - :term:`LND` -:term:`CMEPS`

This chapter details the build and run options for each supported configuration. 
Click on the Configuration Name in :numref:`Table %s <UFS-configurations-documented>` 
to go to that section. Each configuration includes sample code for setting ``CMAKE_FLAGS`` and ``CCPP_SUITES``. 
Additionally, there is a list of preferred physics suites, examples of ``nems.configure`` files, 
and links to information on other input files required to run the model. 

============
Background
============

Each RT configuration file (located in the ``tests/tests`` `directory <https://github.com/ufs-community/ufs-weather-model/tree/develop/tests/tests>`__) sets default variables by calling setup functions from ``tests/default_vars.sh`` (see defaults `here <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/default_vars.sh>`__). Then, the RT configuration file sets test-specific variables; these values will override the defaults. For example, the ``control_c48`` test file sets a list of files that it will use, calls the ``export_fv3`` function from ``default_vars.sh``, and then exports test-specific variables. An excerpt is included below (``...`` indicates omitted lines): 

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

``default_vars.sh`` contains eight functions that set defaults for different types of tests. :numref:`Table %s <def-funcs>` describes what each function does. Multiple functions may be called in a given RT. 

.. _def-funcs:

.. list-table:: *default_vars.sh functions*
   :widths: 10 70
   :header-rows: 1
   
   * - Function Name
     - Description
   * - export_fv3
     - Export variables to the FV3 default values (first common variables, then model specific ones). Different machines may have different defaults for some variables. 
   * - export_cpl
     - Export variables to the default values for *coupled* / S2S configurations. 
   * - export_35d_run
     - Export variables to the default values for EMC's weekly coupled benchmark 35d tests (see `rt_35d.conf <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/rt_35d.conf>`__). 
   * - export_datm_cdeps
     - Export variables to the default values for configurations that use the data atmosphere (:term:`DATM`) component. 
   * - export_hafs_datm_cdeps
     - Export variables to the default values for HAFS configurations that use the data atmosphere (DATM) component. 
   * - export_hafs_docn_cdeps
     - Export variables to the default values for HAFS configurations that use the data ocean (:term:`DOCN`) component. 
   * - export_hafs_regional
     - Export variables to the default values for regional HAFS configurations. 
   * - export_hafs
     - Export variables to the default values for HAFS configurations. 

====================================
Atmospheric Model Configurations
====================================

The atmospheric model configurations all use the UFS WM atmospheric component 
and may couple it with other models (e.g., a wave or aerosol model).

.. _atm_documented:

ATM - Standalone Atmospheric Model
=====================================

The standalone atmospheric model (:term:`ATM`) is an :term:`FV3`-based prognostic 
atmospheric model that can be used for short- and medium-range research and operational 
forecasts. In standalone mode, ``ATM`` is not coupled to any other model. 

Current ATM regression tests cover a wide variety of functionality and involve several physics tests. :numref:`Table %s <atm-rts>` contains a selection of ATM-only RTs; it will be expanded to cover the full range of ATM-only supported configurations in time: 

.. _atm-rts:

.. list-table:: *ATM regression test descriptions*
   :widths: 10 40 10 10 15 5
   :header-rows: 1

   * - Test Name
     - Description
     - Physics suite
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

.. _lnd_documented:

LND
=======

The LND configuration couples :term:`DATM`, :term:`CDEPS`, and :term:`CMEPS` with :term:`LND`. These test use default values set in the ``export_datm_cdeps`` function of ``default_vars.sh``. 

.. _lnd-rts:

.. list-table:: *LND regression test descriptions*
   :widths: 10 40 10 10 15 5
   :header-rows: 1

   * - Test Name
     - Description
     - Physics suite
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





