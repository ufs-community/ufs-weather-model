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
   * - :ref:`ATMW <atmw-documented>`
     - Coupled :term:`ATM` and :term:`WW3`
   * - :ref:`ATMAERO <atmaero-documented>`
     - Coupled :term:`ATM` and :term:`GOCART`
   * - :ref:`ATML <atml-documented>`
     - Coupled :term:`ATM` and :term:`LND`
   * - :ref:`LND <lnd-documented>`
     - Coupled :term:`CDEPS` - :term:`DATM` - :term:`LND` -:term:`CMEPS`
   * - :ref:`RRFS <rrfs-documented>`
     - :term:`ATM` with :term:`data assimilation`
   * - :ref:`HAFS <hafs-documented>`
     - Coupled components may include :term:`CDEPS` - :term:`ATM` - :term:`HYCOM` - :term:`WW3` - :term:`MOM6` - :term:`CMEPS`

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

.. _atmw-documented:

ATMW
=======

The ATMW configuration couples :term:`ATM` with :term:`WaveWatch III`. 
These tests use default values set in the ``export_fv3`` function of ``default_vars.sh``.

.. list-table:: *ATMW regression test descriptions*
   :widths: 50 10 30 50 10 10 10 10 10
   :header-rows: 1

   * - Test |nbsp| Name |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp|
     - Description
     - General Physics Parameters
     - Detailed |nbsp| Physics |nbsp| Parameters |nbsp| (see |nbsp| namelist |nbsp| options `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_c_c_p_psuite_nml_desp.html>`__ |nbsp| for variable definitions)
     - Start |nbsp| Date |nbsp| |nbsp| |nbsp| |nbsp|
     - Fcst Length (hours)
     - Output Grid
     - Configuration Files
     - Other
   * - `atmwav_control_noaero_p8 <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/atmwav_control_noaero_p8>`__
     - Compare global control results with previous trunk version
     - **Suite:** CCPP_SUITE="FV3_GFS_v16" :raw-html:`<br/> <br/>`

       **Microphysics:** IMP_PHYSICS=8 :raw-html:`<br/> <br/>`

       **Time Step:** DT_ATMOS=720 :raw-html:`<br/> <br/>`
     
     - **Set to FALSE:** LHEATSTRG, DO_UGWP_V1, DO_GSL_DRAG_LS_BL, DO_GSL_DRAG_TOFD, DO_UGWP_V1_OROG_ONLY, DO_UGWP_V0_NST_ONLY, LDIAG_UGWP, CA_GLOBAL, LANDICE, LGFDLMPRAD, DO_SAT_ADJ, MULTIGRID, USE_CICE_ALB, DO_RRTMGP :raw-html:`<br/> <br/>`
       **Set to TRUE:** USE_MERRA2, LSEASPRAY, DO_UGWP_V0, DO_GSL_DRAG_SS, DO_CA, CA_SGS, CA_TRIGGER, TILEDFIX, CPL, CPLWAV, CPLWAV2ATM, FRAC_GRID, WRITE_NSFLIP, DOGP_CLDOPTICS_LUT, DOGP_LWSCAT, DOGP_SGS_CNV, SATMEDMF :raw-html:`<br/> <br/>`
       **Set to VALUE:** IALB=2, IEMS=2, LSM=2, IOPT_DVEG=4, IOPT_CRS=2, IOPT_RAD=3, IOPT_ALB=1, IOPT_STC=3, IOPT_SFC=3, IOPT_TRS=2, IOPT_DIAG=2, D2_BG_K1=0.20, D2_BG_K2=0.04, PSM_BC=1, DDDMP=0.1, IAER=1011, KNOB_UGWP_VERSION=0, KNOB_UGWP_NSLOPE=1, NCA=1, NCELLS=5, NLIVES=12, NTHRESH=18, NSEED=1, NFRACSEED=0.5, NSPINUP=1, ISEED_CA=12345, FSICL=0, FSICS=0, DNATS=0,  DZ_MIN=6, cap_dbug_flag=0, MIN_SEAICE=0.15, 
     - 2021-03-22 06:00:00
     - 12
     - OUTPUT_GRID=gaussian_grid :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=$INPES_cpl_atmw, JNPES=$JNPES_cpl_atmw, NPZ=127, NPZP=128
     - FIELD_TABLE=field_table_thompson_noaero_tke
       DIAG_TABLE=diag_table_p8_template
       INPUT_NML=cpld_control.nml.IN
       UFS_CONFIGURE=ufs.configure.atmw.IN
       FV3_RUN=control_run.IN
     - RUNTYPE=startup, med_model=cmeps, atm_model=fv3, wav_model=ww3

.. _atmaero-documented:

ATMAERO
=========

The ATMAERO configuration couples :term:`ATM` with :term:`GOCART`. 
These tests use default values set in the ``export_fv3`` function of ``default_vars.sh``.

.. attention:: 
   
   Certain physics-related settings are common to all of the supported RRFS configurations. These values are set in each test's configuration file because they differ from the ``default_vars.sh`` values:

      General Physics Parameters:
          * **Suite:** CCPP_SUITE= `FV3_GFS_v17_p8 <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_g_f_s_v17_p8_page.html>`__
          * **Microphysics:** IMP_PHYSICS=8
          * **Time Step:** DT_ATMOS=720

      Detailed Physics Parameters:
          * **Set to FALSE:** DO_UGWP_V1, DO_GSL_DRAG_LS_BL, DO_GSL_DRAG_TOFD, DO_UGWP_V1_OROG_ONLY, DO_UGWP_V0_NST_ONLY, LDIAG_UGWP, CA_GLOBAL, LANDICE, LGFDLMPRAD, DO_SAT_ADJ, USE_CICE_ALB, DO_RRTMGP
          * **Set to TRUE:** WRITE_DOPOST, CPL, CPLCHM, USE_MERRA2, LSEASPRAY, DO_UGWP_V0, DO_GSL_DRAG_SS, DO_CA, CA_SGS, CA_TRIGGER, TILEDFIX, FRAC_GRID, WRITE_NSFLIP, DOGP_CLDOPTICS_LUT, DOGP_LWSCAT, DOGP_SGS_CNV, SATMEDMF
          * **Set to VALUE:** NSTF_NAME='2,0,0,0,0', atm_model='fv3', chm_model='gocart', DOMAINS_STACK_SIZE=8000000, IALB=2, IEMS=2, LSM=2, IOPT_DVEG=4, IOPT_CRS=2, IOPT_RAD=3, IOPT_ALB=1, IOPT_STC=3, IOPT_SFC=3, IOPT_TRS=2, IOPT_DIAG=2, D2_BG_K1=0.20, D2_BG_K2=0.04, PSM_BC=1, DDDMP=0.1, GWD_OPT=2, KNOB_UGWP_VERSION=0, KNOB_UGWP_NSLOPE=1, NCA=1, NCELLS=5, NLIVES=12, NTHRESH=18, NSEED=1, NFRACSEED=0.5, NSPINUP=1, ISEED_CA=12345, FSICL=0, FSICS=0, DZ_MIN=6, MIN_SEAICE=0.15
   
   The "Detailed Physics Parameters" column in :numref:`Table %s <atmaero-rts>` details physics settings that differ from both the ``default_vars.sh`` values and these ATMAERO-specific defaults. 

.. _atmaero-rts:

.. list-table:: *ATMAERO regression test descriptions*
   :widths: 50 10 50 10 10 10 10 10
   :header-rows: 1

   * - Test |nbsp| Name |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp|
     - Description
     - Detailed |nbsp| Physics |nbsp| Parameters |nbsp| (see |nbsp| namelist |nbsp| options `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_c_c_p_psuite_nml_desp.html>`__ |nbsp| for variable definitions)
     - Start |nbsp| Date |nbsp| |nbsp| |nbsp| |nbsp|
     - Fcst Length (hours)
     - Output Grid
     - Configuration Files
     - Other
   * - `atmaero_control_p8 <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/atmaero_control_p8>`__
     - Compare global results for prognostic aerosols with previous trunk version
     - **Set to FALSE:** LHEATSTRG :raw-html:`<br/> <br/>`
       **Set to TRUE:** ATMAERO default values only :raw-html:`<br/> <br/>`
       **Set to VALUE:** IAER=1011, DNATS=2
     - 2021-03-22 06:00:00
     - 24
     - OUTPUT_GRID=gaussian_grid :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=${INPES_atmaero}, JNPES=${JNPES_atmaero}, NPZ=127, NPZP=128
     - FIELD_TABLE=field_table_thompson_noaero_tke_GOCART
       DIAG_TABLE=diag_table_cpld.IN
       INPUT_NML=ufs.configure.atmaero_esmf.IN
       UFS_CONFIGURE=ufs.configure.atmaero.IN
       FV3_RUN=control_run.IN
     - RESTART_INTERVAL=12 -1
   * - `atmaero_control_p8_rad <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/atmaero_control_p8_rad>`__
     - Compare global results for prognostic aerosols with previous trunk version
     - **Set to FALSE:** ATMAERO values only :raw-html:`<br/> <br/>`
       **Set to TRUE:** LHEATSTRG :raw-html:`<br/> <br/>`
       **Set to VALUE:** IAER=2011, DNATS=2
     - 2021-03-22 06:00:00
     - 24
     - OUTPUT_GRID=gaussian_grid :raw-html:`<br/> <br/>`
       **Grid Parameters**: NPZ=127, NPZP=128
     - FIELD_TABLE=field_table_thompson_noaero_tke_GOCART
       DIAG_TABLE=diag_table_cpld.IN
       INPUT_NML=cpld_control.nml.IN
       UFS_CONFIGURE=ufs.configure.atmaero_esmf.IN
       FV3_RUN=control_run.IN
     - RESTART_INTERVAL=12 -1
   * - `atmaero_control_p8_rad_micro <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/atmaero_control_p8_rad_micro>`__
     - Compare global results for prognostic aerosols with previous trunk version
     - **Set to FALSE:**  :raw-html:`<br/> <br/>`
       **Set to TRUE:** LHEATSTRG :raw-html:`<br/> <br/>`
       **Set to VALUE:** IAER=2011, DNATS=4
     - 2021-03-22 06:00:00
     - 24
     - OUTPUT_GRID=gaussian_grid :raw-html:`<br/> <br/>`
       **Grid Parameters**: NPZ=127, NPZP=128
     - FIELD_TABLE=field_table_thompson_noaero_tke_GOCART
       DIAG_TABLE=diag_table_p8_gocart_micro
       INPUT_NML=merra2_thompson.nml.IN
       UFS_CONFIGURE=ufs.configure.atmaero_esmf.IN
       FV3_RUN=control_run.IN
     - RESTART_INTERVAL='12 -1'

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

.. _hafs-documented:

========================================================
Hurricane Analysis and Reforecast System Configurations
========================================================

The HAFS configuration uses an :term:`DATM`-only configuration.

These tests use the default values set in the ``export_fv3``, ``export_hafs``, ``export_hafs_regional``, ``export_hafs_datm_cdeps``, and ``export_hafs_docn_cdeps`` functions of ``default_vars.sh`` unless other values are explicitly set in a given test file. In all tests, the values in ``export_fv3`` are set first. 

.. note:: 

   ``export_hafs`` calls ``export_hafs_regional``, which calls ``export_hafs_datm_cdeps`` or ``export_hafs_docn_cdeps``, which calls ``export_fv3``. Values from ``export_fv3`` are set first, followed by values in ``export_hafs``, ``export_hafs_regional``, and then values in ``export_hafs_datm_cdeps`` or ``export_hafs_docn_cdeps``. 


.. list-table:: *Default physics-related variables used in the HAFS configurations below*
   :widths: 10 50
   :header-rows: 1

   * - Export Function
     - Variables
   * - export_hafs
     - **Set to FALSE:** S2S, AQM, DATM_CDEPS, DOCN_CDEPS, HYBEDMF, CNVGWD, LTAEROSOL, LHEATSTRG, IS_MOVING_NEST :raw-html:`<br/> <br/>`
       **Set to TRUE:** FV3, HAFS, SATMEDMF, HURR_PBL, DO_GSL_DRAG_LS_BL, DO_GSL_DRAG_SS, DO_GSL_DRAG_TOFD, LRADAR, CPL_IMP_MRG :raw-html:`<br/> <br/>`
       **Set to VALUE:** NTILES=1, IMFSHALCNV=2, IMFDEEPCNV=2, MONINQ_FAC=-1.0, ISATMEDMF=1, IOPT_SFC=1, IOPT_DVEG=2, IOPT_CRS=1, IOPT_RAD=1, IOPT_ALB=2, IOPT_STC=1, LSM=1, IMP_PHYSICS=11, IAER=111, CDMBWD=1.0,1.0,1.0,1.0, FV_CORE_TAU=5., RF_CUTOFF=30.e2, RF_CUTOFF_NEST=50.e2, VORTEX_TRACKER=0, NTRACK=0, MOVE_CD_X=0, MOVE_CD_Y=0, NFHOUT=3, NFHMAX_HF=-1, NFHOUT_HF=3, NSOUT=-1, OUTPUT_FH=-1
   * - export_hafs_regional
     - **Set to FALSE:** S2S, AQM, DOCN_CDEPS, WRITE_DOPOST, USE_COLDSTART, MULTIGRID :raw-html:`<br/> <br/>`
       **Set to TRUE:** FV3, HAFS, CPL, QUILTING, OUTPUT_HISTORY, CPL_IMP_MRG :raw-html:`<br/> <br/>`
       **Set to VALUE:** NTILES=1, FHMAX=6, ENS_NUM=1, DT_ATMOS=900, RESTART_INTERVAL=0, FHROT=0, coupling_interval_fast_sec=0, WRITE_GROUP=1, WRTTASK_PER_GROUP=6, NUM_FILES=2, FILENAME_BASE="'atm' 'sfc'", OUTPUT_GRID="'regional_latlon'", OUTPUT_FILE="'netcdf'", IDEFLATE=0, QUANTIZE_NSD=0, NFHOUT=3, NFHMAX_HF=-1, NFHOUT_HF=3, CEN_LON=-62.0, CEN_LAT=25.0, LON1=-114.5, LAT1=-5.0, LON2=-9.5, LAT2=55.0, DLON=0.03, DLAT=0.03, DIAG_TABLE=diag_table_hafs, FIELD_TABLE=field_table_hafs, WW3OUTDTHR=3, OUTPARS_WAV="WND HS T01 T02 DIR FP DP PHS PTP PDIR UST CHA USP", WAV_CUR='C', med_model=cmeps, pio_rearranger=box, CAP_DBUG_FLAG=0, CPLMODE=hafs, RUNTYPE=startup, MESH_WAV=mesh.hafs.nc, MODDEF_WAV=mod_def.natl_6m 
   * - export_hafs_datm_cdeps
     - **Set to FALSE:** FV3, S2S, AQM, DOCN_CDEPS :raw-html:`<br/> <br/>`
       **Set to TRUE:** HAFS, DATM_CDEPS :raw-html:`<br/> <br/>`
       **Set to VALUE:** NTILES=1, atm_model=datm, DATM_IN_CONFIGURE=datm_in, DATM_STREAM_CONFIGURE=hafs_datm.streams.era5.IN
   * - export_hafs_docn_cdeps
     - **Set to FALSE:** S2S, AQM :raw-html:`<br/> <br/>`
       **Set to TRUE:** FV3, HAFS, DOCN_CDEPS  :raw-html:`<br/> <br/>`
       **Set to VALUE:** NTILES=1, ocn_model=docn, ocn_datamode=sstdata, pio_rearranger=box, DOCN_IN_CONFIGURE=docn_in, DOCN_STREAM_CONFIGURE=hafs_docn.streams.IN

.. _hafs-rts:

.. list-table:: *HAFS regression test descriptions*
   :widths: 50 10 30 50 10 10 10 10 10
   :header-rows: 1

   * - Test |nbsp| Name |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp| |nbsp|
     - Description
     - General Physics Parameters
     - Detailed |nbsp| Physics |nbsp| Parameters |nbsp| (see |nbsp| namelist |nbsp| options `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_c_c_p_psuite_nml_desp.html>`__ |nbsp| for variable definitions)
     - Start |nbsp| Date |nbsp| |nbsp| |nbsp| |nbsp|
     - Fcst Length (hours)
     - Output Grid
     - Configuration Files
     - Other
   * - `rhafs_global_1nest_atm <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_global_1nest_atm>`__
     - Compare HAFS global with 1 nest and atmosphere only results with previous trunk version
     - **Suite:** CCPP_SUITE="FV3_HAFS_v1_gfdlmp_tedmf"

       **Microphysics:** IMP_PHYSICS=11

       **Time Step:** DT_ATMOS=90
     - **Set to FALSE:** MOUNTAIN, WARM_START, FULL_ZS_FILTER, CPLFLX, CPLWAV, CPLWAV2ATM, CPL_IMP_MRG, CMEPS, USE_COLDSTART :raw-html:`<br/> <br/>`
       **Set to TRUE:** EXTERNAL_IC, NGGPS_IC, CPLOCN2ATM, NESTED :raw-html:`<br/> <br/>`
       **Set to VALUE:** 
       See ``export_hafs`` default values.
     - 2020-08-25 12:00:00
     - 6
     - OUTPUT_GRID=global_latlon, OUTPUT_GRID_2=rotated_latlon :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=4, JNPES=5, NPX=97, NPY=97, NPZ=64, NPZP=$(($NPZ + 1)), INPES_NEST02=6, JNPES_NEST02=10, NPX_NEST02=241, NPY_NEST02=241
     - FIELD_TABLE=field_table_hafs
       DIAG_TABLE=diag_table_hafs_template
       INPUT_NML=input_global_hafs.nml.IN
       INPUT_NEST02_NML=input_nest_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm.IN"
       FV3_RUN="hafs_fv3_run.IN"
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.false., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"
   * - `hafs_global_multiple_4nests_atm <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_global_multiple_4nests_atm>`__
     - Compare HAFS global with 4 multiple nests and atmosphere only results with previous trunk version
     - **Suite:** CCPP_SUITE="FV3_HAFS_v1_gfdlmp_tedmf"

       **Microphysics:** IMP_PHYSICS=11

       **Time Step:** DT_ATMOS=90
     - **Set to FALSE:** MOUNTAIN, WARM_START, FULL_ZS_FILTER, CPLFLX, CPLWAV, CPLWAV2ATM, CPL_IMP_MRG, CMEPS, USE_COLDSTART :raw-html:`<br/> <br/>`
       **Set to TRUE:** WRITE_DOPOST, EXTERNAL_IC, NGGPS_IC, CPLOCN2ATM, NESTED :raw-html:`<br/> <br/>`
       **Set to VALUE:**
       Also, see export_hafs default values.
     - 2020-08-25 12:00:00
     - 6
     - OUTPUT_GRID=global_latlon, OUTPUT_GRID_2=regional_latlon, OUTPUT_GRID_3=rotated_latlon, OUTPUT_GRID_4=rotated_latlon, OUTPUT_GRID_5=rotated_latlon :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=4, JNPES=5, NPX=97, NPY=97, NPZ=64, NPZP=$(($NPZ + 1)), INPES_NEST02=6, JNPES_NEST02=10, NPX_NEST02=241, NPY_NEST02=241, INPES_NEST03=6, JNPES_NEST03=10, NPX_NEST03=241, NPY_NEST03=241, INPES_NEST04=6, JNPES_NEST04=10, NPX_NEST04=361, NPY_NEST04=361, INPES_NEST05=6, JNPES_NEST05=10, NPX_NEST05=361, NPY_NEST05=361 
     - FIELD_TABLE=field_table_hafs
       DIAG_TABLE=diag_table_hafs_template
       INPUT_NML=input_global_hafs.nml.IN
       INPUT_NEST02_NML=input_nest_hafs.nml.IN
       INPUT_NEST03_NML=input_nest_hafs.nml.IN
       INPUT_NEST04_NML=input_nest_hafs.nml.IN
       INPUT_NEST05_NML=input_nest_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm.IN"
       FV3_RUN="hafs_fv3_run.IN"
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.false., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"
   * - `hafs_global_storm_following_1nest_atm <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_global_storm_following_1nest_atm>`__
     - Compare HAFS global with 1 storm-following moving nest and atmosphere only results with previous trunk version
     - **Suite:** CCPP_SUITE="FV3_HAFS_v1_gfdlmp_tedmf"

       **Microphysics:** IMP_PHYSICS=11

       **Time Step:** DT_ATMOS=180
     - **Set to FALSE:** MOUNTAIN, WARM_START, FULL_ZS_FILTER, IS_MOVING_NEST=".false.,.true.", CPLFLX, CPLWAV, CPLWAV2ATM, CPL_IMP_MRG, CMEPS, USE_COLDSTART :raw-html:`<br/> <br/>`
       **Set to TRUE:** EXTERNAL_IC, NGGPS_IC, CPLOCN2ATM, NESTED :raw-html:`<br/> <br/>`
       **Set to VALUE:** 
       Also, see export_hafs default values.
     - 2020-08-25 12:00:00
     - 6
     - OUTPUT_GRID=global_latlon, OUTPUT_GRID_2=rotated_latlon :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=4, JNPES=5, NPX=97, NPY=97, NPZ=64, NPZP=$(($NPZ + 1)), INPES_NEST02=6, JNPES_NEST02=10, NPX_NEST02=73, NPY_NEST02=73
     - FIELD_TABLE=field_table_hafs
       DIAG_TABLE=diag_table_hafs_template
       INPUT_NML=input_global_hafs.nml.IN
       INPUT_NEST02_NML=input_nest_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm.IN"
       FV3_RUN="hafs_fv3_run.IN"
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.false., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"
   * - `hafs_regional_1nest_atm <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_regional_1nest_atm>`__
     - Compare HAFS regional with 1 nest and atmosphere only results with previous trunk version
     - **Suite:** CCPP_SUITE="FV3_HAFS_v1_gfdlmp_tedmf"

       **Microphysics:** IMP_PHYSICS=11

       **Time Step:** DT_ATMOS=90
     - **Set to FALSE:** MOUNTAIN, WARM_START, FULL_ZS_FILTER, CPLFLX, CPLWAV, CPLWAV2ATM, CPL_IMP_MRG, CMEPS, USE_COLDSTART :raw-html:`<br/> <br/>`
       **Set to TRUE:** EXTERNAL_IC, NGGPS_IC, REGIONAL, CPLOCN2ATM, NESTED :raw-html:`<br/> <br/>`
       **Set to VALUE:**
       Also, see export_hafs default values.
     - 2020-08-25 12:00:00
     - 6
     - OUTPUT_GRID=rotated_latlon, OUTPUT_GRID_2=rotated_latlon :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=6, JNPES=10, NPX=241, NPY=241, NPZ=64, NPZP=$(($NPZ + 1)), INPES_NEST02=6, JNPES_NEST02=10, NPX_NEST02=361, NPY_NEST02=361
     - FIELD_TABLE=field_table_hafs
       DIAG_TABLE=diag_table_hafs_template
       INPUT_NML=input_regional_hafs.nml.IN
       INPUT_NEST02_NML=input_nest_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm.IN"
       FV3_RUN="hafs_fv3_run.IN"
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.false., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"
   * - `hafs_regional_atm <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_regional_atm>`__
     - Compare HAFS regional atmosphere only results with previous trunk version
     - **Suite:** CCPP_SUITE="FV3_HAFS_v1_gfdlmp_tedmf"

       **Microphysics:** IMP_PHYSICS=11

       **Time Step:** DT_ATMOS=180
     - **Set to FALSE:** MOUNTAIN, WARM_START, FULL_ZS_FILTER, CPLFLX, CPLWAV, CPLWAV2ATM, CPL_IMP_MRG, CMEPS, USE_COLDSTART :raw-html:`<br/> <br/>`
       **Set to TRUE:** EXTERNAL_IC, NGGPS_IC, REGIONAL, CPLOCN2ATM :raw-html:`<br/> <br/>`
       **Set to VALUE:** 
       Also, see export_hafs default values.
     - 2019-08-29 00:00:00
     - 6
     - OUTPUT_GRID=regional_latlon :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=20, JNPES=12, NPX=721, NPY=601, NPZ=91, NPZP=$(($NPZ + 1))
     - FIELD_TABLE=field_table_hafs
       DIAG_TABLE=diag_table_hafs_template
       INPUT_NML=input_regional_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm.IN"
       FV3_RUN="hafs_fv3_run.IN"
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.false., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"
   * - `hafs_regional_atm_ocn <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_regional_atm_ocn>`__
     - Compare HAFS regional atmosphere-ocean coupled HYCOM results with previous trunk version
     - **Suite:** CCPP_SUITE="FV3_HAFS_v1_gfdlmp_tedmf_nonsst"

       **Microphysics:** IMP_PHYSICS=11

       **Time Step:** DT_ATMOS=180
     - **Set to FALSE:** MOUNTAIN, WARM_START, FULL_ZS_FILTER, CPLWAV, CPLWAV2ATM, CDEPS_DOCN :raw-html:`<br/> <br/>`
       **Set to TRUE:** EXTERNAL_IC, NGGPS_IC, REGIONAL, CPLFLX, CPLOCN2ATM, CPL_IMP_MRG :raw-html:`<br/> <br/>`
       **Set to VALUE:** 
       Also, see export_hafs_regional then export_hafs default values.
     - 2019-08-29 00:00:00
     - 6
     - OUTPUT_GRID=regional_latlon :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=20, JNPES=12, NPX=721, NPY=601, NPZ=91, NPZP=$(($NPZ + 1))
     - FIELD_TABLE=field_table_hafs
       DIAG_TABLE=diag_table_hafs_template
       INPUT_NML=input_regional_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm_ocn.IN"
       FV3_RUN="hafs_fv3_run.IN hycom_hat10_run.IN"
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.false., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"
   * - `hafs_regional_atm_ocn_wav <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_regional_atm_ocn_wav>`__
     - Compare HAFS regional atmosphere-ocean-wave coupled results with previous trunk version
     - **Suite:** CCPP_SUITE="FV3_HAFS_v1_gfdlmp_tedmf_nonsst"

       **Microphysics:** IMP_PHYSICS=11

       **Time Step:** DT_ATMOS=180
     - **Set to FALSE:** MOUNTAIN, WARM_START, FULL_ZS_FILTER, CPLWAV2ATM, CDEPS_DOCN :raw-html:`<br/> <br/>`
       **Set to TRUE:** EXTERNAL_IC, NGGPS_IC, REGIONAL, CPLFLX, CPLOCN2ATM, CPLWAV, CPL_IMP_MRG :raw-html:`<br/> <br/>`
       **Set to VALUE:** 
       Also, see export_hafs_regional then export_hafs default values.
     - 2019-08-29 00:00:00
     - 6
     - OUTPUT_GRID=regional_latlon :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=20, JNPES=12, NPX=721, NPY=601, NPZ=91, NPZP=$(($NPZ + 1))
     - FIELD_TABLE=field_table_hafs
       DIAG_TABLE=diag_table_hafs_template
       INPUT_NML=input_regional_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm_ocn_wav.IN"
       FV3_RUN="hafs_fv3_run.IN hycom_hat10_run.IN hafs_ww3_run.IN"
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.false., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"
   * - `hafs_regional_atm_thompson_gfdlsf <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_regional_atm_thompson_gfdlsf>`__
     - Compare the results from HAFS regional atmosphere only using the Thompson microphysics scheme and GFDL surface layer scheme with previous trunk version
     - **Suite:** CCPP_SUITE="FV3_HAFS_v1_thompson_tedmf_gfdlsf"

       **Microphysics:** IMP_PHYSICS=8

       **Time Step:** DT_ATMOS=180
     - **Set to FALSE:** MOUNTAIN, WARM_START, FULL_ZS_FILTER, DO_SAT_ADJ, CPLFLX, CPLWAV, CPLWAV2ATM, CPL_IMP_MRG, CMEPS, USE_COLDSTART :raw-html:`<br/> <br/>`
       **Set to TRUE:** EXTERNAL_IC, NGGPS_IC, REGIONAL, CPLOCN2ATM :raw-html:`<br/> <br/>`
       **Set to VALUE:** 
       Also, see export_hafs default values.
     - 2019-08-29 00:00:00
     - 6
     - OUTPUT_GRID=cubed_sphere_grid :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=20, JNPES=12, NPX=721, NPY=601, NPZ=91, NPZP=$(($NPZ + 1))
     - FIELD_TABLE=field_table_hafs_thompson
       DIAG_TABLE=diag_table_hafs_template
       INPUT_NML=input_regional_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm.IN"
       FV3_RUN="hafs_fv3_run.IN"
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.false., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"
   * - `hafs_regional_atm_wav <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_regional_atm_wav>`__
     - Compare HAFS regional atmosphere-wave coupled results with previous trunk version
     - **Suite:** CCPP_SUITE="FV3_HAFS_v1_gfdlmp_tedmf"

       **Microphysics:** IMP_PHYSICS=11

       **Time Step:** DT_ATMOS=180
     - **Set to FALSE:** MOUNTAIN, WARM_START, FULL_ZS_FILTER, CPLOCN2ATM, CDEPS_DOCN :raw-html:`<br/> <br/>`
       **Set to TRUE:** EXTERNAL_IC, NGGPS_IC, REGIONAL, CPLFLX, CPLWAV, CPLWAV2ATM, CPL_IMP_MRG :raw-html:`<br/> <br/>`
       **Set to VALUE:** 
       Also, see export_hafs_regional then export_hafs default values.
     - 2019-08-29 00:00:00
     - 6
     - OUTPUT_GRID=regional_latlon :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=20, JNPES=12, NPX=721, NPY=601, NPZ=91, NPZP=$(($NPZ + 1))
     - FIELD_TABLE=field_table_hafs
       DIAG_TABLE=diag_table_hafs_template
       INPUT_NML=input_regional_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm_wav.IN"
       FV3_RUN="hafs_fv3_run.IN hafs_ww3_run.IN"
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.false., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"
   * - `hafs_regional_datm_cdeps <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_regional_datm_cdeps>`__
     - Compare HAFS regional coupled CDEPS data atmosphere from ERA5 with regional HYCOM results with previous trunk version
     - N/A: No active atmospheric component
     - **Set to FALSE:** CPLWAV, CDEPS_DOCN :raw-html:`<br/> <br/>`
       **Set to TRUE:** :raw-html:`<br/> <br/>`
       **Set to VALUE:** 
       Also, see export_hafs_datm_cdeps then export_hafs_regional then export_hafs default values.
     - 2019-08-29 00:00:00
     - 24
     - OUTPUT_GRID=regional_latlon :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=$INPES_dflt, JNPES=$JNPES_dflt
     - FIELD_TABLE=field_table_hafs
       DIAG_TABLE=diag_table_hafs
       INPUT_NML=input_regional_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm_ocn.IN"
       FV3_RUN="hafs_datm_cdeps_era5.IN hycom_hat10_run.IN"
       DATM_STREAM_CONFIGURE=hafs_datm.streams.era5.IN
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.false., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"
   * - `hafs_regional_docn <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_regional_docn>`__
     - Compare HAFS regional coupled with regional data ocean from MOM6 results with previous trunk version
     - **Suite:** CCPP_SUITE="FV3_HAFS_v1_gfdlmp_tedmf_nonsst"

       **Microphysics:** IMP_PHYSICS=11

       **Time Step:** DT_ATMOS=180
     - **Set to FALSE:** MOUNTAIN, WARM_START, FULL_ZS_FILTER, CPLWAV, CPLWAV2ATM :raw-html:`<br/> <br/>`
       **Set to TRUE:** EXTERNAL_IC, NGGPS_IC, REGIONAL, CPLFLX, CPLOCN2ATM, CPL_IMP_MRG :raw-html:`<br/> <br/>`
       **Set to VALUE:** 
       Also, see export_hafs_docn_cdeps then export_hafs_regional then export_hafs default values.
     - 2019-08-29 00:00:00
     - 24
     - OUTPUT_GRID=regional_latlon :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=20, JNPES=12, NPX=721, NPY=601, NPZ=91, NPZP=$(($NPZ + 1))
     - FIELD_TABLE=field_table_hafs
       DIAG_TABLE=diag_table_hafs_template
       INPUT_NML=input_regional_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm_docn.IN"
       FV3_RUN="hafs_fv3_run.IN hafs_docn_cdeps_mom6.IN"
       DOCN_STREAM_CONFIGURE=hafs_docn.streams.IN
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.false., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"
   * - `hafs_regional_docn_oisst <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_regional_docn_oisst>`__
     - Compare HAFS regional coupled with global data ocean from OISST results with previous trunk version
     - **Suite:** CCPP_SUITE=FV3_HAFS_v1_gfdlmp_tedmf_nonsst

       **Microphysics:** IMP_PHYSICS=11

       **Time Step:** DT_ATMOS=180
     - **Set to FALSE:** MOUNTAIN, WARM_START, FULL_ZS_FILTER, CPLWAV, CPLWAV2ATM :raw-html:`<br/> <br/>`
       **Set to TRUE:** EXTERNAL_IC, NGGPS_IC, REGIONAL, CPLFLX, CPLOCN2ATM, CPL_IMP_MRG :raw-html:`<br/> <br/>`
       **Set to VALUE:** 
       Also, see export_hafs_docn_cdeps then export_hafs_regional then export_hafs default values.
     - 2019-08-29 00:00:00
     - 6
     - OUTPUT_GRID=regional_latlon :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=20, JNPES=12, NPX=721, NPY=601, NPZ=91, NPZP=$(($NPZ + 1))
     - FIELD_TABLE=field_table_hafs
       DIAG_TABLE=diag_table_hafs_template
       INPUT_NML=input_regional_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm_docn.IN"
       FV3_RUN="hafs_fv3_run.IN hafs_docn_cdeps_oisst.IN"
       DOCN_STREAM_CONFIGURE=hafs_docn.streams.IN
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.true., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"
   * - `hafs_regional_specified_moving_1nest_atm <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_regional_specified_moving_1nest_atm>`__
     - Compare HAFS regional with 1 specified moving nest and atmosphere only results with previous trunk version
     - **Suite:** CCPP_SUITE="FV3_HAFS_v1_gfdlmp_tedmf"

       **Microphysics:** IMP_PHYSICS=11

       **Time Step:** DT_ATMOS=180
     - **Set to FALSE:** MOUNTAIN, WARM_START, FULL_ZS_FILTER, IS_MOVING_NEST=".false.,.true.", CPLFLX, CPLWAV, CPLWAV2ATM, CPL_IMP_MRG, CMEPS, USE_COLDSTART :raw-html:`<br/> <br/>`
       **Set to TRUE:** WRITE_DOPOST, EXTERNAL_IC, NGGPS_IC, REGIONAL, CPLOCN2ATM :raw-html:`<br/> <br/>`
       **Set to VALUE:** 
       Also, see export_hafs default values.
     - 2020-08-25 12:00:00
     - 6
     - OUTPUT_GRID=rotated_latlon, OUTPUT_GRID_2=rotated_latlon_moving :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=6, JNPES=10, NPX=241, NPY=241, NPZ=64, NPZP=$(($NPZ + 1)), INPES_NEST02=6, JNPES_NEST02=10, NPX_NEST02=361, NPY_NEST02=361
     - FIELD_TABLE=field_table_hafs
       DIAG_TABLE=diag_table_hafs_template
       INPUT_NML=input_regional_hafs.nml.IN
       INPUT_NEST02_NML=input_nest_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm.IN"
       FV3_RUN="hafs_fv3_run.IN"
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.false., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"
   * - `hafs_regional_storm_following_1nest_atm <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_regional_storm_following_1nest_atm>`__
     - Compare HAFS regional with 1 storm-following moving nest and atmosphere only results with previous trunk version
     - **Suite:** CCPP_SUITE="FV3_HAFS_v1_gfdlmp_tedmf"

       **Microphysics:** IMP_PHYSICS=11

       **Time Step:** DT_ATMOS=180
     - **Set to FALSE:** MOUNTAIN, WARM_START, FULL_ZS_FILTER, IS_MOVING_NEST=".false.,.true.", CPLFLX, CPLWAV, CPLWAV2ATM, CPL_IMP_MRG, CMEPS, USE_COLDSTART :raw-html:`<br/> <br/>`
       **Set to TRUE:** EXTERNAL_IC, NGGPS_IC, REGIONAL, CPLOCN2ATM :raw-html:`<br/> <br/>`
       **Set to VALUE:** 
       Also, see export_hafs default values.
     - 2020-08-25 12:00:00
     - 6
     - OUTPUT_GRID=rotated_latlon, OUTPUT_GRID_2=rotated_latlon_moving :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=6, JNPES=10, NPX=241, NPY=241, NPZ=64, NPZP=$(($NPZ + 1)), INPES_NEST02=6, JNPES_NEST02=10, NPX_NEST02=361, NPY_NEST02=361
     - FIELD_TABLE=field_table_hafs
       DIAG_TABLE=diag_table_hafs_template
       INPUT_NML=input_regional_hafs.nml.IN
       INPUT_NEST02_NML=input_nest_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm.IN"
       FV3_RUN="hafs_fv3_run.IN"
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.false., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"
   * - `hafs_regional_storm_following_1nest_atm_ocn <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_regional_storm_following_1nest_atm_ocn>`__
     - Compare HAFS regional with 1 storm-following moving nest and atmosphere-ocean coupled results with previous trunk version
     - **Suite:** CCPP_SUITE="FV3_HAFS_v1_gfdlmp_tedmf_nonsst"

       **Microphysics:** IMP_PHYSICS=11

       **Time Step:** DT_ATMOS=180
     - **Set to FALSE:** MOUNTAIN, WARM_START, FULL_ZS_FILTER, IS_MOVING_NEST=".false.,.true.", CPLWAV, CPLWAV2ATM, USE_COLDSTART, CDEPS_DOCN :raw-html:`<br/> <br/>`
       **Set to TRUE:** EXTERNAL_IC, NGGPS_IC, REGIONAL, CPLFLX, CPLOCN2ATM, CPL_IMP_MRG :raw-html:`<br/> <br/>`
       **Set to VALUE:** 
       Also, see export_hafs_regional default values then export_hafs.
     - 2020-08-25 12:00:00
     - 6
     - OUTPUT_GRID=regional_latlon, OUTPUT_GRID_2=regional_latlon_moving :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=6, JNPES=10, NPX=241, NPY=241, NPZ=64, NPZP=$(($NPZ + 1)), INPES_NEST02=6, JNPES_NEST02=10, NPX_NEST02=361, NPY_NEST02=361
     - FIELD_TABLE=field_table_hafs
       DIAG_TABLE=diag_table_hafs_template
       INPUT_NML=input_regional_hafs.nml.IN
       INPUT_NEST02_NML=input_nest_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm_ocn.IN"
       FV3_RUN="hafs_fv3_run.IN hycom_hat10_run.IN"
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.false., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"
   * - `hafs_regional_storm_following_1nest_atm_ocn_debug <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_regional_storm_following_1nest_atm_ocn_debug>`__
     - Compare HAFS regional with 1 storm-following moving nest and atmosphere-ocean coupled results with previous trunk version
     - **Suite:** CCPP_SUITE="FV3_HAFS_v1_gfdlmp_tedmf_nonsst"

       **Microphysics:** IMP_PHYSICS=11

       **Time Step:** DT_ATMOS=180
     - **Set to FALSE:** MOUNTAIN, WARM_START, FULL_ZS_FILTER, IS_MOVING_NEST=".false.,.true.", CPLWAV, CPLWAV2ATM, USE_COLDSTART, CDEPS_DOCN :raw-html:`<br/> <br/>`
       **Set to TRUE:** EXTERNAL_IC, NGGPS_IC, REGIONAL, CPLFLX, CPLOCN2ATM, CPL_IMP_MRG :raw-html:`<br/> <br/>`
       **Set to VALUE:** 
       Also, see export_hafs_regional default values then export_hafs.
     - 2020-08-25 12:00:00
     - 6
     - OUTPUT_GRID=regional_latlon, OUTPUT_GRID_2=regional_latlon_moving :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=6, JNPES=10, NPX=241, NPY=241, NPZ=64, NPZP=$(($NPZ + 1)), INPES_NEST02=6, JNPES_NEST02=10, NPX_NEST02=361, NPY_NEST02=361
     - FIELD_TABLE=field_table_hafs
       DIAG_TABLE=diag_table_hafs_template
       INPUT_NML=input_regional_hafs.nml.IN
       INPUT_NEST02_NML=input_nest_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm_ocn.IN"
       FV3_RUN="hafs_fv3_run.IN hycom_hat10_run.IN"
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.false., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"
   * - `hafs_regional_storm_following_1nest_atm_ocn_debug <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_regional_storm_following_1nest_atm_ocn_debug>`__
     - Compare HAFS regional with 1 storm-following moving nest and atmosphere-ocean coupled results with previous trunk version
     - **Suite:** CCPP_SUITE="FV3_HAFS_v1_gfdlmp_tedmf_nonsst"

       **Microphysics:** IMP_PHYSICS=11

       **Time Step:** DT_ATMOS=180
     - **Set to FALSE:** MOUNTAIN, WARM_START, FULL_ZS_FILTER, IS_MOVING_NEST=".false.,.true.", CPLWAV2ATM, USE_COLDSTART, CDEPS_DOCN :raw-html:`<br/> <br/>`
       **Set to TRUE:** EXTERNAL_IC, NGGPS_IC, REGIONAL, CPLFLX, CPLOCN2ATM, CPL_IMP_MRG :raw-html:`<br/> <br/>`
       **Set to VALUE:** 
       Also, see export_hafs_regional default values then export_hafs.
     - 2020-08-25 12:00:00
     - 6
     - OUTPUT_GRID=regional_latlon, OUTPUT_GRID_2=regional_latlon_moving :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=$INPES_thrd, JNPES=$JNPES_thrd, INPES=6, JNPES=10, NPX=241, NPY=241, NPZ=64, NPZP=$(($NPZ + 1))
     - FIELD_TABLE=field_table_hafs
       DIAG_TABLE=diag_table_hafs_template
       INPUT_NML=input_regional_hafs.nml.IN
       INPUT_NEST02_NML=input_nest_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm_ocn.IN"
       FV3_RUN="hafs_fv3_run.IN hycom_hat10_run.IN"
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.false., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"
   * - `hafs_regional_storm_following_1nest_atm_ocn_wav <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_regional_storm_following_1nest_atm_ocn_wav>`__
     - Compare HAFS regional with 1 storm-following moving nest and atmosphere-ocean-wave coupled results with previous trunk version
     - **Suite:** CCPP_SUITE="FV3_HAFS_v1_gfdlmp_tedmf_nonsst"

       **Microphysics:** IMP_PHYSICS=11

       **Time Step:** DT_ATMOS=180
     - **Set to FALSE:** MOUNTAIN, WARM_START, FULL_ZS_FILTER, IS_MOVING_NEST=".false.,.true.", CPLWAV2ATM, USE_COLDSTART, CDEPS_DOCN :raw-html:`<br/> <br/>`
       **Set to TRUE:** EXTERNAL_IC, NGGPS_IC, REGIONAL, CPLFLX, CPLOCN2ATM, CPLWAV, CPL_IMP_MRG :raw-html:`<br/> <br/>`
       **Set to VALUE:** 
       Also, see export_hafs_regional default values then export_hafs.
     - 2020-08-25 12:00:00
     - 6
     - OUTPUT_GRID=rotated_latlon, OUTPUT_GRID_2=rotated_latlon :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=6, JNPES=10, NPX=241, NPY=241, NPZ=64, NPZP=$(($NPZ + 1)), INPES_NEST02=6, JNPES_NEST02=10, NPX_NEST02=361, NPY_NEST02=361
     - FIELD_TABLE=field_table_hafs
       DIAG_TABLE=diag_table_hafs_template
       INPUT_NML=input_regional_hafs.nml.IN
       INPUT_NEST02_NML=input_nest_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm_ocn_wav.IN"
       FV3_RUN="hafs_fv3_run.IN hycom_hat10_run.IN hafs_ww3_run.IN"
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.false., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"
   * - `hafs_regional_telescopic_2nests_atm <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/hafs_regional_telescopic_2nests_atm>`__
     - Compare HAFS regional with two telescopic nests and atmosphere only results with previous trunk version
     - **Suite:** CCPP_SUITE="FV3_HAFS_v1_gfdlmp_tedmf"

       **Microphysics:** IMP_PHYSICS=11

       **Time Step:** DT_ATMOS=90
     - **Set to FALSE:** MOUNTAIN, WARM_START, FULL_ZS_FILTER, CPLFLX, CPLWAV, CPLWAV2ATM, CMEPS, USE_COLDSTART :raw-html:`<br/> <br/>`
       **Set to TRUE:** EXTERNAL_IC, NGGPS_IC, REGIONAL, CPLOCN2ATM :raw-html:`<br/> <br/>`
       **Set to VALUE:** 
       Also, see export_hafs default values.
     - 2020-08-25 12:00:00
     - 6
     - OUTPUT_GRID=rotated_latlon, OUTPUT_GRID_2=lambert_conformal, OUTPUT_GRID_3=regional_latlon :raw-html:`<br/> <br/>`
       **Grid Parameters**: INPES=6, JNPES=10, NPX=241, NPY=241, NPZ=64, NPZP=$(($NPZ + 1)), INPES_NEST02=6, JNPES_NEST02=10, NPX_NEST02=361, NPY_NEST02=361, INPES_NEST03=6, JNPES_NEST03=10, NPX_NEST03=361, NPY_NEST03=361
     - FIELD_TABLE=field_table_hafs
       DIAG_TABLE=diag_table_hafs_template
       INPUT_NML=input_regional_hafs.nml.IN
       INPUT_NEST02_NML=input_nest_hafs.nml.IN
       INPUT_NEST03_NML=input_nest_hafs.nml.IN
       MODEL_CONFIGURE="model_configure_hafs.IN"
       UFS_CONFIGURE="ufs.configure.hafs_atm.IN"
       FV3_RUN="hafs_fv3_run.IN"
     - RESTART_INTERVAL=1, atm_omp_num_threads=2, WARM_START=.false., READ_INCREMENT=.false., RES_LATLON_DYNAMICS="'fv3_increment.nc'"

**Sample** ``CMAKE_FLAGS`` **Setting**

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=HAFS"

**Supported Physics Suites**

.. list-table:: *Physics suites used in the HAFS configurations above*
   :widths: 10 50
   :header-rows: 1

   * - Physics Suite
     - Description
   * - FV3_HAFS_v1_gfdlmp_tedmf
     - The FV3_HAFS_v1_gfdlmp_tedmf physics suite is described in the :term:`CCPP` documentation `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/rap_suite_page.html>`__.
   * - FV3_HAFS_v1_gfdlmp_tedmf_nonsst
     - The FV3_HAFS_v1_gfdlmp_tedmf_nonsst physics suite is described in the CCPP documentation `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_h_r_r_r_suite_page.html>`__.
   * - FV3_HAFS_v1_thompson_tedmf_gfdlsf
     - The FV3_HAFS_v1_thompson_tedmf_gfdlsf physics suite is described in the CCPP documentation `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_r_r_f_s_v1beta_page.html>`__.

