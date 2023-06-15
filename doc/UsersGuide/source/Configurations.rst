.. _Configurations:

*************************
Configurations
*************************

The UFS Weather Model (WM) can be run in any of several configurations, from a single-component atmospheric 
model to a fully coupled model with multiple earth system components (e.g., atmosphere, ocean, sea-ice, land, and 
mediator). This chapter documents a few of the currently supported configurations. For a full list of 
supported configurations, view the `rt.conf <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/rt.conf>`__ 
and `rt.gnu.conf <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/rt_gnu.conf>`__ files. 

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
Additionally, there is a list of preferred physics suites, examples of ``nems.configure`` files, 
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
     - Physics Suite (see namelist options `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_c_c_p_psuite_nml_desp.html>`__)
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

Input files required for ATM configurations can be viewed in :numref:`Section %s <atm-io>`
or in the `UFS WM RT Data Bucket <https://registry.opendata.aws/noaa-ufs-regtests/>`__. 
Information on ``nems.configure`` files is available in :numref:`Section %s <nems-conf>`,
and a sample ATM ``nems.configure`` file (``nems.configure.atm.IN``) is available 
`here <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/parm/nems.configure.atm.IN>`__.


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
     - Physics Suite (see namelist options `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_c_c_p_psuite_nml_desp.html>`__)
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

Input files required for ATML configurations can be viewed in :numref:`Section %s (ATM) <atm-io>` 
and :numref:`Section %s (LND) <lnd-io>` or in the `UFS WM RT Data Bucket <https://registry.opendata.aws/noaa-ufs-regtests/>`__. 
Information on ``nems.configure`` files is available in :numref:`Section %s <nems-conf>`,
and a sample ATML ``nems.configure`` file (``nems.configure.atm_lnd.IN``) is available 
`here <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/parm/nems.configure.atm_lnd.IN>`__.


.. _rrfs-documented:

=======================================
Rapid Refresh Forecast System (RRFS)
=======================================

The RRFS configurations use an :term:`ATM`-only configuration on a high-resolution 
regional grid with data assimilation capabilities. 
These tests use default values set in the ``export_fv3`` function of ``default_vars.sh``. 

Current RRFS regression tests cover a wide variety of functionality and involve several 
physics tests. :numref:`Table %s <rrfs-rts>` contains RTs for RRFS functionality; 
in other words, : 


.. _rrfs-rts:

.. list-table:: *RRFS regression test descriptions*
   :widths: 10 10 10 10 10 10 10 10 10
   :header-rows: 1

   * - Test Name
     - Description
     - Physics Suite (see namelist options `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_c_c_p_psuite_nml_desp.html>`__)
     - DT_ATMOS
     - Start Date
     - Forecast Length (hours)
     - Output Grid
     - Configuration Files
     - Notes
   * - `rrfs_v1beta <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/rrfs_v1beta>`__ / `rrfs_v1beta_debug <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/rrfs_v1beta_debug>`__
     - Compare RRFS_v1beta / rrfs_v1beta_debug results with previous trunk version
     - FV3_RRFS_v1beta; IMP_PHYSICS=8
     - 300
     - 2021-03-22 06:00:00
     - 24
     - gaussian_grid
     - INPUT_NML=rap.nml.IN | FIELD_TABLE=field_table_thompson_aero_tke | DIAG_TABLE=diag_table_rap_noah
     - Notes: DO_SAT_ADJ=.false. ; LRADAR=.true. ; LTAEROSOL=.true. ; IALB=2; IEMS=2; HYBEDMF=.false. ; DO_MYNNEDMF=.true. ; DO_MYNNSFCLAY=.true. ; DO_DEEP=.false. ; SHAL_CNV=.false. ; IMFSHALCNV=-1 ; IMFDEEPCNV=-1 ; LHEATSTRG=.false. ; LSM=2 ; LSOIL_LSM=4
   * - `rrfs_v1nssl <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/rrfs_v1nssl>`__
     - Compare RRFS_v1nssl results with previous trunk version
     - FV3_RRFS_v1nssl; IMP_PHYSICS=17
     - 300
     - 2021-03-22 06:00:00
     - 24
     - gaussian_grid
     - INPUT_NML=rap.nml.IN | FIELD_TABLE=field_table_nssl_tke | DIAG_TABLE=diag_table_rap_noah; NSSL_CCN_ON=.true.
     - Notes: NSSL_HAIL_ON=.true. ; NSSL_INVERTCCN=.true. ; DO_SAT_ADJ=.false. ; LTAEROSOL=.false. ; IALB=2 ; IEMS=2 ; HYBEDMF=.false. ; DO_MYNNEDMF=.true. ; DO_MYNNSFCLAY=.true. ; DO_DEEP=.false. ; SHAL_CNV=.false. ; IMFSHALCNV=-1 ; IMFDEEPCNV=-1 ; LHEATSTRG=.false. ; LSM=2 ; LSOIL_LSM=4
   * - `rrfs_v1nssl_nohailnoccn <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/rrfs_v1nssl_nohailnoccn>`__
     - Compare RRFS_v1nssl_nohailnoccn results with previous trunk version
     - FV3_RRFS_v1nssl; IMP_PHYSICS=17
     - 300
     - 2021-03-22 06:00:00
     - 24
     - gaussian_grid
     - FV3_RUN=control_run.IN | INPUT_NML=rap.nml.IN | FIELD_TABLE=field_table_nssl_nohailnoccn_tke | DIAG_TABLE=diag_table_rap_noah
     - Notes: DNATS=0 ; NWAT=6 ; NSSL_CCN_ON=.false. ; NSSL_HAIL_ON=.false. ; NSSL_INVERTCCN=.true. ; DO_SAT_ADJ=.false. ; LTAEROSOL=.false. ; IALB=2 ; IEMS=2 ; HYBEDMF=.false. ; DO_MYNNEDMF=.true. ; DO_MYNNSFCLAY=.true. ; DO_DEEP=.false. ; SHAL_CNV=.false. ; IMFSHALCNV=-1 ; IMFDEEPCNV=-1 ; LHEATSTRG=.false. ; LSM=2 ; LSOIL_LSM=4
   * - `rrfs_conus13km_hrrr_warm <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/rrfs_conus13km_hrrr_warm>`__
     - HRRR physics on 13km domain, control run
     - FV3_HRRR; IMP_PHYSICS=8
     - 120
     - 2021-05-12 16:00:00
     - 2
     - lambert_conformal
     - FV3_RUN=rrfs_warm_run.IN | INPUT_NML=rrfs_conus13km_hrrr.nml.IN | FIELD_TABLE=field_table_thompson_aero_tke | DIAG_TABLE=diag_table_hrrr | MODEL_CONFIGURE=model_configure_rrfs_conus13km.IN
     - Notes: LKM=1; IOPT_LAKE=2; SFCLAY_COMPUTE_FLUX=.true.; IALB=2; ICLIQ_SW=2; IEMS=2; IOVR=3; KICE=9; LSM=3; LSOIL_LSM=9; DO_MYNNSFCLAY=.true.; DO_MYNNEDMF=.true.; DO_MYJPBL=.true; HYBEDMF=.false.; SHAL_CNV=.false.; DO_SAT_ADJ=.false.; DO_DEEP=.false. ; MAKE_NH=.false.; NA_INIT=0; DNATS=0; EXTERNAL_IC=.false.; NGGPS_IC=.false.; MOUNTAIN=.true.; WARM_START=.true.; READ_INCREMENT=.false.; RES_LATLON_DYNAMICS="'fv3_increment.nc'"; FHZERO=1.0; LDIAG3D=.false.; QDIAG3D=.false.; FHCYC=0.0; IAER=1011; LHEATSTRG=.false.; RANDOM_CLDS=.false.; CNVCLD=.false.; IMFSHALCNV=-1; IMFDEEPCNV=-1; CDMBWD='3.5,1.0'; DO_SPPT=.false.; DO_SHUM=.false.; DO_SKEB=.false.; LNDP_TYPE=0; N_VAR_LNDP=0; GWD_OPT=3; DO_UGWP_V0=.false.; DO_UGWP_V0_OROG_ONLY=.false.; DO_GSL_DRAG_LS_BL=.true.; DO_GSL_DRAG_SS=.true.; DO_GSL_DRAG_TOFD=.true.; DO_UGWP_V1=.false.; DO_UGWP_V1_OROG_ONLY=.false.; FRAC_ICE=.true.
   * - `rrfs_conus13km_hrrr_warm_debug <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/rrfs_conus13km_hrrr_warm_debug>`__
     - HRRR physics on 13km domain, debug run
     - FV3_HRRR; IMP_PHYSICS=8
     - 120
     - 2021-05-12 16:00:00
     - 1
     - lambert_conformal
     - 
     - Notes: 
   * - `rrfs_conus13km_hrrr_warm_restart_mismatch <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/rrfs_conus13km_hrrr_warm_restart_mismatch>`__
     - HRRR physics on 13km domain, control run
     - FV3_HRRR; IMP_PHYSICS=8
     - 120
     - 2021-05-12 16:00:00
     - 2
     - lambert_conformal
     - 
     - Notes: 
   * - `rrfs_smoke_conus13km_hrrr_warm <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/rrfs_smoke_conus13km_hrrr_warm>`__
     - HRRR smoke physics on 13km domain, control run
     - FV3_HRRR; IMP_PHYSICS=8
     - 120
     - 2021-05-12 16:00:00
     - 2
     - lambert_conformal
     - 
     - Notes: 
   * - `rrfs_smoke_conus13km_hrrr_warm_2threads <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/rrfs_smoke_conus13km_hrrr_warm_2threads>`__
     - HRRR smoke physics on 13km domain, different threads
     - FV3_HRRR; IMP_PHYSICS=8
     - 120
     - 2021-05-12 16:00:00
     - 2
     - lambert_conformal
     - 
     - Notes: 
   * - `rrfs_smoke_conus13km_hrrr_warm_debug <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/rrfs_smoke_conus13km_hrrr_warm_debug>`__
     - HRRR smoke physics on 13km domain, control run
     - FV3_HRRR; IMP_PHYSICS=8
     - 120
     - 2021-05-12 16:00:00
     - 1
     - lambert_conformal
     - 
     - Notes: 
   * - `rrfs_smoke_conus13km_hrrr_warm_debug_2threads <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/rrfs_smoke_conus13km_hrrr_warm_debug_2threads>`__
     - HRRR smoke physics on 13km domain, control run
     - FV3_HRRR; IMP_PHYSICS=8
     - 120
     - 2021-05-12 16:00:00
     - 1
     - lambert_conformal
     - 
     - Notes: 
   * - `rrfs_smoke_conus13km_radar_tten_warm <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/tests/rrfs_smoke_conus13km_radar_tten_warm>`__
     - HRRR smoke physics on 13km domain with radar-derived temperature tendencies
     - FV3_HRRR; IMP_PHYSICS=8
     - 120
     - 2021-05-12 16:00:00
     - 2
     - lambert_conformal
     - 
     - Notes: 

.. COMMENT: What are PEs??? And check rrfs_conus13km_hrrr_warm_restart_mismatch description. It's the same as the rrfs_conus13km_hrrr_warm description

**Sample** ``CMAKE_FLAGS`` **Setting**

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATM -DCCPP_SUITES=FV3_RAP,FV3_RAP_sfcdiff,FV3_HRRR,FV3_HRRR_flake,FV3_RRFS_v1beta,FV3_RRFS_v1nssl -D32BIT=ON"

.. COMMENT: Edit this section! 

**Supported Physics Suites**

.. COMMENT: Edit this section! 

.. list-table:: *Physics suites used in the RRFS configurations above*
   :widths: 10 50
   :header-rows: 1

   * - Physics Suite
     - Description
   * - FV3_RAP
     - The FV3_RAP physics suite is described in the :term:`CCPP` documentation `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/rap_suite_page.html>`__.
   * - FV3_RAP_sfcdiff
     - The FV3_RAP_sfcdiff physics suite is described in the CCPP documentation `here <>`__. 
   * - FV3_HRRR
     - The FV3_HRRR physics suite is described in the CCPP documentation `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_h_r_r_r_suite_page.html>`__.
   * - FV3_HRRR_flake
     - The FV3_HRRR_flake physics suite is described in the CCPP documentation `here <>`__.
   * - FV3_RRFS_v1beta 
     - The FV3_RRFS_v1beta physics suite is described in the CCPP documentation `here <https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_r_r_f_s_v1beta_page.html>`__.
   * - FV3_RRFS_v1nssl
     - The FV3_RRFS_v1nssl physics suite is described in the CCPP documentation `here <>`__.

**Additional Information**

.. COMMENT: Edit this section! 

Input files required for RRFS ATM configurations can be viewed in :numref:`Table %s <rrfs-files>`
or in the `UFS WM RT Data Bucket <https://registry.opendata.aws/noaa-ufs-regtests/>`__. 

.. COMMENT: Edit this section before adding!

   Information on ``nems.configure`` files is available in :numref:`Section %s <nems-conf>`,
   and a sample RRFS ``nems.configure`` file (``nems.configure.atm.IN``) is available 
   `here <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/parm/nems.configure.atm.IN>`__.

.. _rrfs-files:

.. list-table:: Files Required for RRFS RTs
   :widths: 50 10 10 10 10 100
   :header-rows: 1

   * - Tests
     - sfcf*.nc
     - atmf*.nc
     - GFSFLX.GrbF*
     - GFSPRS.GrbF*
     - Other
   * - rrfs_v1beta
     - sfcf000.nc
       sfcf009.nc
       sfcf012.nc
     - atmf000.nc
       atmf009.nc
       atmf012.nc
     - GFSFLX.GrbF00
       GFSFLX.GrbF09
       GFSFLX.GrbF12
     - GFSPRS.GrbF00
       GFSPRS.GrbF09
       GFSPRS.GrbF12
     - coupler.research
       
       fv_core.res.nc
       
       fv_core.res.tile[`1\nobreakdash-6`].nc &#8209
       
       fv_srf_wnd.res.tile[1\\nobreakdash-6].nc
       
       fv_tracer.res.tile[1\nobreakdash-6].nc
       
       phy_data.tile[:math:`1\\nobreakdash-6`].nc        
       
       sfc_data.tile[:math:`1\nobreakdash-6`].nc
   * - rrfs_v1beta_debug
       rrfs_conus13km_hrrr_warm_debug
       rrfs_smoke_conus13km_hrrr_warm_debug
       rrfs_smoke_conus13km_hrrr_warm_debug_2threads
     - sfcf000.nc
       sfcf001.nc
     - atmf000.nc
       atmf001.nc
     - 
     - 
     - 
   * - rrfs_v1nssl
       rrfs_v1nssl_nohailnoccn
     - sfcf000.nc
       sfcf009.nc
       sfcf012.nc
     - atmf000.nc
       atmf009.nc
       atmf012.nc
     - GFSFLX.GrbF00
       GFSFLX.GrbF09
       GFSFLX.GrbF12
     - GFSPRS.GrbF00
       GFSPRS.GrbF09
       GFSPRS.GrbF12
     - 
   * - rrfs_conus13km_hrrr_warm
       rrfs_smoke_conus13km_hrrr_warm
       rrfs_smoke_conus13km_hrrr_warm_2threads
       rrfs_conus13km_radar_tten_warm
     - sfcf000.nc
       sfcf001.nc
       sfcf002.nc
     - atmf000.nc
       atmf001.nc
       atmf002.nc
     - atmf000.nc
       atmf009.nc
       atmf012.nc
     - 
     - 
   * - rrfs_conus13km_hrrr_warm_restart_mismatch
     - sfcf002.nc
     - atmf002.nc
     - 
     - 
     - 

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

Input files required for LND configurations can be viewed in :numref:`Section %s (LND) <lnd-io>` 
or in the `UFS WM RT Data Bucket <https://registry.opendata.aws/noaa-ufs-regtests/>`__. 
Information on ``nems.configure`` files is available in :numref:`Section %s <nems-conf>`,
and a sample ATML ``nems.configure`` file (``nems.configure.atm_lnd.IN``) is available 
`here <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/parm/nems.configure.atm_lnd.IN>`__.


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





