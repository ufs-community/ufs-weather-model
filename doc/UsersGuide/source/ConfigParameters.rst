.. _ConfigParams:

******************************************
Configuration Parameters
******************************************

=================================
Build Configuration Parameters
=================================

.. _dapp:

Configuration Options
=========================

``-DAPP``:
   Sets the :term:`WM` configuration to build. 
   Valid values: ``ATM``, ``ATMW``, ``ATMAERO``, ``ATMAQ``, ``S2S``, ``S2SA``, ``S2SW``, ``S2SWA``, ``NG-GODAS``, ``HAFS``, ``HAFSW``, ``HAFS-ALL``


.. _suites:

Physics Options
=======================

``-DCCPP_SUITES``:
   Sets the physics suites that will be made available when the :term:`WM` is built. 
   
   Physics suites supported in regression testing:
   
   | ``FV3_GFS_cpld_rasmgshocnsstnoahmp_ugwp``
   | ``FV3_GFS_v15p2``
   | ``FV3_GFS_v15_thompson_mynn``
   | ``FV3_GFS_v15_thompson_mynn_lam3km``
   | ``FV3_GFS_v16``
   | ``FV3_GFS_v16_csawmg``
   | ``FV3_GFS_v16_fv3wam``
   | ``FV3_GFS_v16_noahmp``
   | ``FV3_GFS_v16_ras``
   | ``FV3_GFS_v16_ugwpv1``
   | ``FV3_GFS_v17_p8``
   | ``FV3_GFS_v17_p8_rrtmgp``
   | ``FV3_GFS_v17_coupled_p8``
   | ``FV3_GFS_v17_coupled_p8_sfcocn``
   | ``FV3_HAFS_v0_gfdlmp_tedmf``
   | ``FV3_HAFS_v0_gfdlmp_tedmf_nonsst``
   | ``FV3_HAFS_v0_thompson_tedmf_gfdlsf``
   | ``FV3_HRRR``
   | ``FV3_HRRR_smoke``
   | ``FV3_RAP``
   | ``FV3_RAP_RRTMGP``
   | ``FV3_RAP_sfcdiff``
   | ``FV3_RRFS_v1beta``
   | ``FV3_RRFS_v1nssl``

   Other valid values: 

   | ``FV3_CPT_v0``
   | ``FV3_GFS_2017``
   | ``FV3_GFS_2017_csawmg``
   | ``FV3_GFS_2017_csawmgshoc``
   | ``FV3_GFS_2017_gfdlmp``
   | ``FV3_GFS_2017_gfdlmp_noahmp``
   | ``FV3_GFS_2017_gfdlmp_regional``
   | ``FV3_GFS_2017_gfdlmp_regional_c768``
   | ``FV3_GFS_2017_h2ophys``
   | ``FV3_GFS_2017_myj``
   | ``FV3_GFS_2017_ntiedtke``
   | ``FV3_GFS_2017_ozphys_2015``
   | ``FV3_GFS_2017_sas``
   | ``FV3_GFS_2017_satmedmf``
   | ``FV3_GFS_2017_satmedmfq``
   | ``FV3_GFS_2017_shinhong``
   | ``FV3_GFS_2017_stretched``
   | ``FV3_GFS_2017_ysu``
   | ``FV3_GFS_cpld_rasmgshoc``
   | ``FV3_GFS_cpld_rasmgshocnsst``
   | ``FV3_GFS_cpld_rasmgshocnsst_flake``
   | ``FV3_GFS_cpld_rasmgshocnsst_ugwp``
   | ``FV3_GFS_cpldnst_rasmgshoc``
   | ``FV3_GFS_rasmgshoc``
   | ``FV3_GFS_v15``
   | ``FV3_GFS_v15_gf``
   | ``FV3_GFS_v15_gf_thompson``
   | ``FV3_GFS_v15_mynn``
   | ``FV3_GFS_v15_ras``
   | ``FV3_GFS_v15_rasmgshoc``
   | ``FV3_GFS_v15_thompson``
   | ``FV3_GFS_v15p2_no_nsst``
   | ``FV3_GFS_v15plus``
   | ``FV3_GFS_v15plusras``
   | ``FV3_GFS_v16_coupled``
   | ``FV3_GFS_v16_coupled_noahmp``
   | ``FV3_GFS_v16_coupled_nsstNoahmp``
   | ``FV3_GFS_v16_coupled_nsstNoahmpUGWPv1``
   | ``FV3_GFS_v16_coupled_p8``
   | ``FV3_GFS_v16_coupled_p8_sfcocn``
   | ``FV3_GFS_v16_couplednsst``
   | ``FV3_GFS_v16_flake``
   | ``FV3_GFS_v16_no_nsst``
   | ``FV3_GFS_v16_nsstNoahmpUGWPv1``
   | ``FV3_GFS_v16_p8``
   | ``FV3_GFS_v16_thompson``
   | ``FV3_GFSv17alp_cpldnsstrasnoahmp``
   | ``FV3_GFSv17alp_cpldnsstrasugwpnoahmp``
   | ``FV3_GFSv17alp_cpldnsstsasugwpnoahmp``
   | ``FV3_GFSv17alpha_cpldnsstras``
   | ``FV3_GFSv17alpha_cpldnsstras_flake``
   | ``FV3_GFSv17alpha_cpldnsstras_ugwp``
   | ``FV3_GFSv17alpha_cpldnsstrasnoshal``
   | ``FV3_GFSv17alpha_cpldnsstsas``
   | ``FV3_GFSv17alpha_cpldnsstsas_ugwp``
   | ``FV3_GFSv17alpha_ras``
   | ``FV3_GFSv17alpha_ras_flake``
   | ``FV3_GFSv17alpha_ras_ugwp``
   | ``FV3_GFSv17alpha_sas``
   | ``FV3_RAP_cires_ugwp``
   | ``FV3_RAP_flake``
   | ``FV3_RAP_noah``
   | ``FV3_RAP_noah_sfcdiff_cires_ugwp``
   | ``FV3_RAP_noah_sfcdiff_ugwpv1``
   | ``FV3_RAP_noah_sfcdiff_unified_ugwp``
   | ``FV3_RAP_unified_ugwp``
   | ``FV3_RRFS_v1alpha``

.. _other-build-options:

Other Build Options
=======================

``-DCMEPS_AOFLUX``: (Default: OFF)
   Enables atmosphere-ocean flux calculation in mediator. 
   Valid values: ``ON`` | ``OFF``

   .. COMMENT: But when/why would you do this?

``-DDEBUG``: (Default: OFF)
   Enables DEBUG mode.
   Valid values: ``ON`` | ``OFF``

   .. COMMENT: And what extras does DEBUG mode provide (that VERBOSE) doesn't?

``-D32BIT``: (Default: OFF)
   Enables 32-bit, single precision arithmetic in dycore and fast physics.
   Valid values: ``ON`` | ``OFF``

   .. COMMENT: But when/why would you do this?

``-DCCPP_32BIT``: (Default: OFF)
   Enables 32-bit, single precision arithmetic in slow physics.
   Valid values: ``ON`` | ``OFF``

   .. COMMENT: But when/why would you do this?

``-DMOVING_NEST``: (Default: OFF)
   Enables moving nest code.
   Valid values: ``ON`` | ``OFF``

   .. COMMENT: But what does that mean? When/why is the moving nest used?

``-DMULTI_GASES``: (Default: OFF)
   Enable ``MULTI_GASES``. 
   Valid values: ``ON`` | ``OFF``

   .. COMMENT: But what does this DO?! And when/why is it used?


.. COMMENT: Add any of the following options with -D in front???
      set(AVX2            ON  CACHE BOOL "Enable AVX2 instruction set")
      set(AVX             OFF CACHE BOOL "Enable AVX-I instruction set")
      set(SIMDMULTIARCH   OFF CACHE BOOL "Enable multi-target SIMD instruction sets")
      set(INLINE_POST     OFF CACHE BOOL "Enable inline post")
      set(OPENMP          ON  CACHE BOOL "Enable OpenMP threading")
      set(PARALLEL_NETCDF OFF CACHE BOOL "Enable parallel NetCDF")
      set(JEDI_DRIVER     OFF CACHE BOOL "Enable JEDI as top level driver")


