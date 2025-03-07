.. _ConfigParams:

******************************************
Configuration Parameters
******************************************

The UFS Weather Model build parameters are set in :wm-repo:`CMakeLists.txt <blob/develop/CMakeLists.txt>` or in the ``CMakeLists.txt`` file of one of its subrepositories. 

=================================
Build Configuration Parameters
=================================

.. _dapp:

Configuration Options
=========================

``-DAPP``:
   Sets the :term:`WM` configuration to build. 
   Valid values: ``ATM``, ``ATMW``, ``ATMAERO``, ``ATMAQ``, ``ATMWM``, ``ATML``, ``ATMF``, ``ATM_DS2S``, ``ATM_DS2S-PCICE``, ``LND``, ``LND-LM4``, ``S2S``, ``S2SA``, ``S2SW``, ``S2SWA``, ``S2SL``, ``S2SWL``, ``S2SWAL``, ``NG-GODAS``, ``HAFS``, ``HAFSW``, ``HAFS-MOM6``, ``HAFS-MOM6W``, ``HAFS-ALL``


.. _suites:

Physics Options
=======================

``-DCCPP_SUITES``:
   Sets the physics suites that will be made available when the :term:`WM` is built. 
   
   Physics suites supported in regression testing:
   
   | ``FV3_GFS_v17_coupled_p8``
   | ``FV3_GFS_v17_coupled_p8_sfcocn``
   | ``FV3_GFS_v17_coupled_p8_ugwpv1`` 
   | ``FV3_GFS_v17_p8``
   | ``FV3_GFS_v17_p8_mynn``
   | ``FV3_GFS_v17_p8_rrtmgp``
   | ``FV3_GFS_v17_p8_ugwpv1``
   | ``FV3_GFS_v16``
   | ``FV3_GFS_v16_csawmg``
   | ``FV3_GFS_v16_flake``
   | ``FV3_GFS_v16_ras``
   | ``FV3_GFS_v15p2``
   | ``FV3_GFS_v15_thompson_mynn_lam3km``
   | ``FV3_global_nest_v1``
   | ``FV3_HAFS_v1_gfdlmp_tedmf``
   | ``FV3_HAFS_v1_gfdlmp_tedmf_nonsst``
   | ``FV3_HAFS_v1_thompson``
   | ``FV3_HAFS_v1_thompson_nonsst``
   | ``FV3_HAFS_v1_thompson_tedmf_gfdlsf``
   | ``FV3_HRRR``
   | ``FV3_HRRR_c3``
   | ``FV3_HRRR_gf``
   | ``FV3_RAP``
   | ``FV3_RAP_cires_ugwp``
   | ``FV3_RAP_clm_lake``
   | ``FV3_RAP_flake``
   | ``FV3_RAP_noah``
   | ``FV3_RAP_noah_sfcdiff_cires_ugwp``
   | ``FV3_RAP_sfcdiff``
   | ``FV3_RAP_unified_ugwp``
   | ``FV3_RRFS_v1beta``
   | ``FV3_RRFS_v1nssl``
   | ``FV3_WoFS_v0``

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


