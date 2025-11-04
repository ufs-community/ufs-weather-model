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
   Valid values: ``ATM``, ``ATMW``, ``ATMAERO``, ``ATMAQ``, ``ATML``, ``ATMF``, ``ATMMPAS``, ``ATM_DS2S``, ``ATM_DS2S-PCICE``, ``LND``, ``LND-LM4``, ``S2S``, ``S2SA``, ``S2SW``, ``S2SWA``, ``S2SL``, ``S2SWL``, ``S2SWAL``, ``NG-GODAS``, ``HAFS``, ``HAFSW``, ``HAFS-MOM6``, ``HAFS-MOM6W``, ``HAFS-ALL``


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
   | ``FV3_GFS_v17_p8_ugwpv1_tempo``
   | ``FV3_GFS_v16``
   | ``FV3_GFS_v16_csawmg``
   | ``FV3_GFS_v16_flake``
   | ``FV3_GFS_v16_ras``
   | ``FV3_GFS_v15_thompson_mynn_lam3km``
   | ``FV3_global_nest_v1``
   | ``FV3_HAFS_v1_gfdlmp_tedmf``
   | ``FV3_HAFS_v1_gfdlmpv3_tedmf``
   | ``FV3_HAFS_v1_gfdlmp_tedmf_nonsst``
   | ``FV3_HAFS_v1_thompson``
   | ``FV3_HAFS_v1_thompson_nonsst``
   | ``FV3_HAFS_v1_thompson_tedmf_gfdlsf``
   | ``FV3_HRRR``
   | ``FV3_HRRR_c3``
   | ``FV3_HRRR_gf``
   | ``FV3_ideal_mp_nssl``
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
   | ``MPAS_RRFS``

.. _other-build-options:

Other Build Options
=======================

``-D32BIT`` (Default: OFF)  
   Enable 32BIT (single precision arithmetic in dycore and fast physics).  
   Valid values: ``ON`` | ``OFF``  

``-DCCPP_32BIT`` (Default: OFF)  
   Enable CCPP_32BIT (single precision arithmetic in dycore and slow physics).  
   Valid values: ``ON`` | ``OFF``  

``-DAVX2`` (Default: ON)  
   Enable AVX2 instruction set.  
   Valid values: ``ON`` | ``OFF``  

``-DAVX`` (Default: OFF)  
   Enable AVX-I instruction set.  
   Valid values: ``ON`` | ``OFF``  

``-DSIMDMULTIARCH`` (Default: OFF)  
   Enable multi-target SIMD instruction sets.  
   Valid values: ``ON`` | ``OFF``  

``-DDEBUG`` (Default: OFF)  
   Enable DEBUG mode.  
   Valid values: ``ON`` | ``OFF``  

``-DDISABLE_FMA`` (Default: OFF)  
   Disable Fused Multiply-Add instructions (workaround needed for AMD EPYC).  
   Valid values: ``ON`` | ``OFF``  

``-DINLINE_POST`` (Default: ON)  
   Enable inline post.  
   Valid values: ``ON`` | ``OFF``  

``-DMULTI_GASES`` (Default: OFF)  
   Enable MULTI_GASES.  
   Valid values: ``ON`` | ``OFF``  

``-DMOVING_NEST`` (Default: OFF)  
   Enable moving nest code.  
   Valid values: ``ON`` | ``OFF``  

``-DREGIONAL_MOM6`` (Default: OFF)  
   Enable Regional MOM6.  
   Valid values: ``ON`` | ``OFF``  

``-DOPENMP`` (Default: ON)  
   Enable OpenMP threading.  
   Valid values: ``ON`` | ``OFF``  

``-DPARALLEL_NETCDF`` (Default: OFF)  
   Enable parallel NetCDF.  
   Valid values: ``ON`` | ``OFF``  

``-DJEDI_DRIVER`` (Default: OFF)  
   Enable JEDI as top level driver.  
   Valid values: ``ON`` | ``OFF``  

``-DCMEPS_AOFLUX`` (Default: OFF)  
   Enable atmosphere-ocean flux calculation in mediator.  
   Valid values: ``ON`` | ``OFF``  

``-DPDLIB`` (Default: OFF)  
   Enable domain decomposition in WW3 via PDLIB with BT1.  
   Valid values: ``ON`` | ``OFF``  

``-DPDLIB_BT4`` (Default: OFF)  
   Enable domain decomposition in WW3 via PDLIB with BT4.  
   Valid values: ``ON`` | ``OFF``  

``-DCDEPS_INLINE`` (Default: OFF)  
   Enable CDEPS inline capability.  
   Valid values: ``ON`` | ``OFF``  

``-DHYDRO`` (Default: OFF)  
   Enable hydrostatic set.  
   Valid values: ``ON`` | ``OFF``  

``BUILD_WITH_IFI`` (Default: OFF)  
   Build UPP with In-Flight Icing (IFI) library if present.  
   Valid values: ``ON`` | ``OFF``  

``REQUIRE_IFI`` (Default: OFF)  
   Abort if libIFI is not found ; enables BUILD_WITH_IFI=ON.  
   Valid values: ``ON`` | ``OFF``  

``INTERNAL_IFI`` (Default: OFF)  
   Compile with IFI inside the executable, instead of a library.  
   Valid values: ``ON`` | ``OFF``  

