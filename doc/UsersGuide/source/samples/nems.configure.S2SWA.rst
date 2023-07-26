:orphan:

*******************************************************************
Sample ``nems.configure`` File for the ``S2SWA`` WM Configuration
*******************************************************************

.. code-block:: console

	# EARTH #
	EARTH_component_list: MED ATM CHM OCN ICE WAV
	EARTH_attributes::
	  Verbosity = 0
	::

	# MED #
	MED_model:                      cmeps
	MED_petlist_bounds:             0 287
	::

	# ATM #
	ATM_model:                      fv3
	ATM_petlist_bounds:             0 311
	ATM_attributes::
	  Verbosity = 0
	  DumpFields = false
	  ProfileMemory = false
	  OverwriteSlice = true
	::

	# CHM #
	CHM_model:                      gocart
	CHM_petlist_bounds:             0 287
	CHM_attributes::
	  Verbosity = 0
	::

	# OCN #
	OCN_model:                      mom6
	OCN_petlist_bounds:             312 431
	OCN_attributes::
	  Verbosity = 0
	  DumpFields = false
	  ProfileMemory = false
	  OverwriteSlice = true
	  mesh_ocn = mesh.mx025.nc
	::

	# ICE #
	ICE_model:                      cice6
	ICE_petlist_bounds:             432 479
	ICE_attributes::
	  Verbosity = 0
	  DumpFields = false
	  ProfileMemory = false
	  OverwriteSlice = true
	  mesh_ice = mesh.mx025.nc
	  stop_n = 6
	  stop_option = nhours
	  stop_ymd = -999
	::

	# WAV #
	WAV_model:                      ww3
	WAV_petlist_bounds:             480 559
	WAV_attributes::
	  Verbosity = 0
	  OverwriteSlice = false
	::

	# CMEPS warm run sequence
	runSeq::
	@1800
	   MED med_phases_prep_ocn_avg
	   MED -> OCN :remapMethod=redist
	   OCN -> WAV
	   WAV -> OCN :srcMaskValues=1
	   OCN
	   @300
	     MED med_phases_prep_atm
	     MED med_phases_prep_ice
	     MED -> ATM :remapMethod=redist
	     MED -> ICE :remapMethod=redist
	     WAV -> ATM :srcMaskValues=1
	     ATM -> WAV
	     ICE -> WAV
	     ATM phase1
	     ATM -> CHM
	     CHM
	     CHM -> ATM
	     ATM phase2
	     ICE
	     WAV
	     ATM -> MED :remapMethod=redist
	     MED med_phases_post_atm
	     ICE -> MED :remapMethod=redist
	     MED med_phases_post_ice
	     MED med_phases_prep_ocn_accum
	   @
	   OCN -> MED :remapMethod=redist
	   MED med_phases_post_ocn
	   MED med_phases_restart_write
	@
	::

	# CMEPS variables

	DRIVER_attributes::
	::

	MED_attributes::
	      ATM_model = fv3
	      ICE_model = cice6
	      OCN_model = mom6
	      history_n = 1
	      history_option = nhours
	      history_ymd = -999
	      coupling_mode = nems_frac
	      history_tile_atm = 384
	::
	ALLCOMP_attributes::
	      ScalarFieldCount = 2
	      ScalarFieldIdxGridNX = 1
	      ScalarFieldIdxGridNY = 2
	      ScalarFieldName = cpl_scalars
	      start_type = startup
	      restart_dir = RESTART/
	      case_name = ufs.cpld
	      restart_n = 6
	      restart_option = nhours
	      restart_ymd = -999
	      dbug_flag = 0
	      use_coldstart = false
	      use_mommesh = true
	      eps_imesh = 1.0e-1
	      stop_n = 6
	      stop_option = nhours
	      stop_ymd = -999
	::


