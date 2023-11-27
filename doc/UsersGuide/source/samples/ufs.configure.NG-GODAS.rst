:orphan:

***********************************************************************
Sample ``nems.configure`` File for the ``NG-GODAS`` WM Configuration
***********************************************************************

.. code-block:: console

	# EARTH #
	EARTH_component_list: MED ATM OCN ICE
	EARTH_attributes::
	  Verbosity = 0
	::

	# MED #
	MED_model:                      cmeps
	MED_petlist_bounds:             0 11
	  Verbosity = 5
	  dbug_flag = 5

	::

	# ATM #
	ATM_model:                      datm
	ATM_petlist_bounds:             0 11
	ATM_attributes::
	  Verbosity = 0
	  DumpFields = false
	  mesh_atm  = DATM_INPUT/cfsr_mesh.nc
	  diro = "."
	  logfile = atm.log
     stop_n = 24
     stop_option = nhours
     stop_ymd = -999
     write_restart_at_endofrun = .true.
	::

	# OCN #
	OCN_model:                      mom6
	OCN_petlist_bounds:             12 27
	OCN_attributes::
	  Verbosity = 0
	  DumpFields = false
	  ProfileMemory = false
	  OverwriteSlice = true
	  mesh_ocn = mesh.mx100.nc
	::

	# ICE #
	ICE_model:                      cice6
	ICE_petlist_bounds:             28 39
	ICE_attributes::
	  Verbosity = 0
	  DumpFields = false
	  ProfileMemory = false
	  OverwriteSlice = true
	  mesh_ice = mesh.mx100.nc
	  stop_n = 12
	  stop_option = nhours
	  stop_ymd = -999
	::

	# CMEPS concurrent warm run sequence

	runSeq::
	@3600
	   MED med_phases_prep_ocn_avg
	   MED -> OCN :remapMethod=redist
	   OCN
	   @900
	     MED med_phases_prep_ice
	     MED -> ICE :remapMethod=redist
	     ATM
	     ICE
	     ATM -> MED :remapMethod=redist
	     MED med_phases_post_atm
	     ICE -> MED :remapMethod=redist
	     MED med_phases_post_ice
	     MED med_phases_aofluxes_run
	     MED med_phases_prep_ocn_accum
	   @
	   OCN -> MED :remapMethod=redist
	   MED med_phases_post_ocn
	   MED med_phases_restart_write
	@
	::

	# CMEPS variables

	DRIVER_attributes::
	      mediator_read_restart = false
	::
	MED_attributes::
	      ATM_model = datm
	      ICE_model = cice6
	      OCN_model = mom6
	      history_n = 1
	      history_option = nhours
	      history_ymd = -999
	      coupling_mode = nems_orig_data
	::
	ALLCOMP_attributes::
	      ScalarFieldCount = 3
	      ScalarFieldIdxGridNX = 1
	      ScalarFieldIdxGridNY = 2
	      ScalarFieldIdxNextSwCday = 3
	      ScalarFieldName = cpl_scalars
	      start_type = startup
	      restart_dir = RESTART/
	      case_name = DATM_CFSR
	      restart_n = 12
	      restart_option = nhours
	      restart_ymd = -999
	      dbug_flag = 0
	      use_coldstart = false
	      use_mommesh = true
	      coldair_outbreak_mod = .false.
	      flds_wiso = .false.
	      flux_convergence = 0.0
	      flux_max_iteration = 2
	      ocn_surface_flux_scheme = 0
	      orb_eccen = 1.e36
	      orb_iyear = 2000
	      orb_iyear_align = 2000
	      orb_mode = fixed_year
	      orb_mvelp = 1.e36
	      orb_obliq = 1.e36
	::



