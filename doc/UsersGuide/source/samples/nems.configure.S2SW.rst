:orphan:

*******************************************************************
Sample ``nems.configure`` File for the ``S2SW`` WM Configuration
*******************************************************************

.. code-block:: console

	# EARTH #
	EARTH_component_list: MED ATM OCN ICE WAV
	EARTH_attributes::
	  Verbosity = 0
	::

	# MED #
	MED_model:                      cmeps
	MED_petlist_bounds:             0 143
	::

	# ATM #
	ATM_model:                      fv3
	ATM_petlist_bounds:             0 149
	ATM_attributes::
	::

	# OCN #
	OCN_model:                      mom6
	OCN_petlist_bounds:             150 179
	OCN_attributes::
	  mesh_ocn = mesh.mx100.nc
	::

	# ICE #
	ICE_model:                      cice6
	ICE_petlist_bounds:             180 191
	ICE_attributes::
	  mesh_ice = mesh.mx100.nc
	::

	# WAV #
	WAV_model:                      ww3
	WAV_petlist_bounds:             192 395
	WAV_attributes::
	::

	# CMEPS warm run sequence
	runSeq::
	@3600
	   MED med_phases_prep_ocn_avg
	   MED -> OCN :remapMethod=redist
	   OCN -> WAV
	   WAV -> OCN :srcMaskValues=1
	   OCN
	   @900
	     MED med_phases_prep_atm
	     MED med_phases_prep_ice
	     MED -> ATM :remapMethod=redist
	     MED -> ICE :remapMethod=redist
	     WAV -> ATM :srcMaskValues=1
	     ATM -> WAV
	     ICE -> WAV
	     ATM
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

	::
	MED_attributes::
	      ATM_model = fv3
	      ICE_model = cice6
	      OCN_model = mom6
	      history_n = 1
	      history_option = nhours
	      history_ymd = -999
	      coupling_mode = nems_orig
	::
	ALLCOMP_attributes::
	      ScalarFieldCount = 2
	      ScalarFieldIdxGridNX = 1
	      ScalarFieldIdxGridNY = 2
	      ScalarFieldName = cpl_scalars
	      start_type = startup
	      restart_dir = RESTART/
	      case_name = ufs.cpld
	      restart_n = 24
	      restart_option = nhours
	      restart_ymd = -999
	      dbug_flag = 0
	      use_coldstart = false
	      use_mommesh = true
	::
