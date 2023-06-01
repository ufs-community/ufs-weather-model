:orphan:

*******************************************************************
Sample ``nems.configure`` File for the ``S2S`` WM Configuration
*******************************************************************

.. code-block:: console

        # EARTH #
        EARTH_component_list: MED ATM CHM OCN ICE WAV
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
          Verbosity = 0
          DumpFields = false
          ProfileMemory = false
          OverwriteSlice = true
        ::
        
        # OCN #
        OCN_model:                      mom6
        OCN_petlist_bounds:             150 269
        OCN_attributes::
          Verbosity = 0
          DumpFields = false
          ProfileMemory = false
          OverwriteSlice = true
          mesh_ocn = mesh.mx025.nc
        ::
        
        # ICE #
        ICE_model:                      cice6
        ICE_petlist_bounds:             270 317
        ICE_attributes::
          Verbosity = 0
          DumpFields = false
          ProfileMemory = false
          OverwriteSlice = true
          mesh_ice = mesh.mx025.nc
          stop_n = 840
          stop_option = nhours
          stop_ymd = -999
        ::
        
        # CMEPS warm run sequence
        runSeq::
        @720
           MED med_phases_prep_ocn_avg
           MED -> OCN :remapMethod=redist
           OCN
           @720
             MED med_phases_aofluxes_run
             MED med_phases_prep_atm
             MED med_phases_prep_ice
             MED -> ATM :remapMethod=redist
             MED -> ICE :remapMethod=redist
             ATM
             ICE
             ATM -> MED :remapMethod=redist
             MED med_phases_post_atm
             ICE -> MED :remapMethod=redist
             MED med_phases_post_ice
             MED med_phases_prep_ocn_accum
           @
           OCN -> MED :remapMethod=redist
           MED med_phases_post_ocn
           MED med_phases_restart_write
           MED med_phases_history_write
        @
        ::
        
        # CMEPS variables
        
        DRIVER_attributes::
        ::
        
        MED_attributes::
              ATM_model = fv3
              ICE_model = cice6
              OCN_model = mom6
              history_n = 3
              history_option = nhours
              history_ymd = -999
              coupling_mode = nems_frac_aoflux
              history_tile_atm = 96
              aoflux_grid = 'xgrid'
              aoflux_code = 'ccpp'
              aoflux_ccpp_suite = 'FV3_sfc_ocean'
              ccpp_restart_interval = -1
              ccpp_ini_mosaic_file = 'INPUT/C96_mosaic.nc'
              ccpp_input_dir = 'INPUT/'
              ccpp_ini_file_prefix = 'INPUT/sfc_data.tile'
              ccpp_nstf_name = 2,1,0,0,0
              ccpp_ini_read = true
        ::
        ALLCOMP_attributes::
              ScalarFieldCount = 2
              ScalarFieldIdxGridNX = 1
              ScalarFieldIdxGridNY = 2
              ScalarFieldName = cpl_scalars
              start_type = startup
              restart_dir = RESTART/
              case_name = ufs.cpld
              restart_n = 12
              restart_option = nhours
              restart_ymd = -999
              dbug_flag = 0
              use_coldstart = false
              use_mommesh = true
              eps_imesh = 1.0e-1
              stop_n = 840
              stop_option = nhours
              stop_ymd = -999
        ::

.. note:: The *aoflux_grid* option is used to select the grid/mesh to perform atmosphere-ocean flux calculation. The possible options are *xgrid* (exchange grid), *agrid* (atmosphere model grid) and *ogrid* (ocean model grid).

.. note:: The *aoflux_code* option is used to define the algorithm that will be used to calculate atmosphere-ocean fluxes. The possible options are *cesm* and *ccpp*. If *ccpp* is selected then the suite file provided in the *aoflux_ccpp_suite* option is used to calculate atmosphere-ocean fluxes through the use of CCPP host model.

