:orphan:

*******************************************************************
Sample ``nems.configure`` File for the ``HAFS`` WM Configuration
*******************************************************************

.. code-block:: console

        # EARTH #
        EARTH_component_list: ATM OCN MED

        # MED #
        MED_model:                      cmeps
        MED_petlist_bounds:             1340 1399
        MED_attributes::
          coupling_mode = hafs
          system_type = ufs
          normalization = none
          merge_type = copy
          ATM_model = fv3
          OCN_model = hycom
          history_ymd = -999
          ScalarFieldCount = 0
          ScalarFieldIdxGridNX = 0
          ScalarFieldIdxGridNY = 0
          ScalarFieldName = cpl_scalars
        ::

        # ATM #
        ATM_model:                      fv3
        ATM_petlist_bounds:             0000 1339
        ATM_attributes::
          Verbosity = 1
          Diagnostic = 0
        ::

        # OCN #
        OCN_model:                      hycom
        OCN_petlist_bounds:             1340 1399
        OCN_attributes::
          Verbosity = 1
          Diagnostic = 0
          cdf_impexp_freq = 3
          cpl_hour = 0
          cpl_min = 0
          cpl_sec = 360
          base_dtg = 2020082512
          merge_import = .true.
          skip_first_import = .true.
          hycom_arche_output = .false.
          hyc_esmf_exp_output = .true.
          hyc_esmf_imp_output = .true.
          import_diagnostics = .false.
          import_setting = flexible
          hyc_impexp_file = nems.configure
          espc_show_impexp_minmax = .true.
          ocean_start_dtg = 43702.50000
          start_hour = 0
          start_min = 0
          start_sec = 0
          end_hour = 12
          end_min = 0
          end_sec = 0
        ::

        DRIVER_attributes::
          start_type = startup
        ::

        ALLCOMP_attributes::
          mediator_read_restart = false
        ::

        # CMEPS cold run sequence

        runSeq::
        @360
          ATM -> MED :remapMethod=redist
          MED med_phases_post_atm
          OCN -> MED :remapMethod=redist
          MED med_phases_post_ocn
          MED med_phases_prep_atm
          MED med_phases_prep_ocn_accum
          MED med_phases_prep_ocn_avg
          MED -> ATM :remapMethod=redist
          MED -> OCN :remapMethod=redist
          ATM
          OCN
        @
        ::

        # HYCOM field coupling configuration (location set by hyc_impexp_file)

        ocn_export_fields::
          'sst'     'sea_surface_temperature'   'K'
          'mask'    'ocean_mask'                '1'
        ::

        ocn_import_fields::
          'taux10'  'mean_zonal_moment_flx_atm' 'N_m-2'
          'tauy10'  'mean_merid_moment_flx_atm' 'N_m-2'
          'prcp'    'mean_prec_rate'            'kg_m-2_s-1'
          'swflxd'  'mean_net_sw_flx'           'W_m-2'
          'lwflxd'  'mean_net_lw_flx'           'W_m-2'
          'mslprs'  'inst_pres_height_surface'  'Pa'
          'sensflx' 'mean_sensi_heat_flx'       'W_m-2'
          'latflx'  'mean_laten_heat_flx'       'W_m-2'
        ::


