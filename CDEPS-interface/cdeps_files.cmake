# CDEPS generates these with genf90.pl acting on templates .in in CDEPS/share
# For the UFS, they are pre-generated and committed to the UFS repository
list(APPEND ufs_cdeps_share_files
  ufs/cdeps_share/shr_assert_mod.F90
  ufs/cdeps_share/shr_frz_mod.F90
  ufs/cdeps_share/shr_infnan_mod.F90
)

list(APPEND cdeps_share_files
  CDEPS/share/glc_elevclass_mod.F90
  CDEPS/share/shr_abort_mod.F90
  CDEPS/share/shr_assert.h
  CDEPS/share/shr_cal_mod.F90
  CDEPS/share/shr_const_mod.F90
  CDEPS/share/shr_file_mod.F90
  CDEPS/share/shr_kind_mod.F90
  CDEPS/share/shr_log_mod.F90
  CDEPS/share/shr_nl_mod.F90
  CDEPS/share/shr_orb_mod.F90
  CDEPS/share/shr_precip_mod.F90
  CDEPS/share/shr_strconvert_mod.F90
  CDEPS/share/shr_string_mod.F90
  CDEPS/share/shr_sys_mod.F90
  CDEPS/share/shr_timer_mod.F90
  CDEPS/share/shr_file_mod.F90
  CDEPS/share/shr_nl_mod.F90
)

list(APPEND cdeps_streams_files
  CDEPS/streams/dshr_methods_mod.F90
  CDEPS/streams/dshr_strdata_mod.F90
  CDEPS/streams/dshr_stream_mod.F90
  CDEPS/streams/dshr_tinterp_mod.F90
)

list(APPEND cdeps_dshr_files
  CDEPS/dshr/dshr_dfield_mod.F90
  CDEPS/dshr/dshr_fldlist_mod.F90
  CDEPS/dshr/dshr_mod.F90
)

list(APPEND cdeps_datm_files
  CDEPS/datm/atm_comp_nuopc.F90
  CDEPS/datm/datm_datamode_cfsr_mod.F90
  CDEPS/datm/datm_datamode_clmncep_mod.F90
  CDEPS/datm/datm_datamode_core2_mod.F90
  CDEPS/datm/datm_datamode_cplhist_mod.F90
  CDEPS/datm/datm_datamode_era5_mod.F90
  CDEPS/datm/datm_datamode_gefs_mod.F90
  CDEPS/datm/datm_datamode_gfs_mod.F90
  CDEPS/datm/datm_datamode_gfs_hafs_mod.F90
  CDEPS/datm/datm_datamode_jra_mod.F90
  CDEPS/datm/datm_datamode_simple_mod.F90
)

list(APPEND cdeps_dice_files
  CDEPS/dice/dice_datamode_ssmi_mod.F90
  CDEPS/dice/dice_flux_atmice_mod.F90
  CDEPS/dice/ice_comp_nuopc.F90
)

list(APPEND cdeps_dlnd_files
  CDEPS/dlnd/lnd_comp_nuopc.F90
)

list(APPEND cdeps_docn_files
  CDEPS/docn/docn_datamode_aquaplanet_mod.F90
  CDEPS/docn/docn_datamode_copyall_mod.F90
  CDEPS/docn/docn_datamode_iaf_mod.F90
  CDEPS/docn/docn_datamode_som_mod.F90
  CDEPS/docn/docn_datamode_cplhist_mod.F90
  CDEPS/docn/docn_import_data_mod.F90
  CDEPS/docn/ocn_comp_nuopc.F90
)

list(APPEND cdeps_drof_files
  CDEPS/drof/rof_comp_nuopc.F90
)

list(APPEND cdeps_dwav_files
  CDEPS/dwav/wav_comp_nuopc.F90
)
