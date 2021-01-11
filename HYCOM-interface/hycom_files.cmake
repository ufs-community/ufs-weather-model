list(APPEND _hycom_src_files
  HYCOM/mod_dimensions.F90
  HYCOM/mod_xc.F90
  HYCOM/mod_za.F90
  HYCOM/mod_cb_arrays.F90
  HYCOM/mod_pipe.F90
  HYCOM/mod_incupd.F90
  HYCOM/mod_floats.F90
  HYCOM/mod_stokes.F90
  HYCOM/mod_tides.F90
  HYCOM/mod_mean.F90
  HYCOM/mod_archiv.F90
  HYCOM/mod_tsadvc.F90
  HYCOM/mod_momtum.F90
  HYCOM/mod_barotp.F90
  HYCOM/mod_asselin.F90
  HYCOM/mod_restart.F90
  HYCOM/mod_hycom.F90

  HYCOM/bigrid.F90
  HYCOM/blkdat.F90
  HYCOM/cnuity.F90
  HYCOM/convec.F90
  HYCOM/diapfl.F90
  HYCOM/dpthuv.F90
  HYCOM/dpudpv.F90
  HYCOM/forfun.F90
  HYCOM/geopar.F90
  HYCOM/hybgen.F90
  HYCOM/icloan.F90
  HYCOM/inicon.F90
  HYCOM/inigiss.F90
  HYCOM/inikpp.F90
  HYCOM/inimy.F90
  HYCOM/latbdy.F90
  HYCOM/matinv.F90
  HYCOM/mxkprf.F90
  HYCOM/mxkrt.F90
  HYCOM/mxkrtm.F90
  HYCOM/mxpwp.F90
  HYCOM/overtn.F90
  HYCOM/poflat.F90
  HYCOM/prtmsk.F90
  HYCOM/psmoo.F90
  HYCOM/thermf.F90
  HYCOM/trcupd.F90
  HYCOM/machine.F90
  HYCOM/wtime.F90
  HYCOM/machi_c.c
  HYCOM/isnan.F90
  HYCOM/s8gefs.F90
)

list(APPEND _hycom_nuopc_src_files
  HYCOM/NUOPC/HYCOM_OceanComp.F90
  HYCOM/NUOPC/HYCOM_ESMF_Extensions.F90
  HYCOM/NUOPC/hycom_couple.F90
  HYCOM/NUOPC/read_impexp_config_mod.F90
  HYCOM/NUOPC/impexpField_cdf_mod.F90
  HYCOM/NUOPC/export_from_hycom_tiled.F90
  HYCOM/NUOPC/hycom_read_latlon.F90
  HYCOM/NUOPC/ocn_comp_NUOPC.F90
)

list(APPEND _hycom_offline_src_files
  HYCOM/hycom.F90
)
