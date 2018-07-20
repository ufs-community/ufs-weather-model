module physics_abstraction_layer

  use GFS_typedefs,    only: init_type         =>  GFS_init_type,       &
                             control_type      =>  GFS_control_type,    &
                             statein_type      =>  GFS_statein_type,    &
                             stateout_type     =>  GFS_stateout_type,   &
                             sfcprop_type      =>  GFS_sfcprop_type,    &
                             coupling_type     =>  GFS_coupling_type,   &
                             grid_type         =>  GFS_grid_type,       &
                             tbd_type          =>  GFS_tbd_type,        &
                             cldprop_type      =>  GFS_cldprop_type,    &
                             radtend_type      =>  GFS_radtend_type,    &
                             intdiag_type      =>  GFS_diag_type,       &
                             interstitial_type =>  GFS_interstitial_type

#ifdef CCPP
  use GFS_driver,      only: initialize       =>  GFS_initialize,       &
                             time_vary_step   =>  GFS_time_vary_step,   &
                             finalize         =>  GFS_finalize
#else
  use GFS_driver,      only: initialize       =>  GFS_initialize,       &
                             time_vary_step   =>  GFS_time_vary_step,   &
                             radiation_step1  =>  GFS_radiation_driver, &
                             physics_step1    =>  GFS_physics_driver,   &
                             physics_step2    =>  GFS_stochastic_driver,&
                             finalize         =>  GFS_finalize
#endif

!----------------------
!  public physics types
!----------------------
  public  init_type
  public  control_type
  public  statein_type
  public  stateout_type
  public  sfcprop_type
  public  coupling_type
  public  grid_type
  public  tbd_type
  public  cldprop_type
  public  radtend_type
  public  intdiag_type
  public  interstitial_type

!--------------------------
!  public physics functions
!--------------------------
  public  initialize
#ifndef CCPP
  public  time_vary_step
  public  radiation_step1
  public  physics_step1
  public  physics_step2
#endif
  public  finalize

CONTAINS

end module physics_abstraction_layer
