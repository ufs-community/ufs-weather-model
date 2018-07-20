module IPD_driver

  use IPD_typedefs,               only: IPD_init_type,                       &
                                        IPD_control_type,  IPD_data_type,    &
                                        IPD_diag_type,     IPD_restart_type, &
                                        IPD_interstitial_type
#ifdef CCPP
  use physics_abstraction_layer,  only: initialize,        time_vary_step,   &
                                        finalize
#else
  use physics_abstraction_layer,  only: initialize,        time_vary_step,   &
                                        radiation_step1,   physics_step1,    &
                                        physics_step2,     finalize
#endif

  use physics_diag_layer,         only: diag_populate

  use physics_restart_layer,      only: restart_populate

  implicit none

!------------------------------------------------------!
!  IPD containers                                      !
!------------------------------------------------------!
!  type(GFS_control_type)              :: IPD_Control  !
!  type(IPD_data_type)     allocatable :: IPD_Data(:)  !
!  type(IPD_diag_type),                :: IPD_Diag(:)  !
!  type(IPD_restart_type),             :: IPD_Restart  !
!------------------------------------------------------!

!----------------
! Public Entities
!----------------
! functions
  public IPD_initialize
  public IPD_setup_step
#ifndef CCPP
  public IPD_radiation_step
  public IPD_physics_step1
  public IPD_physics_step2
#endif
  public IPD_finalize

  CONTAINS
!*******************************************************************************************


!----------------
!  IPD initialize 
!----------------
  subroutine IPD_initialize (IPD_control, IPD_Data, IPD_Diag, IPD_Restart, IPD_Interstitial, IPD_init_parm)
    type(IPD_control_type), intent(inout)      :: IPD_Control
    type(IPD_data_type),    intent(inout)      :: IPD_Data(:)
    type(IPD_diag_type),    intent(inout)      :: IPD_Diag(:)
    type(IPD_restart_type), intent(inout)      :: IPD_Restart
    type(IPD_interstitial_type), intent(inout) :: IPD_Interstitial(:)
    type(IPD_init_type),    intent(in)         :: IPD_init_parm

    !--- initialize the physics suite
    call initialize (IPD_Control, IPD_Data(:)%Statein, IPD_Data(:)%Stateout,      &
                     IPD_Data(:)%Sfcprop, IPD_Data(:)%Coupling, IPD_Data(:)%Grid, &
                     IPD_Data(:)%Tbd, IPD_Data(:)%Cldprop, IPD_Data(:)%Radtend,   &
                     IPD_Data(:)%Intdiag, IPD_Interstitial(:), IPD_init_parm)


    !--- populate/associate the Diag container elements
    call diag_populate (IPD_Diag(:), IPD_control, IPD_Data%Statein, IPD_Data(:)%Stateout, &
                        IPD_Data(:)%Sfcprop, IPD_Data(:)%Coupling, IPD_Data(:)%Grid,      &
                        IPD_Data(:)%Tbd, IPD_Data(:)%Cldprop, IPD_Data(:)%Radtend,        &
                        IPD_Data(:)%Intdiag, IPD_init_parm)


    !--- allocate and populate/associate the Restart container elements
    call restart_populate (IPD_Restart, IPD_control, IPD_Data(:)%Statein, IPD_Data(:)%Stateout, &
                           IPD_Data(:)%Sfcprop, IPD_Data(:)%Coupling, IPD_Data(:)%Grid,         &
                           IPD_Data(:)%Tbd, IPD_Data(:)%Cldprop, IPD_Data(:)%Radtend,           &
                           IPD_Data(:)%Intdiag, IPD_init_parm)

  end subroutine IPD_initialize


!---------------------------------------------
!  IPD setup step
!    surface data cycling, random streams, etc
!---------------------------------------------
  subroutine IPD_setup_step (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart)
    type(IPD_control_type), intent(inout) :: IPD_Control
    type(IPD_data_type),    intent(inout) :: IPD_Data(:)
    type(IPD_diag_type),    intent(inout) :: IPD_Diag(:)
    type(IPD_restart_type), intent(inout) :: IPD_Restart

    call time_vary_step (IPD_Control, IPD_Data(:)%Statein, IPD_Data(:)%Stateout,      &
                         IPD_Data(:)%Sfcprop, IPD_Data(:)%Coupling, IPD_Data(:)%Grid, &
                         IPD_Data(:)%Tbd, IPD_Data(:)%Cldprop, IPD_Data(:)%Radtend,   &
                         IPD_Data(:)%Intdiag)

  end subroutine IPD_setup_step


#ifndef CCPP
!--------------------
!  IPD radiation step
!--------------------
  subroutine IPD_radiation_step (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart)
    type(IPD_control_type), intent(inout) :: IPD_Control
    type(IPD_data_type),    intent(inout) :: IPD_Data
    type(IPD_diag_type),    intent(inout) :: IPD_Diag(:)
    type(IPD_restart_type), intent(inout) :: IPD_Restart

    call radiation_step1 (IPD_control, IPD_Data%Statein, IPD_Data%Stateout,   &
                          IPD_Data%Sfcprop, IPD_Data%Coupling, IPD_Data%Grid, &
                          IPD_Data%Tbd, IPD_Data%Cldprop, IPD_Data%Radtend,   &
                          IPD_Data%Intdiag)

  end subroutine IPD_radiation_step


!-------------------
!  IPD physics step1
!-------------------
  subroutine IPD_physics_step1 (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart)
    type(IPD_control_type), intent(inout) :: IPD_Control
    type(IPD_data_type),    intent(inout) :: IPD_Data
    type(IPD_diag_type),    intent(inout) :: IPD_Diag(:)
    type(IPD_restart_type), intent(inout) :: IPD_Restart

    call physics_step1 (IPD_control, IPD_Data%Statein, IPD_Data%Stateout,   &
                        IPD_Data%Sfcprop, IPD_Data%Coupling, IPD_Data%Grid, &
                        IPD_Data%Tbd, IPD_Data%Cldprop, IPD_Data%Radtend,   &
                        IPD_Data%Intdiag)

  end subroutine IPD_physics_step1


!-------------------
!  IPD physics step2
!-------------------
  subroutine IPD_physics_step2 (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart)
    type(IPD_control_type), intent(inout) :: IPD_Control
    type(IPD_data_type),    intent(inout) :: IPD_Data
    type(IPD_diag_type),    intent(inout) :: IPD_Diag(:)
    type(IPD_restart_type), intent(inout) :: IPD_Restart

    call physics_step2 (IPD_control, IPD_Data%Statein, IPD_Data%Stateout,   &
                        IPD_Data%Sfcprop, IPD_Data%Coupling, IPD_Data%Grid, &
                        IPD_Data%Tbd, IPD_Data%Cldprop, IPD_Data%Radtend,   &
                        IPD_Data%Intdiag)

  end subroutine IPD_physics_step2
#endif

!----------------
!  IPD finalize
!----------------
  subroutine IPD_finalize ()
    !--- finalize the physics suite
    call finalize ()
  end subroutine IPD_finalize

end module IPD_driver
