module noah_driver


  use noah_loop, only: noah_loop_init, noah_loop_run
  
  implicit none

  public :: noah_loop_drv
  
contains

  subroutine init_driver()

    implicit none
    !private

    integer                 ::   im         ! horiz dimension
    integer                 ::   km         ! vertical soil layer

    integer                 ::   isot       ! sfc soil type data source 
    integer                 ::   ivegsrc    ! sfc veg type data source
    
    character(len=128)      ::   errmsg     ! error messaging added to ccpp
    integer                 ::   errflg     ! error messaging added to ccpp


    ! quick test, should be read in from nml
    isot = 1
    ivegsrc = 1
    
    call noah_loop_init(0, isot, ivegsrc, 0 , errmsg, errflg)
    
  end subroutine init_driver

  subroutine noah_loop_drv()

  end subroutine noah_loop_drv

end module noah_driver
