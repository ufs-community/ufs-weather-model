module noah_driver

  use machine, only: kind_phys
  
  use noah_loop, only: noah_loop_init, noah_loop_run
  use noah_type_mod
  
  
  
  implicit none
  private
  
  type (noah_type), public        ::   noah_pubinst
  
  public :: noah_loop_drv, init_driver
  
contains

  subroutine init_driver(procbounds)

    use proc_bounds, only : procbounds_type
    
    type(procbounds_type),  intent(in)    :: procbounds

    ! local
    !type (noah_type)        ::   noah
    integer                 ::   im         ! horiz dimension
    integer                 ::   km         ! vertical soil layer

    integer                 ::   isot       ! sfc soil type data source 
    integer                 ::   ivegsrc    ! sfc veg type data source
    
    character(len=128)      ::   errmsg     ! error messaging added to ccpp
    integer                 ::   errflg     ! error messaging added to ccpp

    im = procbounds%im
    
    ! quick test, should be read in from nml
    isot = 1
    ivegsrc = 1

    call noah_pubinst%Create(im)
    
    call noah_loop_init(0, isot, ivegsrc, 0 , errmsg, errflg)
    
    
  end subroutine init_driver

  subroutine noah_loop_drv(procbounds, noah_pubinst)

    use proc_bounds, only : procbounds_type

    implicit none

    type(procbounds_type),  intent(in)    :: procbounds 
    type(noah_type),        intent(inout)    :: noah_pubinst ! land model's variable type
    
    ! local
    integer                 :: i, de, gridbeg, gridend, im
    real(kind_phys)         :: foodata(procbounds%im)
    !
    
    associate(foodata => noah_pubinst%model%foo_atm2lndfield  )

      ! first test
      !write(*,*) 'NLP test: ', noah_pubinst%model%soiltyp ! get all zeros, good
      
      de = procbounds%de
      im = procbounds%im
      gridbeg = procbounds%gridbeg
      gridend = procbounds%gridend

      write(*,*) 'NLP1: ', de, gridbeg, gridend, im, size(noah_pubinst%model%foo_atm2lndfield)

      ! foodata(1:gridend-gridbeg+1) = noah_pubinst%model%foo_atm2lndfield(gridbeg:gridend)
      ! foodata = noah_pubinst%model%foo_atm2lndfield(gridbeg:gridend)
      ! foodata = noah_pubinst%model%foo_atm2lndfield


      ! do i = 1,im
      !    write(*,*) 'NLP2: ', de, gridbeg, gridend, foodata(i)
      ! end do

    end associate

  end subroutine noah_loop_drv

end module noah_driver
