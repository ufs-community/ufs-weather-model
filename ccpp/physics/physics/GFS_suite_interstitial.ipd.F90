!> \file GFS_suite_interstitial.f90
!!  Contains code related to more than one scheme in the GFS physics suite.

      module GFS_suite_interstitial_1

      contains

      subroutine GFS_suite_interstitial_1_init ()
      end subroutine GFS_suite_interstitial_1_init

      subroutine GFS_suite_interstitial_1_finalize()
      end subroutine GFS_suite_interstitial_1_finalize

      subroutine GFS_suite_interstitial_1_run (Model, Grid, tottracer, trc_shft, tracers, ntk, skip_macro, &
                                               clw, cnvc, cnvw, errmsg, errflg)

        use machine,               only: kind_phys
        use GFS_typedefs,          only: GFS_control_type, GFS_grid_type

        implicit none

        type(GFS_control_type),           intent(in) :: Model
        type(GFS_grid_type),              intent(in) :: Grid
        integer,                          intent(out) :: tottracer, trc_shft, tracers, ntk
        logical, dimension(size(Grid%xlon,1)), intent(out) :: skip_macro
        real(kind=kind_phys), allocatable, intent(out) :: clw(:,:,:), cnvc(:,:), cnvw(:,:)

        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        tottracer = 0            ! no convective transport of tracers
        if (Model%trans_trac .or. Model%cscnv) then
          if (Model%ntcw > 0) then
            if (Model%ntoz < Model%ntcw) then
              trc_shft = Model%ntcw + Model%ncld - 1
            else
              trc_shft = Model%ntoz
            endif
          elseif (Model%ntoz > 0) then
            trc_shft = Model%ntoz
          else
            trc_shft = 1
          endif

          tracers   = Model%ntrac - trc_shft
          tottracer = tracers
          if (Model%ntoz > 0) tottracer = tottracer + 1  ! ozone is added separately
        endif
        if (Model%ntke > 0) ntk = Model%ntke - trc_shft + 3

        skip_macro = .false.

        allocate ( clw(size(Grid%xlon,1),Model%levs,tottracer+2) )
        if (Model%imfdeepcnv >= 0 .or. Model%imfshalcnv > 0) then
          allocate (cnvc(size(Grid%xlon,1),Model%levs), cnvw(size(Grid%xlon,1),Model%levs))
        endif

      end subroutine GFS_suite_interstitial_1_run

    end module GFS_suite_interstitial_1

    module GFS_suite_interstitial_2

    contains

    subroutine GFS_suite_interstitial_2_init ()
    end subroutine GFS_suite_interstitial_2_init

    subroutine GFS_suite_interstitial_2_finalize()
    end subroutine GFS_suite_interstitial_2_finalize

    subroutine GFS_suite_interstitial_2_run (Model, Grid, Sfcprop, Statein, Diag, rhbbot, rhpbl, rhbtop, frain, islmsk, &
                                             work1, work2, dudt, dvdt, dtdt, dtdtc, dqdt, errmsg, errflg)

      use machine,               only: kind_phys
      use physcons,              only: dxmin, dxinv
      use GFS_typedefs,          only: GFS_control_type, GFS_grid_type, GFS_sfcprop_type, GFS_statein_type, GFS_diag_type

      implicit none

      type(GFS_control_type),           intent(in) :: Model
      type(GFS_grid_type),              intent(in) :: Grid
      type(GFS_sfcprop_type),           intent(in) :: Sfcprop
      type(GFS_statein_type),           intent(in) :: Statein
      type(GFS_diag_type),              intent(inout) :: Diag

      real(kind=kind_phys), intent(out) :: rhbbot, rhpbl, rhbtop, frain
      integer, dimension(size(Grid%xlon,1)), intent(out) :: islmsk
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(out)  :: work1, work2
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(out) :: dudt, dvdt, dtdt, dtdtc
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs,Model%ntrac), intent(out) ::  dqdt

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: i

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      rhbbot = Model%crtrh(1)
      rhpbl  = Model%crtrh(2)
      rhbtop = Model%crtrh(3)

      frain = Model%dtf / Model%dtp

      do i = 1, size(Grid%xlon,1)
        islmsk(i)   = nint(Sfcprop%slmsk(i))
        work1(i)   = (log(Grid%area(i)) - dxmin) * dxinv
        work1(i)   = max(0.0, min(1.0,work1(i)))
        work2(i)   = 1.0 - work1(i)
        Diag%psurf(i)   = Statein%pgr(i)
      end do

      dudt(:,:)  = 0.
      dvdt(:,:)  = 0.
      dtdt(:,:)  = 0.
      dtdtc(:,:) = 0.
      dqdt(:,:,:) = 0.

    end subroutine GFS_suite_interstitial_2_run

  end module GFS_suite_interstitial_2

  module GFS_suite_interstitial_3

  contains

  subroutine GFS_suite_interstitial_3_init ()
  end subroutine GFS_suite_interstitial_3_init

  subroutine GFS_suite_interstitial_3_finalize()
  end subroutine GFS_suite_interstitial_3_finalize

  subroutine GFS_suite_interstitial_3_run (Model, Grid, Statein, Radtend, xcosz, adjsfcdsw, adjsfcdlw, adjsfculw, xmu, &
                                           Diag, kcnv, hflx, evap, errmsg, errflg)

    use machine,               only: kind_phys
    use GFS_typedefs,          only: GFS_control_type, GFS_grid_type, GFS_statein_type, GFS_radtend_type, GFS_diag_type

    implicit none

    type(GFS_control_type),           intent(in)    :: Model
    type(GFS_grid_type),              intent(in)    :: Grid
    type(GFS_statein_type),           intent(in)    :: Statein
    type(GFS_radtend_type),           intent(in)    :: Radtend
    type(GFS_diag_type),              intent(inout) :: Diag

    integer, dimension(size(Grid%xlon,1)), intent(out) :: kcnv
    real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in) :: xcosz, adjsfcdsw, adjsfcdlw, adjsfculw, xmu
    real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(out) :: hflx, evap

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(kind=kind_phys), parameter :: czmin   = 0.0001      ! cos(89.994)

    integer :: i, k

    real(kind=kind_phys) :: tem1

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (Model%lssav) then      !  --- ...  accumulate/save output variables

  !  --- ...  sunshine duration time is defined as the length of time (in mdl output
  !           interval) that solar radiation falling on a plane perpendicular to the
  !           direction of the sun >= 120 w/m2

      do i = 1, size(Grid%xlon,1)
        if ( xcosz(i) >= czmin ) then   ! zenth angle > 89.994 deg
          tem1 = adjsfcdsw(i) / xcosz(i)
          if ( tem1 >= 120.0 ) then
            Diag%suntim(i) = Diag%suntim(i) + Model%dtf
          endif
        endif
      enddo

  !  --- ...  sfc lw fluxes used by atmospheric model are saved for output

      Diag%dlwsfc(:) = Diag%dlwsfc(:) +   adjsfcdlw(:)*Model%dtf
      Diag%ulwsfc(:) = Diag%ulwsfc(:) +   adjsfculw(:)*Model%dtf
      Diag%psmean(:) = Diag%psmean(:) + Statein%pgr(:)*Model%dtf        ! mean surface pressure

      if (Model%ldiag3d) then
        if (Model%lsidea) then
          Diag%dt3dt(:,:,1) = Diag%dt3dt(:,:,1) + Radtend%lwhd(:,:,1)*Model%dtf
          Diag%dt3dt(:,:,2) = Diag%dt3dt(:,:,2) + Radtend%lwhd(:,:,2)*Model%dtf
          Diag%dt3dt(:,:,3) = Diag%dt3dt(:,:,3) + Radtend%lwhd(:,:,3)*Model%dtf
          Diag%dt3dt(:,:,4) = Diag%dt3dt(:,:,4) + Radtend%lwhd(:,:,4)*Model%dtf
          Diag%dt3dt(:,:,5) = Diag%dt3dt(:,:,5) + Radtend%lwhd(:,:,5)*Model%dtf
          Diag%dt3dt(:,:,6) = Diag%dt3dt(:,:,6) + Radtend%lwhd(:,:,6)*Model%dtf
        else
          do k = 1, Model%levs
            Diag%dt3dt(:,k,1) = Diag%dt3dt(:,k,1) + Radtend%htrlw(:,k)*Model%dtf
            Diag%dt3dt(:,k,2) = Diag%dt3dt(:,k,2) + Radtend%htrsw(:,k)*Model%dtf*xmu(:)
          enddo
        endif
      endif
    endif    ! end if_lssav_block

    kcnv(:)   = 0

    hflx(:)       = 0.0
    evap(:)       = 0.0

    Diag%t1(:)      = Statein%tgrs(:,1)
    Diag%q1(:)      = Statein%qgrs(:,1,1)
    Diag%u1(:)      = Statein%ugrs(:,1)
    Diag%v1(:)      = Statein%vgrs(:,1)

  end subroutine GFS_suite_interstitial_3_run

  end module GFS_suite_interstitial_3

  module GFS_suite_update_stateout

  contains

  subroutine GFS_suite_update_stateout_init ()
  end subroutine GFS_suite_update_stateout_init

  subroutine GFS_suite_update_stateout_finalize()
  end subroutine GFS_suite_update_stateout_finalize

  subroutine GFS_suite_update_stateout_run (Statein, Model, Grid, dudt, dvdt, dtdt, dqdt, Stateout, errmsg, errflg)

    use machine,               only: kind_phys
    use GFS_typedefs,          only: GFS_control_type, GFS_statein_type, GFS_grid_type, GFS_stateout_type

    implicit none

    type(GFS_control_type),           intent(in)    :: Model
    type(GFS_statein_type),           intent(in)    :: Statein
    type(GFS_grid_type),              intent(in)    :: Grid
    type(GFS_stateout_type),          intent(inout) :: Stateout

    real(kind=kind_phys), dimension(size(Grid%xlon,1), Model%levs), intent(in) :: dudt, dvdt, dtdt
    real(kind=kind_phys), dimension(size(Grid%xlon,1), Model%levs, Model%ntrac), intent(in) :: dqdt

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    Stateout%gt0(:,:)   = Statein%tgrs(:,:) + dtdt(:,:) * Model%dtp
    Stateout%gu0(:,:)   = Statein%ugrs(:,:) + dudt(:,:) * Model%dtp
    Stateout%gv0(:,:)   = Statein%vgrs(:,:) + dvdt(:,:) * Model%dtp
    Stateout%gq0(:,:,:) = Statein%qgrs(:,:,:) + dqdt(:,:,:) * Model%dtp

  end subroutine GFS_suite_update_stateout_run

end module GFS_suite_update_stateout

module GFS_suite_interstitial_4

contains

subroutine GFS_suite_interstitial_4_init ()
end subroutine GFS_suite_interstitial_4_init

subroutine GFS_suite_interstitial_4_finalize()
end subroutine GFS_suite_interstitial_4_finalize

subroutine GFS_suite_interstitial_4_run (Model, Grid, Statein, rhbbot, rhbtop, work1, work2, clw, cnvc, cnvw, &
                                         ktop, kbot, rhc, errmsg, errflg)

  use machine,               only: kind_phys
  use GFS_typedefs,          only: GFS_control_type, GFS_grid_type, GFS_statein_type
  use physcons,              only: rhc_max

  implicit none

  type(GFS_control_type),           intent(in)    :: Model
  type(GFS_grid_type),              intent(in)    :: Grid
  type(GFS_statein_type),           intent(in)    :: Statein

  real(kind=kind_phys), intent(in)                                           :: rhbbot, rhbtop
  real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in)             :: work1, work2
  real(kind=kind_phys), intent(inout)                                        :: clw(:,:,:), cnvc(:,:), cnvw(:,:)
  integer, dimension(size(Grid%xlon,1)), intent(inout)                       :: ktop, kbot
  real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(out) :: rhc

  character(len=*), intent(out) :: errmsg
  integer,          intent(out) :: errflg

  integer :: i,k
  real(kind=kind_phys) :: tem

  ! Initialize CCPP error handling variables
  errmsg = ''
  errflg = 0

  clw(:,:,1) = 0.0
  clw(:,:,2) = -999.9
  if ((Model%imfdeepcnv >= 0) .or. (Model%imfshalcnv > 0)) then
    cnvc(:,:)  = 0.0
    cnvw(:,:)  = 0.0
  endif

  ktop(:)  = 1
  kbot(:)  = Model%levs
  rhc(:,:) = 0.0

  if (Model%ntcw > 0) then
    do k=1,Model%levs
      do i=1, size(Grid%xlon,1)
        tem      = rhbbot - (rhbbot-rhbtop) * (1.0-Statein%prslk(i,k))
        tem      = rhc_max * work1(i) + tem * work2(i)
        rhc(i,k) = max(0.0, min(1.0,tem))
      enddo
    enddo
  endif

end subroutine GFS_suite_interstitial_4_run

end module GFS_suite_interstitial_4

module GFS_suite_interstitial_5

contains

subroutine GFS_suite_interstitial_5_init ()
end subroutine GFS_suite_interstitial_5_init

subroutine GFS_suite_interstitial_5_finalize()
end subroutine GFS_suite_interstitial_5_finalize

subroutine GFS_suite_interstitial_5_run (clw, cnvc, cnvw, errmsg, errflg)

  use machine,               only: kind_phys

  implicit none

  real(kind=kind_phys), allocatable, intent(inout) :: clw(:,:,:), cnvc(:,:), cnvw(:,:)

  character(len=*), intent(out) :: errmsg
  integer,          intent(out) :: errflg

  ! Initialize CCPP error handling variables
  errmsg = ''
  errflg = 0

  deallocate (clw)
  if (allocated(cnvc)) deallocate(cnvc)
  if (allocated(cnvw)) deallocate(cnvw)

end subroutine GFS_suite_interstitial_5_run

end module GFS_suite_interstitial_5
