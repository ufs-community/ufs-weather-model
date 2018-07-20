!>\file GFS_rad_time_vary.f90
!!  Contains code related to GFS physics suite setup (radiation part of time_vary_step)
      module GFS_rad_time_vary
      contains

!>\defgroup GFS_rad_time_vary GFS RRTMG Update 
!!\ingroup RRTMG
!! @{
!! \section arg_table_GFS_rad_time_vary_init Argument Table
!!
      subroutine GFS_rad_time_vary_init
      end subroutine GFS_rad_time_vary_init

!> \section arg_table_GFS_rad_time_vary_run Argument Table
!! | local_name        | standard_name                                          | long_name                                                                     | units    | rank |  type                 |   kind    | intent | optional |
!! |-------------------|--------------------------------------------------------|-------------------------------------------------------------------------------|----------|------|-----------------------|-----------|--------|----------|
!! | Model             | FV3-GFS_Control_type                                   | Fortran DDT containing FV3-GFS model control parameters                       | DDT      |    0 | GFS_control_type      |           | inout  | F        |
!! | Statein           | FV3-GFS_Statein_type                                   | Fortran DDT containing FV3-GFS prognostic state data in from dycore           | DDT      |    0 | GFS_statein_type      |           | in     | F        |
!! | Tbd               | FV3-GFS_Tbd_type                                       | Fortran DDT containing FV3-GFS data not yet assigned to a defined container   | DDT      |    0 | GFS_tbd_type          |           | inout  | F        |
!! | errmsg            | error_message                                          | error message for error handling in CCPP                                      | none     |    0 | character             | len=*     | out    | F        |
!! | errflg            | error_flag                                             | error flag for error handling in CCPP                                         | flag     |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_rad_time_vary_run (Model, Statein, Tbd, errmsg, errflg)

      use physparam,                 only: ipsd0, ipsdlim, iaerflg
      use mersenne_twister,          only: random_setseed, random_index, random_stat
      use machine,                   only: kind_phys
      use GFS_typedefs,              only: GFS_statein_type,   &
                                           GFS_control_type,   &
                                           GFS_grid_type,      &
                                           GFS_tbd_type
      use GFS_radupdate,             only: GFS_radupdate_run
      use radcons,                   only: qmin, con_100

      implicit none

      type(GFS_control_type), intent(inout) :: Model
      type(GFS_statein_type), intent(in)    :: Statein
      type(GFS_tbd_type),     intent(inout) :: Tbd
      character(len=*),       intent(out) :: errmsg
      integer,                intent(out) :: errflg

      !--- local variables
      type (random_stat) :: stat
      integer :: ix, nb, j, i, nblks, ipseed
      integer :: numrdm(Model%cnx*Model%cny*2)

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (Model%lsswr .or. Model%lslwr) then

        if (Tbd%blkno==1) then

          call GFS_radupdate_run (Model%idat, Model%jdat, Model%fhswr, Model%dtf, Model%lsswr, &
                        Model%me, Model%slag, Model%sdec, Model%cdec, Model%solcon,            &
                        Model%ictm, Model%isol )
        endif

        !--- set up random seed index in a reproducible way for entire cubed-sphere face (lat-lon grid)
        if ((Model%isubc_lw==2) .or. (Model%isubc_sw==2)) then
          ipseed = mod(nint(con_100*sqrt(Model%sec)), ipsdlim) + 1 + ipsd0
          call random_setseed (ipseed, stat)
          call random_index (ipsdlim, numrdm, stat)
    
          !--- set the random seeds for each column in a reproducible way
          ix = 0
          nb = 1
          ! DH* TODO - this could be sped up by saving jsc, jec, isc, iec in Tbd (for example)
          ! and looping just over them; ix would then run from 1 to blksz(nb); one could also
          ! use OpenMP to speed up this loop *DH
          do j = 1,Model%ny
            do i = 1,Model%nx
              ix = ix + 1
              if (ix .gt. Model%blksz(nb)) then
                ix = 1
                nb = nb + 1
              endif
              if (nb == Tbd%blkno) then
                !--- for testing purposes, replace numrdm with '100'
                Tbd%icsdsw(ix) = numrdm(i+Model%isc-1 + (j+Model%jsc-2)*Model%cnx)
                Tbd%icsdlw(ix) = numrdm(i+Model%isc-1 + (j+Model%jsc-2)*Model%cnx + Model%cnx*Model%cny)
              endif
            enddo
          enddo
        endif  ! isubc_lw and isubc_sw

        if (Model%num_p3d == 4) then
          if (Model%kdt == 1) then
            Tbd%phy_f3d(:,:,1) = Statein%tgrs
            Tbd%phy_f3d(:,:,2) = max(qmin,Statein%qgrs(:,:,1))
            Tbd%phy_f3d(:,:,3) = Statein%tgrs
            Tbd%phy_f3d(:,:,4) = max(qmin,Statein%qgrs(:,:,1))
            Tbd%phy_f2d(:,1)   = Statein%prsi(:,1)
            Tbd%phy_f2d(:,2)   = Statein%prsi(:,1)
          endif
        endif

      endif

  end subroutine GFS_rad_time_vary_run
 
!> \section arg_table_GFS_rad_time_vary_finalize Argument Table
!!
  subroutine GFS_rad_time_vary_finalize()
  end subroutine GFS_rad_time_vary_finalize
!! @}
  end module GFS_rad_time_vary
