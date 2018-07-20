!> \file GFS_phys_time_vary.f90
!!  Contains code related to GFS physics suite setup (physics part of time_vary_step)

    module GFS_phys_time_vary_1

      contains

      subroutine GFS_phys_time_vary_1_init ()
      end subroutine GFS_phys_time_vary_1_init

      subroutine GFS_phys_time_vary_1_finalize()
      end subroutine GFS_phys_time_vary_1_finalize

!> \section arg_table_GFS_phys_time_vary_1_run Argument Table
!! | local_name     | standard_name                                          | long_name                                                               | units    | rank |  type                 |   kind    | intent | optional |
!! |----------------|--------------------------------------------------------|-------------------------------------------------------------------------|----------|------|-----------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                   | Fortran DDT containing FV3-GFS model control parameters                 | DDT      |    0 | GFS_control_type      |           | inout  | F        |
!! | errmsg         | error_message                                          | error message for error handling in CCPP                                | none     |    0 | character             | len=*     | out    | F        |
!! | errflg         | error_flag                                             | error flag for error handling in CCPP                                   | flag     |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_phys_time_vary_1_run (Model, errmsg, errflg)

        use machine,               only: kind_phys
        use GFS_typedefs,          only: GFS_control_type

        implicit none

        type(GFS_control_type),           intent(inout) :: Model
        character(len=*),                 intent(out)   :: errmsg
        integer,                          intent(out)   :: errflg

        real(kind=kind_phys), parameter :: con_24  =   24.0_kind_phys
        real(kind=kind_phys), parameter :: con_hr  = 3600.0_kind_phys
        real(kind=kind_phys) :: rinc(5)

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        !--- Model%jdat is being updated directly inside of FV3GFS_cap.F90
        !--- update calendars and triggers
        rinc(1:5)   = 0
        call w3difdat(Model%jdat,Model%idat,4,rinc)
        Model%sec = rinc(4)

        Model%phour = Model%sec/con_hr
        !--- set current bucket hour
        Model%zhour = Model%phour
        Model%fhour = (Model%sec + Model%dtp)/con_hr
        Model%kdt   = nint((Model%sec + Model%dtp)/Model%dtp)

        Model%ipt    = 1
        Model%lprnt  = .false.
        Model%lssav  = .true.

        !--- radiation triggers
        Model%lsswr  = (mod(Model%kdt, Model%nsswr) == 1)
        Model%lslwr  = (mod(Model%kdt, Model%nslwr) == 1)

        !--- set the solar hour based on a combination of phour and time initial hour
        Model%solhr  = mod(Model%phour+Model%idate(1),con_24)

        if ((Model%debug) .and. (Model%me == Model%master)) then
          print *,'   sec ', Model%sec
          print *,'   kdt ', Model%kdt
          print *,' nsswr ', Model%nsswr
          print *,' nslwr ', Model%nslwr
          print *,' nscyc ', Model%nscyc
          print *,' lsswr ', Model%lsswr
          print *,' lslwr ', Model%lslwr
          print *,' fhour ', Model%fhour
          print *,' phour ', Model%phour
          print *,' solhr ', Model%solhr
        endif

      end subroutine GFS_phys_time_vary_1_run

    end module GFS_phys_time_vary_1

    module GFS_phys_time_vary_2

      contains

      subroutine GFS_phys_time_vary_2_init ()
      end subroutine GFS_phys_time_vary_2_init

      subroutine GFS_phys_time_vary_2_finalize()
      end subroutine GFS_phys_time_vary_2_finalize

!> \section arg_table_GFS_phys_time_vary_2_run Argument Table
!! | local_name     | standard_name                                          | long_name                                                               | units    | rank |  type                 |   kind    | intent | optional |
!! |----------------|--------------------------------------------------------|-------------------------------------------------------------------------|----------|------|-----------------------|-----------|--------|----------|
!! | Grid           | FV3-GFS_Grid_type_all_blocks                           | Fortran DDT containing FV3-GFS grid and interpolation related data      | DDT      |    1 | GFS_grid_type         |           | in     | F        |
!! | Model          | FV3-GFS_Control_type                                   | Fortran DDT containing FV3-GFS model control parameters                 | DDT      |    0 | GFS_control_type      |           | inout  | F        |
!! | Tbd            | FV3-GFS_Tbd_type_all_blocks                            | Fortran DDT containing FV3-GFS miscellaneous data                       | DDT      |    1 | GFS_tbd_type          |           | inout  | F        |
!! | Sfcprop        | FV3-GFS_Sfcprop_type_all_blocks                        | Fortran DDT containing FV3-GFS surface fields                           | DDT      |    1 | GFS_sfcprop_type      |           | inout  | F        |
!! | Cldprop        | FV3-GFS_Cldprop_type_all_blocks                        | Fortran DDT containing FV3-GFS cloud fields                             | DDT      |    1 | GFS_cldprop_type      |           | inout  | F        |
!! | Diag           | FV3-GFS_Diag_type_all_blocks                           | Fortran DDT containing FV3-GFS fields targeted for diagnostic output    | DDT      |    1 | GFS_diag_type         |           | inout  | F        |
!! | errmsg         | error_message                                          | error message for error handling in CCPP                                | none     |    0 | character             | len=*     | out    | F        |
!! | errflg         | error_flag                                             | error flag for error handling in CCPP                                   | flag     |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_phys_time_vary_2_run (Grid, Model, Tbd, Sfcprop, Cldprop, Diag, errmsg, errflg)

        use mersenne_twister, only: random_setseed, random_number
        use machine,               only: kind_phys
        use physcons,              only: dxmin, dxinv
        use GFS_typedefs,          only: GFS_control_type, GFS_grid_type, &
                                         GFS_tbd_type, GFS_sfcprop_type,  &
                                         GFS_cldprop_type, GFS_diag_type

        implicit none

        type(GFS_grid_type),              intent(in)    :: Grid(:)
        type(GFS_control_type),           intent(inout) :: Model
        type(GFS_tbd_type),               intent(inout) :: Tbd(:)
        type(GFS_sfcprop_type),           intent(inout) :: Sfcprop(:)
        type(GFS_cldprop_type),           intent(inout) :: Cldprop(:)
        type(GFS_diag_type),              intent(inout) :: Diag(:)
        character(len=*),                 intent(out)   :: errmsg
        integer,                          intent(out)   :: errflg

        real(kind=kind_phys), parameter :: con_hr  = 3600.0_kind_phys
        real(kind=kind_phys), parameter :: con_99  =   99.0_kind_phys
        real(kind=kind_phys), parameter :: con_100 =  100.0_kind_phys

        integer :: i, j, k, iseed, iskip, ix, nb, nblks
        real(kind=kind_phys) :: wrk(1)
        real(kind=kind_phys) :: rannie(Model%cny)
        real(kind=kind_phys) :: rndval(Model%cnx*Model%cny*Model%nrcm)

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        nblks = size(Model%blksz)

        !--- switch for saving convective clouds - cnvc90.f
        !--- aka Ken Campana/Yu-Tai Hou legacy
        if ((mod(Model%kdt,Model%nsswr) == 0) .and. (Model%lsswr)) then
          !--- initialize,accumulate,convert
          Model%clstp = 1100 + min(Model%fhswr/con_hr,Model%fhour,con_99)
        elseif (mod(Model%kdt,Model%nsswr) == 0) then
          !--- accumulate,convert
          Model%clstp = 0100 + min(Model%fhswr/con_hr,Model%fhour,con_99)
        elseif (Model%lsswr) then
          !--- initialize,accumulate
          Model%clstp = 1100
        else
          !--- accumulate
          Model%clstp = 0100
        endif

        !--- random number needed for RAS and old SAS and when cal_pre=.true.
        if ( ((Model%imfdeepcnv <= 0) .or. (Model%cal_pre)) .and. (Model%random_clds) ) then
          iseed = mod(con_100*sqrt(Model%fhour*con_hr),1.0d9) + Model%seed0
          call random_setseed(iseed)
          call random_number(wrk)
          do i = 1,Model%cnx*Model%nrcm
            iseed = iseed + nint(wrk(1)) * i
            call random_setseed(iseed)
            call random_number(rannie)
            rndval(1+(i-1)*Model%cny:i*Model%cny) = rannie(1:Model%cny)
          enddo

          ! DH* TODO - this could be sped up by saving jsc, jec, isc, iec in Tbd (for example)
          ! and looping just over them; ix would then run from 1 to blksz(nb); one could also
          ! use OpenMP to speed up this loop or the inside loops *DH
          do k = 1,Model%nrcm
            iskip = (k-1)*Model%cnx*Model%cny
            ix = 0
            nb = 1
            do j = 1,Model%ny
              do i = 1,Model%nx
                ix = ix + 1
                if (ix .gt. Model%blksz(nb)) then
                  ix = 1
                  nb = nb + 1
                endif
                Tbd(nb)%rann(ix,k) = rndval(i+Model%isc-1 + (j+Model%jsc-2)*Model%cnx + iskip)
              enddo
            enddo
          enddo
        endif  ! imfdeepcnv, cal_re, random_clds

        !--- o3 interpolation
        if (Model%ntoz > 0) then
! DH* OpenMP?
          do nb = 1, nblks
            call ozinterpol (Model%me, Model%blksz(nb), Model%idate, Model%fhour, &
                             Grid(nb)%jindx1_o3, Grid(nb)%jindx2_o3,              &
                             Tbd(nb)%ozpl, Grid(nb)%ddy_o3)
          enddo
        endif

        !--- h2o interpolation
        if (Model%h2o_phys) then
! DH* OpenMP?
          do nb = 1, nblks
            call h2ointerpol (Model%me, Model%blksz(nb), Model%idate, Model%fhour, &
                              Grid(nb)%jindx1_h, Grid(nb)%jindx2_h,                &
                              Tbd(nb)%h2opl, Grid(nb)%ddy_h)
          enddo
        endif

        !--- repopulate specific time-varying sfc properties for AMIP/forecast runs
        if (Model%nscyc >  0) then
          if (mod(Model%kdt,Model%nscyc) == 1) THEN
            call gcycle (nblks, Model, Grid(:), Sfcprop(:), Cldprop(:))
          endif
        endif

        !--- determine if diagnostics buckets need to be cleared
        if (mod(Model%kdt,Model%nszero) == 1) then
! DH* OpenMP?
          do nb = 1,nblks
            call Diag(nb)%rad_zero  (Model)
            call Diag(nb)%phys_zero (Model)
        !!!!  THIS IS THE POINT AT WHICH DIAG%ZHOUR NEEDS TO BE UPDATED
          enddo
        endif

      end subroutine GFS_phys_time_vary_2_run

    end module GFS_phys_time_vary_2
