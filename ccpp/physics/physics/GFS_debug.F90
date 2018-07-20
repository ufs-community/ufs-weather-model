!> \file GFS_debug.F90

    module GFS_diagtoscreen

      private
 
      public GFS_diagtoscreen_init, GFS_diagtoscreen_run, GFS_diagtoscreen_finalize

      interface print_var
        module procedure print_int_0d
        module procedure print_real_0d
        module procedure print_real_1d
        module procedure print_real_2d
        module procedure print_real_3d
      end interface

      integer, parameter :: ISTART = 1
      integer, parameter :: IEND = 11

      integer, parameter :: KSTART = 1
      integer, parameter :: KEND = 11

      contains

      subroutine GFS_diagtoscreen_init ()
      end subroutine GFS_diagtoscreen_init

      subroutine GFS_diagtoscreen_finalize ()
      end subroutine GFS_diagtoscreen_finalize

!> \section arg_table_GFS_diagtoscreen_run Argument Table
!! | local_name     | standard_name                                          | long_name                                               | units         | rank | type                  |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                   | derived type GFS_control_type in FV3                    | DDT           |    0 | GFS_control_type      |           | in     | F        |
!! | Statein        | FV3-GFS_Statein_type                                   | derived type GFS_statein_type in FV3                    | DDT           |    0 | GFS_statein_type      |           | in     | F        |
!! | Stateout       | FV3-GFS_Stateout_type                                  | derived type GFS_stateout_type in FV3                   | DDT           |    0 | GFS_stateout_type     |           | in     | F        |
!! | Sfcprop        | FV3-GFS_Sfcprop_type                                   | derived type GFS_sfcprop_type in FV3                    | DDT           |    0 | GFS_sfcprop_type      |           | in     | F        |
!! | Coupling       | FV3-GFS_Coupling_type                                  | derived type GFS_coupling_type in FV3                   | DDT           |    0 | GFS_coupling_type     |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                      | derived type GFS_grid_type in FV3                       | DDT           |    0 | GFS_grid_type         |           | in     | F        |
!! | Tbd            | FV3-GFS_Tbd_type                                       | derived type GFS_tbd_type in FV3                        | DDT           |    0 | GFS_tbd_type          |           | in     | F        |
!! | Cldprop        | FV3-GFS_Cldprop_type                                   | derived type GFS_cldprop_type in FV3                    | DDT           |    0 | GFS_cldprop_type      |           | in     | F        |
!! | Radtend        | FV3-GFS_Radtend_type                                   | derived type GFS_radtend_type in FV3                    | DDT           |    0 | GFS_radtend_type      |           | in     | F        |
!! | Diag           | FV3-GFS_Diag_type                                      | derived type GFS_diag_type in FV3                       | DDT           |    0 | GFS_diag_type         |           | in     | F        |
!! | Interstitial   | FV3-GFS_Interstitial_type                              | derived type GFS_interstitial_type in FV3               | DDT           |    0 | GFS_interstitial_type |           | in     | F        |
!! | errmsg         | error_message                                          | error message for error handling in CCPP                | none          |    0 | character             | len=*     | out    | F        |
!! | errflg         | error_flag                                             | error flag for error handling in CCPP                   | flag          |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_diagtoscreen_run (Model, Statein, Stateout, Sfcprop, Coupling, &
                                       Grid, Tbd, Cldprop, Radtend, Diag, Interstitial, &
                                       errmsg, errflg)

#ifdef MPI
         use mpi
#endif
#ifdef OPENMP
         use omp_lib
#endif
         use machine,               only: kind_phys
         use GFS_typedefs,          only: GFS_control_type, GFS_statein_type,  &
                                          GFS_stateout_type, GFS_sfcprop_type, &
                                          GFS_coupling_type, GFS_grid_type,    &
                                          GFS_tbd_type, GFS_cldprop_type,      &
                                          GFS_radtend_type, GFS_diag_type,     &
                                          GFS_interstitial_type

         implicit none

         !--- interface variables
         type(GFS_control_type),   intent(in   ) :: Model
         type(GFS_statein_type),   intent(in   ) :: Statein
         type(GFS_stateout_type),  intent(in   ) :: Stateout
         type(GFS_sfcprop_type),   intent(in   ) :: Sfcprop
         type(GFS_coupling_type),  intent(in   ) :: Coupling
         type(GFS_grid_type),      intent(in   ) :: Grid
         type(GFS_tbd_type),       intent(in   ) :: Tbd
         type(GFS_cldprop_type),   intent(in   ) :: Cldprop
         type(GFS_radtend_type),   intent(in   ) :: Radtend
         type(GFS_diag_type),      intent(in   ) :: Diag
         type(GFS_interstitial_type), intent(in) :: Interstitial
         character(len=*),           intent(out) :: errmsg
         integer,                    intent(out) :: errflg

         !--- local variables
         integer :: impi, iomp, ierr
         integer :: mpirank,mpisize
         integer :: omprank,ompsize

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

#ifdef MPI
         call MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, mpisize, ierr)
#else
         mpirank = 0
         mpisize = 1
#endif
#ifdef OPENMP
         omprank = OMP_GET_THREAD_NUM()
         ompsize = OMP_GET_NUM_THREADS()
#else
         omprank = 0
         ompsize = 1
#endif

#ifdef OPENMP
!$OMP BARRIER
#endif
#ifdef MPI
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

         do impi=0,mpisize-1
             do iomp=0,ompsize-1
                 if (mpirank==impi .and. omprank==iomp) then
                     call print_var(mpirank,omprank,Tbd%blkno, 'Sfcprop%tsfc',         Sfcprop%tsfc)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Statein%tgrs',         Statein%tgrs)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Statein%ugrs',         Statein%ugrs)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Statein%vgrs',         Statein%vgrs)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Statein%qgrs',         Statein%qgrs)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Statein%prsl',         Statein%prsl)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Statein%pgr',          Statein%pgr)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Diag%topfsw%upfxc',    Diag%topfsw%upfxc)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Diag%topfsw%dnfxc',    Diag%topfsw%dnfxc)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Diag%topfsw%upfx0',    Diag%topfsw%upfx0)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Diag%topflw%upfxc',    Diag%topflw%upfxc)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Diag%topflw%upfx0',    Diag%topflw%upfx0)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Sfcprop%tg3',          Sfcprop%tg3)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Sfcprop%smc',          Sfcprop%smc)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Sfcprop%stc',          Sfcprop%stc)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Sfcprop%t2m',          Sfcprop%t2m)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Sfcprop%q2m',          Sfcprop%q2m)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Sfcprop%tref',         Sfcprop%tref)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Tbd%htlwc',            Tbd%htlwc)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Tbd%htlw0',            Tbd%htlw0)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Tbd%htswc',            Tbd%htswc)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Tbd%htsw0',            Tbd%htsw0)
                     call print_var(mpirank,omprank,Tbd%blkno, 'Model%sec',            Model%sec)
                 end if
#ifdef OPENMP
!$OMP BARRIER
#endif
             end do
#ifdef MPI
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
         end do

#ifdef OPENMP
!$OMP BARRIER
#endif
#ifdef MPI
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

      end subroutine GFS_diagtoscreen_run

      subroutine print_int_0d(mpirank,omprank,blkno,name,var)

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          integer, intent(in) :: var

          write(0,'(2a,3i4,i10)') 'XXX: ', trim(name), mpirank, omprank, blkno, var

      end subroutine print_int_0d

      subroutine print_real_0d(mpirank,omprank,blkno,name,var)

          use machine,               only: kind_phys

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          real(kind_phys), intent(in) :: var

          write(0,'(2a,3i4,e14.6)') 'XXX: ', trim(name), mpirank, omprank, blkno, var

      end subroutine print_real_0d

      subroutine print_real_1d(mpirank,omprank,blkno,name,var)

          use machine,               only: kind_phys

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          real(kind_phys), intent(in) :: var(:)

          integer :: i

          do i=ISTART,min(IEND,size(var(:)))
              write(0,'(2a,3i4,i6,e14.6)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, var(i)
          end do

      end subroutine print_real_1d

      subroutine print_real_2d(mpirank,omprank,blkno,name,var)

          use machine,               only: kind_phys

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          real(kind_phys), intent(in) :: var(:,:)
          
          integer :: k, i

          do i=ISTART,min(IEND,size(var(:,1)))
              do k=KSTART,min(KEND,size(var(1,:)))
                  write(0,'(2a,3i4,2i6,e14.6)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, k, var(i,k)
              end do
          end do

      end subroutine print_real_2d

      subroutine print_real_3d(mpirank,omprank,blkno,name,var)

          use machine,               only: kind_phys

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          real(kind_phys), intent(in) :: var(:,:,:)

          integer :: k, i, l

          do i=ISTART,min(IEND,size(var(:,1,1)))
              do k=KSTART,min(KEND,size(var(1,:,1)))
                  do l=1,size(var(1,1,:))
                      write(0,'(2a,3i4,3i6,e14.6)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, k, l, var(i,k,l)
                  end do
              end do
          end do

      end subroutine print_real_3d

    end module GFS_diagtoscreen


    module GFS_interstitialtoscreen

      private
 
      public GFS_interstitialtoscreen_init, GFS_interstitialtoscreen_run, GFS_interstitialtoscreen_finalize

      contains

      subroutine GFS_interstitialtoscreen_init ()
      end subroutine GFS_interstitialtoscreen_init

      subroutine GFS_interstitialtoscreen_finalize ()
      end subroutine GFS_interstitialtoscreen_finalize

!> \section arg_table_GFS_interstitialtoscreen_run Argument Table
!! | local_name     | standard_name                                          | long_name                                               | units         | rank | type                  |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                   | derived type GFS_control_type in FV3                    | DDT           |    0 | GFS_control_type      |           | in     | F        |
!! | Statein        | FV3-GFS_Statein_type                                   | derived type GFS_statein_type in FV3                    | DDT           |    0 | GFS_statein_type      |           | in     | F        |
!! | Stateout       | FV3-GFS_Stateout_type                                  | derived type GFS_stateout_type in FV3                   | DDT           |    0 | GFS_stateout_type     |           | in     | F        |
!! | Sfcprop        | FV3-GFS_Sfcprop_type                                   | derived type GFS_sfcprop_type in FV3                    | DDT           |    0 | GFS_sfcprop_type      |           | in     | F        |
!! | Coupling       | FV3-GFS_Coupling_type                                  | derived type GFS_coupling_type in FV3                   | DDT           |    0 | GFS_coupling_type     |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                      | derived type GFS_grid_type in FV3                       | DDT           |    0 | GFS_grid_type         |           | in     | F        |
!! | Tbd            | FV3-GFS_Tbd_type                                       | derived type GFS_tbd_type in FV3                        | DDT           |    0 | GFS_tbd_type          |           | in     | F        |
!! | Cldprop        | FV3-GFS_Cldprop_type                                   | derived type GFS_cldprop_type in FV3                    | DDT           |    0 | GFS_cldprop_type      |           | in     | F        |
!! | Radtend        | FV3-GFS_Radtend_type                                   | derived type GFS_radtend_type in FV3                    | DDT           |    0 | GFS_radtend_type      |           | in     | F        |
!! | Diag           | FV3-GFS_Diag_type                                      | derived type GFS_diag_type in FV3                       | DDT           |    0 | GFS_diag_type         |           | in     | F        |
!! | Interstitial   | FV3-GFS_Interstitial_type                              | derived type GFS_interstitial_type in FV3               | DDT           |    0 | GFS_interstitial_type |           | in     | F        |
!! | errmsg         | error_message                                          | error message for error handling in CCPP                | none          |    0 | character             | len=*     | out    | F        |
!! | errflg         | error_flag                                             | error flag for error handling in CCPP                   | flag          |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_interstitialtoscreen_run (Model, Statein, Stateout, Sfcprop, Coupling, &
                                           Grid, Tbd, Cldprop, Radtend, Diag, Interstitial, &
                                           errmsg, errflg)

#ifdef MPI
         use mpi
#endif
#ifdef OPENMP
         use omp_lib
#endif
         use machine,               only: kind_phys
         use GFS_typedefs,          only: GFS_control_type, GFS_statein_type,  &
                                          GFS_stateout_type, GFS_sfcprop_type, &
                                          GFS_coupling_type, GFS_grid_type,    &
                                          GFS_tbd_type, GFS_cldprop_type,      &
                                          GFS_radtend_type, GFS_diag_type,     &
                                          GFS_interstitial_type

         implicit none

         !--- interface variables
         type(GFS_control_type),   intent(in   ) :: Model
         type(GFS_statein_type),   intent(in   ) :: Statein
         type(GFS_stateout_type),  intent(in   ) :: Stateout
         type(GFS_sfcprop_type),   intent(in   ) :: Sfcprop
         type(GFS_coupling_type),  intent(in   ) :: Coupling
         type(GFS_grid_type),      intent(in   ) :: Grid
         type(GFS_tbd_type),       intent(in   ) :: Tbd
         type(GFS_cldprop_type),   intent(in   ) :: Cldprop
         type(GFS_radtend_type),   intent(in   ) :: Radtend
         type(GFS_diag_type),      intent(in   ) :: Diag
         type(GFS_interstitial_type), intent(in) :: Interstitial
         character(len=*),           intent(out) :: errmsg
         integer,                    intent(out) :: errflg

         !--- local variables
         integer :: impi, iomp, ierr
         integer :: mpirank,mpisize
         integer :: omprank,ompsize

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

#ifdef MPI
         call MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, mpisize, ierr)
#else
         mpirank = 0
         mpisize = 1
#endif
#ifdef OPENMP
         omprank = OMP_GET_THREAD_NUM()
         ompsize = OMP_GET_NUM_THREADS()
#else
         omprank = 0
         ompsize = 1
#endif

#ifdef OPENMP
!$OMP BARRIER
#endif
#ifdef MPI
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

         do impi=0,mpisize-1
             do iomp=0,ompsize-1
                 if (mpirank==impi .and. omprank==iomp) then
                     call Interstitial%mprint(mpirank,omprank,Tbd%blkno)
                 end if
#ifdef OPENMP
!$OMP BARRIER
#endif
             end do
#ifdef MPI
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
         end do

#ifdef OPENMP
!$OMP BARRIER
#endif
#ifdef MPI
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

      end subroutine GFS_interstitialtoscreen_run

    end module GFS_interstitialtoscreen


    module GFS_barrier

      private
 
      public GFS_barrier_init, GFS_barrier_run, GFS_barrier_finalize

  contains

      subroutine GFS_barrier_init ()
      end subroutine GFS_barrier_init

      subroutine GFS_barrier_finalize ()
      end subroutine GFS_barrier_finalize

!> \section arg_table_GFS_barrier_run Argument Table
!! | local_name     | standard_name                                          | long_name                                               | units         | rank | type                  |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                   | derived type GFS_control_type in FV3                    | DDT           |    0 | GFS_control_type      |           | in     | F        |
!! | errmsg         | error_message                                          | error message for error handling in CCPP                | none          |    0 | character             | len=*     | out    | F        |
!! | errflg         | error_flag                                             | error flag for error handling in CCPP                   | flag          |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_barrier_run (Model, errmsg, errflg)

#ifdef MPI
         use mpi
#endif
#ifdef OPENMP
         use omp_lib
#endif
         use machine,               only: kind_phys
         use GFS_typedefs,          only: GFS_control_type

         implicit none

         !--- interface variables
         type(GFS_control_type),   intent(in)  :: Model
         character(len=*),         intent(out) :: errmsg
         integer,                  intent(out) :: errflg

         !--- local variables
         integer :: impi, iomp, ierr
         integer :: mpirank,mpisize
         integer :: omprank,ompsize

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

#ifdef MPI
         call MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, mpisize, ierr)
#else
         mpirank = 0
         mpisize = 1
#endif
#ifdef OPENMP
         omprank = OMP_GET_THREAD_NUM()
         ompsize = OMP_GET_NUM_THREADS()
#else
         omprank = 0
         ompsize = 1
#endif

#ifdef OPENMP
!$OMP BARRIER
#endif
#ifdef MPI
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
         ! Keep this for flushing output to disk
         call sleep(1)
#ifdef MPI
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

      end subroutine GFS_barrier_run

    end module GFS_barrier

    module GFS_abort

      private
 
      public GFS_abort_init, GFS_abort_run, GFS_abort_finalize

  contains

      subroutine GFS_abort_init ()
      end subroutine GFS_abort_init

      subroutine GFS_abort_finalize ()
      end subroutine GFS_abort_finalize

!> \section arg_table_GFS_abort_run Argument Table
!! | local_name     | standard_name                                          | long_name                                               | units         | rank | type                  |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                   | derived type GFS_control_type in FV3                    | DDT           |    0 | GFS_control_type      |           | in     | F        |
!! | errmsg         | error_message                                          | error message for error handling in CCPP                | none          |    0 | character             | len=*     | out    | F        |
!! | errflg         | error_flag                                             | error flag for error handling in CCPP                   | flag          |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_abort_run (Model, errmsg, errflg)

#ifdef MPI
         use mpi
#endif
#ifdef OPENMP
         use omp_lib
#endif
         use machine,               only: kind_phys
         use GFS_typedefs,          only: GFS_control_type

         implicit none

         !--- interface variables
         type(GFS_control_type),   intent(in)  :: Model
         character(len=*),         intent(out) :: errmsg
         integer,                  intent(out) :: errflg

         !--- local variables
         integer :: impi, iomp, ierr
         integer :: mpirank,mpisize
         integer :: omprank,ompsize

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

#ifdef MPI
         call MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, mpisize, ierr)
#else
         mpirank = 0
         mpisize = 1
#endif
#ifdef OPENMP
         omprank = OMP_GET_THREAD_NUM()
         ompsize = OMP_GET_NUM_THREADS()
#else
         omprank = 0
         ompsize = 1
#endif

         errflg = 1
         errmsg = 'Abort requested by user in GFS_abort_run'

#ifdef OPENMP
!$OMP BARRIER
#endif
#ifdef MPI
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

      end subroutine GFS_abort_run

    end module GFS_abort
