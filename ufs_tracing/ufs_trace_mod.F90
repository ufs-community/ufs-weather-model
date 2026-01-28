module ufs_trace_mod

  use mpi_f08

  implicit none

  private

  logical :: initialized = .false.
  integer :: unit = -1

  public ufs_trace_init
  public ufs_trace
  public ufs_trace_finalize

contains

  subroutine ufs_trace_init()

    integer :: ierr
    integer :: rank
    character(len=32) :: fname

    if (.not.initialized) then
      call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
      write(fname,'("ufs_trace_",I8.8,".trace")') rank

      open(newunit=unit, file=trim(fname))
      initialized = .true.
    end if

  end subroutine ufs_trace_init


  subroutine ufs_trace(component, routine, ph)
    character(len=*), intent(in) :: component, routine, ph

    character(len=*), parameter :: GFMT = '("{""name"":""",A ,""", &
                                                ""ph"":""",A ,""", &
                                                ""ts"":""",I0,""", &
                                               ""pid"":""",A ,""", &
                                               ""tid"":""",A ,"""},")'

    integer(kind=8) :: cur_time

    cur_time = MPI_Wtime() * 1000000

    write(unit, FMT=GFMT) routine, ph, cur_time, "1", component

  end subroutine ufs_trace

  subroutine ufs_trace_finalize()

    integer :: ierr
    character(len=32) :: fname

    if (.not.initialized) then
       return
    end if

    close(unit)

  end subroutine ufs_trace_finalize

end module ufs_trace_mod
