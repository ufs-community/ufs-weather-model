!> \file memcheck.F90
!!  Contains code to check memory usage with/without CCPP.

    module memcheck

#ifndef CCPP
      use, intrinsic :: iso_c_binding,                                   &
                        only: c_int32_t, c_char, c_null_char
#endif

      use machine, only: kind_phys

      implicit none

      private
 
      public memcheck_init, memcheck_run, memcheck_finalize

#ifndef CCPP
      ! In temporary external library libmemcheck.a
      interface
          integer(c_int32_t)                                             &
          function ccpp_memory_usage_c                                   &
                     (mpicomm, str, lstr)                                &
                     bind(c, name='ccpp_memory_usage_c')
              import :: c_char, c_int32_t
              integer(c_int32_t), value, intent(in) :: mpicomm
              character(kind=c_char), dimension(*)  :: str
              integer(c_int32_t), value, intent(in) :: lstr
          end function ccpp_memory_usage_c
      end interface
#endif

      ! Can use larger time frame to track memory leaks
      real(kind_phys), parameter :: SECONDS_ELAPSED_MIN = 3500.0
      real(kind_phys), parameter :: SECONDS_ELAPSED_MAX = 3700.0

      contains

      subroutine memcheck_init ()
      end subroutine memcheck_init

      subroutine memcheck_finalize ()
      end subroutine memcheck_finalize

!> \section arg_table_memcheck_run Argument Table
!! | local_name      | standard_name                                          | long_name                                               | units         | rank | type      |    kind   | intent | optional |
!! |-----------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | seconds_elapsed | seconds_elapsed_since_model_initialization             | seconds elapsed since model initialization              | s             |    0 | real      | kind_phys | in     | F        |
!! | block_number    | block_number                                           | for explicit data blocking: block number of this block  | index         |    0 | integer   |           | in     | F        |
!! | errmsg          | error_message                                          | error message for error handling in CCPP                | none          |    0 | character | len=*     | out    | F        |
!! | errflg          | error_flag                                             | error flag for error handling in CCPP                   | flag          |    0 | integer   |           | out    | F        |
!!
      subroutine memcheck_run (seconds_elapsed, block_number, errmsg, errflg)

#ifdef MPI
         use mpi
#endif
#ifdef OPENMP
         use omp_lib
#endif
#ifdef CCPP
         use ccpp_api, only: ccpp_memory_usage
#endif

         implicit none

         !--- interface variables
         real(kind=kind_phys),       intent(in)  :: seconds_elapsed
         integer,                    intent(in)  :: block_number
         character(len=*),           intent(out) :: errmsg
         integer,                    intent(out) :: errflg

         !--- local variables
         integer :: impi, ierr
         integer :: mpirank, mpisize, mpicomm ! mpicomm will become input argument
         integer :: ompthread
         character(len=1024) :: memory_usage

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (seconds_elapsed < SECONDS_ELAPSED_MIN .or. &
            seconds_elapsed > SECONDS_ELAPSED_MAX) return

         if (block_number>1) return

#ifdef MPI
         call MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, mpisize, ierr)
         mpicomm = MPI_COMM_WORLD
#else
         mpirank = 0
         mpisize = 1
         mpicomm = 0
#endif

#ifdef OPENMP
         ompthread = OMP_GET_THREAD_NUM()
#else
         ompthread = 0
#endif
 
         ierr = ccpp_memory_usage(mpicomm, memory_usage)

         ! Output ordered by MPI rank
         do impi=0,mpisize-1
            if (mpirank==impi) then
                write(0,'(a)') trim(memory_usage)
            end if
#ifdef MPI
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
         end do

#ifndef CCPP
! Copied from ccpp_memory.F90 for non-CCPP builds

      contains

         function ccpp_memory_usage(mpicomm, memory_usage) result(ierr)

            implicit none

            ! Interface variables
            integer, intent(in)                          :: mpicomm
            character(len=*), intent(out)                :: memory_usage
            ! Function return value
            integer                                      :: ierr
            ! Local variables
            character(len=len(memory_usage),kind=c_char) :: memory_usage_c
            integer                                      :: i

            ierr = ccpp_memory_usage_c(mpicomm, memory_usage_c, len(memory_usage_c))
            if (ierr /= 0) then
                write(memory_usage,fmt='(a)') "An error occurred in the call to ccpp_memory_usage_c in ccpp_memory_usage"
                return
            end if

            memory_usage = memory_usage_c(1:index(memory_usage_c, c_null_char)-1)

         end function ccpp_memory_usage

#endif

      end subroutine memcheck_run

    end module memcheck
