!> \file FV3_test.F90

    module FV3_test

      private
 
      public FV3_test_init, FV3_test_run, FV3_test_finalize

      contains

!> \section arg_table_FV3_test_init Argument Table
!! | local_name | standard_name              | long_name                                      | units | rank | type      |    kind   | intent | optional |
!! |------------|----------------------------|------------------------------------------------|-------|------|-----------|-----------|--------|----------|
!! | errmsg     | error_message              | error message for error handling in CCPP       | none  |    0 | character | len=*     | out    | F        |
!! | errflg     | error_flag                 | error flag for error handling in CCPP          | flag  |    0 | integer   |           | out    | F        |
!!
      subroutine FV3_test_init (errmsg, errflg)

         implicit none

         !--- interface variables
         character(len=*), intent(inout) :: errmsg
         integer,          intent(  out) :: errflg

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

      end subroutine FV3_test_init

!> \section arg_table_FV3_test_finalize Argument Table
!! | local_name | standard_name              | long_name                                      | units | rank | type      |    kind   | intent | optional |
!! |------------|----------------------------|------------------------------------------------|-------|------|-----------|-----------|--------|----------|
!! | errmsg     | error_message              | error message for error handling in CCPP       | none  |    0 | character | len=*     | out    | F        |
!! | errflg     | error_flag                 | error flag for error handling in CCPP          | flag  |    0 | integer   |           | out    | F        |
!!
      subroutine FV3_test_finalize (errmsg, errflg)

         implicit none

         !--- interface variables
         character(len=*), intent(  out) :: errmsg
         integer,          intent(  out) :: errflg

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

      end subroutine FV3_test_finalize

!> \section arg_table_FV3_test_run Argument Table
!! | local_name | standard_name              | long_name                                      | units | rank | type      |    kind   | intent | optional |
!! |------------|----------------------------|------------------------------------------------|-------|------|-----------|-----------|--------|----------|
!! | mpirank    | mpi_rank                   | current MPI-rank                               | index |    0 | integer   |           | in     | F        |
!! | dummy      | FV3_ccpp_integration_dummy | dummy variable to test CCPP integration in FV3 | none  |    0 | integer   |           | inout  | F        |
!! | errmsg     | error_message              | error message for error handling in CCPP       | none  |    0 | character | len=*     | out    | F        |
!! | errflg     | error_flag                 | error flag for error handling in CCPP          | flag  |    0 | integer   |           | out    | F        |
!!
      subroutine FV3_test_run (mpirank, dummy, errmsg, errflg)

#ifdef OPENMP
         use omp_lib
#endif

         implicit none

         !--- interface variables
         integer,          intent(in)    :: mpirank
         integer,          intent(inout) :: dummy
         character(len=*), intent(inout) :: errmsg
         integer,          intent(  out) :: errflg

         !--- local variables
         integer :: ompthread

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

#ifdef OPENMP
         ompthread = OMP_GET_THREAD_NUM()
#else
         ompthread = 0
#endif

         dummy = dummy + ((mpirank+1)*10) + (ompthread+1)

         if (mpirank==0 .and. ompthread==0) then
            write(0,'(a,i0)') 'Called FV3_test_run and set dummy=', dummy
         end if

      end subroutine FV3_test_run

    end module FV3_test
