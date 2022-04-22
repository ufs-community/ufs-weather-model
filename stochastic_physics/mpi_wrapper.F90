module mpi_wrapper

   implicit none

   private

   public :: mype, npes, root, comm, is_rootpe
   public :: mpi_wrapper_initialize, mpi_wrapper_finalize
   public :: mp_reduce_min, mp_reduce_max, mp_reduce_sum
   public :: mp_bcst, mp_alltoall

#include "mpif.h"

   integer, save :: mype = -999
   integer, save :: npes = -999
   integer, save :: root = -999
   integer, save :: comm = -999
   logical, save :: initialized = .false.

   integer :: ierror

   INTERFACE mp_bcst
     MODULE PROCEDURE mp_bcst_i
     MODULE PROCEDURE mp_bcst_r4
     MODULE PROCEDURE mp_bcst_r8
     MODULE PROCEDURE mp_bcst_1d_r4
     MODULE PROCEDURE mp_bcst_1d_r8
     MODULE PROCEDURE mp_bcst_2d_r4
     MODULE PROCEDURE mp_bcst_2d_r8
     MODULE PROCEDURE mp_bcst_3d_r4
     MODULE PROCEDURE mp_bcst_3d_r8
     MODULE PROCEDURE mp_bcst_4d_r4
     MODULE PROCEDURE mp_bcst_4d_r8
     MODULE PROCEDURE mp_bcst_1d_i
     MODULE PROCEDURE mp_bcst_2d_i
     MODULE PROCEDURE mp_bcst_3d_i
     MODULE PROCEDURE mp_bcst_4d_i
   END INTERFACE

   INTERFACE mp_reduce_min
     MODULE PROCEDURE mp_reduce_min_r4
     MODULE PROCEDURE mp_reduce_min_r8
   END INTERFACE

   INTERFACE mp_reduce_max
     MODULE PROCEDURE mp_reduce_max_r4_1d
     MODULE PROCEDURE mp_reduce_max_r4
     MODULE PROCEDURE mp_reduce_max_r8_1d
     MODULE PROCEDURE mp_reduce_max_r8
     MODULE PROCEDURE mp_reduce_max_i
   END INTERFACE

   INTERFACE mp_reduce_sum
     MODULE PROCEDURE mp_reduce_sum_r4
     MODULE PROCEDURE mp_reduce_sum_r4_1d
     MODULE PROCEDURE mp_reduce_sum_r4_1darr
     MODULE PROCEDURE mp_reduce_sum_r4_2darr
     MODULE PROCEDURE mp_reduce_sum_r8
     MODULE PROCEDURE mp_reduce_sum_r8_1d
     MODULE PROCEDURE mp_reduce_sum_r8_1darr
     MODULE PROCEDURE mp_reduce_sum_r8_2darr
     MODULE PROCEDURE mp_reduce_sum_i
     MODULE PROCEDURE mp_reduce_sum_i8
   END INTERFACE

   INTERFACE mp_alltoall
     MODULE PROCEDURE mp_alltoall_r4_1darr
   END INTERFACE

contains

   logical function is_rootpe()
      if (mype==root) then
         is_rootpe = .true.
      else
         is_rootpe = .false.
      end if
   end function is_rootpe

   subroutine mpi_wrapper_initialize(mpiroot, mpicomm)
      integer, intent(in) :: mpiroot, mpicomm
      if (initialized) return
      root = mpiroot
      comm = mpicomm
      call MPI_COMM_RANK(comm, mype, ierror)
      call MPI_COMM_SIZE(comm, npes, ierror)
      initialized = .true.
   end subroutine mpi_wrapper_initialize

   subroutine mpi_wrapper_finalize()
      if (.not.initialized) return
      mype = -999
      npes = -999
      root = -999
      comm = -999
      initialized = .false.
   end subroutine mpi_wrapper_finalize

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_bcst_i :: Call SPMD broadcast 
!
      subroutine mp_bcst_i(q)
         integer, intent(INOUT)  :: q

         call MPI_BCAST(q, 1, MPI_INTEGER, root, comm, ierror)

      end subroutine mp_bcst_i
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_bcst_r4 :: Call SPMD broadcast 
!
      subroutine mp_bcst_r4(q)
         real(kind=4), intent(INOUT)  :: q

         call MPI_BCAST(q, 1, MPI_REAL, root, comm, ierror)

      end subroutine mp_bcst_r4
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_bcst_r8 :: Call SPMD broadcast 
!
      subroutine mp_bcst_r8(q)
         real(kind=8), intent(INOUT)  :: q

         call MPI_BCAST(q, 1, MPI_DOUBLE_PRECISION, root, comm, ierror)

      end subroutine mp_bcst_r8
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_bcst_1d_r4 :: Call SPMD broadcast 
!
      subroutine mp_bcst_1d_r4(q, idim)
         integer, intent(IN)  :: idim
         real(kind=4), intent(INOUT)  :: q(idim)

         call MPI_BCAST(q, idim, MPI_REAL, root, comm, ierror)

      end subroutine mp_bcst_1d_r4
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_bcst_1d_r8 :: Call SPMD broadcast 
!
      subroutine mp_bcst_1d_r8(q, idim)
         integer, intent(IN)  :: idim
         real(kind=8), intent(INOUT)  :: q(idim)

         call MPI_BCAST(q, idim, MPI_DOUBLE_PRECISION, root, comm, ierror)

      end subroutine mp_bcst_1d_r8
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_bcst_2d_r4 :: Call SPMD broadcast 
!
      subroutine mp_bcst_2d_r4(q, idim, jdim)
         integer, intent(IN)  :: idim, jdim
         real(kind=4), intent(INOUT)  :: q(idim,jdim)

         call MPI_BCAST(q, idim*jdim, MPI_REAL, root, comm, ierror)

      end subroutine mp_bcst_2d_r4
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_bcst_2d_r8 :: Call SPMD broadcast 
!
      subroutine mp_bcst_2d_r8(q, idim, jdim)
         integer, intent(IN)  :: idim, jdim
         real(kind=8), intent(INOUT)  :: q(idim,jdim)

         call MPI_BCAST(q, idim*jdim, MPI_DOUBLE_PRECISION, root, comm, ierror)

      end subroutine mp_bcst_2d_r8
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_bcst_3d_r4 :: Call SPMD broadcast 
!
      subroutine mp_bcst_3d_r4(q, idim, jdim, kdim)
         integer, intent(IN)  :: idim, jdim, kdim
         real(kind=4), intent(INOUT)  :: q(idim,jdim,kdim)

         call MPI_BCAST(q, idim*jdim*kdim, MPI_REAL, root, comm, ierror)

      end subroutine mp_bcst_3d_r4
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_bcst_3d_r8 :: Call SPMD broadcast 
!
      subroutine mp_bcst_3d_r8(q, idim, jdim, kdim)
         integer, intent(IN)  :: idim, jdim, kdim
         real(kind=8), intent(INOUT)  :: q(idim,jdim,kdim)

         call MPI_BCAST(q, idim*jdim*kdim, MPI_DOUBLE_PRECISION, root, comm, ierror)

      end subroutine mp_bcst_3d_r8
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!       
!     mp_bcst_4d_r4 :: Call SPMD broadcast 
!
      subroutine mp_bcst_4d_r4(q, idim, jdim, kdim, ldim)
         integer, intent(IN)  :: idim, jdim, kdim, ldim
         real(kind=4), intent(INOUT)  :: q(idim,jdim,kdim,ldim)

         call MPI_BCAST(q, idim*jdim*kdim*ldim, MPI_REAL, root, comm, ierror)
     
      end subroutine mp_bcst_4d_r4
!     
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!       
!     mp_bcst_4d_r8 :: Call SPMD broadcast 
!
      subroutine mp_bcst_4d_r8(q, idim, jdim, kdim, ldim)
         integer, intent(IN)  :: idim, jdim, kdim, ldim
         real(kind=8), intent(INOUT)  :: q(idim,jdim,kdim,ldim)

         call MPI_BCAST(q, idim*jdim*kdim*ldim, MPI_DOUBLE_PRECISION, root, comm, ierror)
     
      end subroutine mp_bcst_4d_r8
!     
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_bcst_3d_i :: Call SPMD broadcast
!
      subroutine mp_bcst_3d_i(q, idim, jdim, kdim)
         integer, intent(IN)  :: idim, jdim, kdim
         integer, intent(INOUT)  :: q(idim,jdim,kdim)

         call MPI_BCAST(q, idim*jdim*kdim, MPI_INTEGER, root, comm, ierror)

      end subroutine mp_bcst_3d_i
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_bcst_1d_i :: Call SPMD broadcast
!
      subroutine mp_bcst_1d_i(q, idim)
         integer, intent(IN)  :: idim
         integer, intent(INOUT)  :: q(idim)

         call MPI_BCAST(q, idim, MPI_INTEGER, root, comm, ierror)

      end subroutine mp_bcst_1d_i
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_bcst_2d_i :: Call SPMD broadcast
!
      subroutine mp_bcst_2d_i(q, idim, jdim)
         integer, intent(IN)  :: idim, jdim
         integer, intent(INOUT)  :: q(idim,jdim)

         call MPI_BCAST(q, idim*jdim, MPI_INTEGER, root, comm, ierror)

      end subroutine mp_bcst_2d_i
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_bcst_4d_i :: Call SPMD broadcast
!
      subroutine mp_bcst_4d_i(q, idim, jdim, kdim, ldim)
         integer, intent(IN)  :: idim, jdim, kdim, ldim
         integer, intent(INOUT)  :: q(idim,jdim,kdim,ldim)

         call MPI_BCAST(q, idim*jdim*kdim*ldim, MPI_INTEGER, root, comm, ierror)

      end subroutine mp_bcst_4d_i
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!       
!     mp_reduce_max_r4_1d :: Call SPMD REDUCE_MAX 
!
      subroutine mp_reduce_max_r4_1d(mymax,npts)
         integer, intent(IN)  :: npts
         real(kind=4), intent(INOUT)  :: mymax(npts)
        
         real(kind=4) :: gmax(npts)
        
         call MPI_ALLREDUCE( mymax, gmax, npts, MPI_REAL, MPI_MAX, &
                             comm, ierror )
      
         mymax = gmax
        
      end subroutine mp_reduce_max_r4_1d
!     
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!       
!     mp_reduce_max_r8_1d :: Call SPMD REDUCE_MAX 
!
      subroutine mp_reduce_max_r8_1d(mymax,npts)
         integer, intent(IN)  :: npts
         real(kind=8), intent(INOUT)  :: mymax(npts)
        
         real(kind=8) :: gmax(npts)
        
         call MPI_ALLREDUCE( mymax, gmax, npts, MPI_DOUBLE_PRECISION, MPI_MAX, &
                             comm, ierror )
      
         mymax = gmax
        
      end subroutine mp_reduce_max_r8_1d
!     
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_reduce_max_r4 :: Call SPMD REDUCE_MAX 
!
      subroutine mp_reduce_max_r4(mymax)
         real(kind=4), intent(INOUT)  :: mymax

         real(kind=4) :: gmax

         call MPI_ALLREDUCE( mymax, gmax, 1, MPI_REAL, MPI_MAX, &
                             comm, ierror )

         mymax = gmax

      end subroutine mp_reduce_max_r4

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_reduce_max_r8 :: Call SPMD REDUCE_MAX 
!
      subroutine mp_reduce_max_r8(mymax)
         real(kind=8), intent(INOUT)  :: mymax

         real(kind=8) :: gmax

         call MPI_ALLREDUCE( mymax, gmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                             comm, ierror )

         mymax = gmax

      end subroutine mp_reduce_max_r8

      subroutine mp_reduce_min_r4(mymin)
         real(kind=4), intent(INOUT)  :: mymin

         real(kind=4) :: gmin

         call MPI_ALLREDUCE( mymin, gmin, 1, MPI_REAL, MPI_MIN, &
                             comm, ierror )

         mymin = gmin

      end subroutine mp_reduce_min_r4

      subroutine mp_reduce_min_r8(mymin)
         real(kind=8), intent(INOUT)  :: mymin

         real(kind=8) :: gmin

         call MPI_ALLREDUCE( mymin, gmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                             comm, ierror )

         mymin = gmin

      end subroutine mp_reduce_min_r8
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_bcst_4d_i :: Call SPMD REDUCE_MAX 
!
      subroutine mp_reduce_max_i(mymax)
         integer, intent(INOUT)  :: mymax

         integer :: gmax

         call MPI_ALLREDUCE( mymax, gmax, 1, MPI_INTEGER, MPI_MAX, &
                             comm, ierror )

         mymax = gmax

      end subroutine mp_reduce_max_i
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_reduce_sum_r4 :: Call SPMD REDUCE_SUM 
!
      subroutine mp_reduce_sum_r4(mysum)
         real(kind=4), intent(INOUT)  :: mysum

         real(kind=4) :: gsum

         call MPI_ALLREDUCE( mysum, gsum, 1, MPI_REAL, MPI_SUM, &
                             comm, ierror )

         mysum = gsum

      end subroutine mp_reduce_sum_r4
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_reduce_sum_r8 :: Call SPMD REDUCE_SUM 
!
      subroutine mp_reduce_sum_r8(mysum)
         real(kind=8), intent(INOUT)  :: mysum

         real(kind=8) :: gsum

         call MPI_ALLREDUCE( mysum, gsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                             comm, ierror )

         mysum = gsum

      end subroutine mp_reduce_sum_r8
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! !
!
!     mp_reduce_sum_r4_1darr :: Call SPMD REDUCE_SUM
!
      subroutine mp_reduce_sum_r4_1darr(mysum, npts)
         integer, intent(in)  :: npts
         real(kind=4), intent(inout)  :: mysum(npts)
         real(kind=4)                 :: gsum(npts)

         gsum = 0.0
         call MPI_ALLREDUCE( mysum, gsum, npts, MPI_REAL, MPI_SUM, &
                             comm, ierror )

         mysum = gsum

      end subroutine mp_reduce_sum_r4_1darr
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! !
!
!     mp_reduce_sum_r4_2darr :: Call SPMD REDUCE_SUM
!
      subroutine mp_reduce_sum_r4_2darr(mysum, npts1,npts2)
         integer, intent(in)  :: npts1,npts2
         real(kind=4), intent(inout)  :: mysum(npts1,npts2)
         real(kind=4)                 :: gsum(npts1,npts2)

         gsum = 0.0
         call MPI_ALLREDUCE( mysum, gsum, npts1*npts2, MPI_REAL, MPI_SUM, &
                             comm, ierror )

         mysum = gsum

      end subroutine mp_reduce_sum_r4_2darr
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_reduce_sum_r4_1d :: Call SPMD REDUCE_SUM 
!
      subroutine mp_reduce_sum_r4_1d(mysum, sum1d, npts)
         integer, intent(in)  :: npts
         real(kind=4), intent(in)     :: sum1d(npts)
         real(kind=4), intent(INOUT)  :: mysum

         real(kind=4) :: gsum
         integer :: i

         mysum = 0.0
         do i=1,npts
            mysum = mysum + sum1d(i)
         enddo 

         call MPI_ALLREDUCE( mysum, gsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                             comm, ierror )

         mysum = gsum

      end subroutine mp_reduce_sum_r4_1d
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_reduce_sum_r8_1d :: Call SPMD REDUCE_SUM 
!
      subroutine mp_reduce_sum_r8_1d(mysum, sum1d, npts)
         integer, intent(in)  :: npts
         real(kind=8), intent(in)     :: sum1d(npts)
         real(kind=8), intent(INOUT)  :: mysum

         real(kind=8) :: gsum
         integer :: i

         mysum = 0.0
         do i=1,npts
            mysum = mysum + sum1d(i)
         enddo 

         call MPI_ALLREDUCE( mysum, gsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                             comm, ierror )

         mysum = gsum

      end subroutine mp_reduce_sum_r8_1d
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! !
!
!     mp_reduce_sum_r8_1darr :: Call SPMD REDUCE_SUM
!
      subroutine mp_reduce_sum_r8_1darr(mysum, npts)
         integer, intent(in)  :: npts
         real(kind=8), intent(inout)  :: mysum(npts)
         real(kind=8)                 :: gsum(npts)

         gsum = 0.0

         call MPI_ALLREDUCE( mysum, gsum, npts, MPI_DOUBLE_PRECISION, &
                             MPI_SUM,                                 &
                             comm, ierror )

         mysum = gsum

      end subroutine mp_reduce_sum_r8_1darr
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! !

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! !
!
!     mp_reduce_sum_r8_2darr :: Call SPMD REDUCE_SUM
!
      subroutine mp_reduce_sum_r8_2darr(mysum, npts1,npts2)
         integer, intent(in)  :: npts1,npts2
         real(kind=8), intent(inout)  :: mysum(npts1,npts2)
         real(kind=8)                 :: gsum(npts1,npts2)

         gsum = 0.0

         call MPI_ALLREDUCE( mysum, gsum, npts1*npts2,      &
                             MPI_DOUBLE_PRECISION, MPI_SUM, &
                             comm, ierror )

         mysum = gsum

      end subroutine mp_reduce_sum_r8_2darr
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! !

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_reduce_sum_i :: Call SPMD REDUCE_SUM 
!
      subroutine mp_reduce_sum_i(mysum)
         integer, intent(INOUT)  :: mysum

         integer :: gsum

         call MPI_ALLREDUCE( mysum, gsum, 1, MPI_INTEGER, MPI_SUM, &
                             comm, ierror )

         mysum = gsum

      end subroutine mp_reduce_sum_i
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     mp_reduce_sum_i8 :: Call SPMD REDUCE_SUM 
!
      subroutine mp_reduce_sum_i8(mysum)
         integer*8, intent(INOUT)  :: mysum

         integer*8 :: gsum

         call MPI_ALLREDUCE( mysum, gsum, 1, MPI_INTEGER8, MPI_SUM, &
                             comm, ierror )

         mysum = gsum

      end subroutine mp_reduce_sum_i8
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! !
!
!     mp_reduce_sum_r8_2darr :: Call SPMD REDUCE_SUM
!
      subroutine mp_alltoall_r4_1darr(sbuf, ssize, sdispl, rbuf, rsize, rdispl)
         real(kind=4), intent(in)    :: sbuf(:)
         real(kind=4), intent(inout) :: rbuf(:)
         integer, intent(in) :: ssize(:),  rsize(:)
         integer, intent(in) :: sdispl(:), rdispl(:)

         call MPI_ALLTOALLV( sbuf, ssize, sdispl, MPI_REAL, &
                             rbuf, rsize, rdispl, MPI_REAL, &
                             comm, ierror )

      end subroutine mp_alltoall_r4_1darr
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! !

end module mpi_wrapper
