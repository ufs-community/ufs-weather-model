      function num_parthds()
#ifdef _OPENMP
#include <omp_lib.h>
!$OMP PARALLEL
      num_parthds=omp_get_num_threads()
!$OMP END PARALLEL
#else
!     num_parthds=8
      num_parthds=1
#endif
      return
      end

!GFDL      function num_parthds()
!GFDL      integer:: number_of_openMP_threads
!GFDL      character(2) :: omp_threads
!GFDL      integer :: stat
!GFDL      call get_environment_variable("OMP_NUM_THREADS",omp_threads)
!GFDL      read(omp_threads,*,iostat=stat)number_of_openMP_threads
!GFDL      num_parthds = number_of_openMP_threads
!GFDL      return
!GFDL      end 

