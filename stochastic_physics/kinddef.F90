module kinddef

      implicit none

      private

      public :: kind_evod, kind_phys
      public :: kind_dbl_prec, kind_qdt_prec
      public :: kind_io4, kind_io8

      integer, parameter :: kind_io4      = 4

      ! DH* TODO - stochastic physics / CA should be using only one of these
      integer, parameter :: kind_evod     = 8
      integer, parameter :: kind_phys     = 8
      integer, parameter :: kind_dbl_prec = 8
      integer, parameter :: kind_io8      = 8

#ifdef NO_QUAD_PRECISION
      integer, parameter :: kind_qdt_prec = 8
#else
      integer, parameter :: kind_qdt_prec = 16
#endif

end module kinddef
