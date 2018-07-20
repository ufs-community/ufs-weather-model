      module h2o_def
      use machine , only : kind_phys
      implicit none

      integer, parameter :: kh2opltc=29

      integer latsh2o, levh2o, timeh2o,  h2o_coeff
      real (kind=kind_phys), allocatable :: h2o_lat(:), h2o_pres(:)
     &,                                     h2o_time(:)
      real (kind=kind_phys), allocatable :: h2oplin(:,:,:,:)

      end module h2o_def
