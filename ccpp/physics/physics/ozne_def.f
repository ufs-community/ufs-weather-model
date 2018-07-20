      module ozne_def
      use machine , only : kind_phys
      implicit none
      
      integer, parameter :: kozpl=28, kozc=48

      integer latsozp, levozp, timeoz, latsozc, levozc, timeozc
     &,       oz_coeff
      real (kind=kind_phys) blatc, dphiozc
      real (kind=kind_phys), allocatable :: oz_lat(:), oz_pres(:)
     &,                                     oz_time(:)
      real (kind=kind_phys), allocatable :: ozplin(:,:,:,:)

      end module ozne_def
