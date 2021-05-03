module proc_bounds

  ! proc and gridcell bounds for simple  mesh with regular decomp

  implicit none

  public procbounds_type
  
  private

  type procbounds_type
     integer :: gridbeg, gridend ! local gridcell range
     integer :: de               ! local de
     integer :: im               ! # gridcells on de
  end type procbounds_type
  
  type(procbounds_type), public :: procbounds

contains


end module proc_bounds
