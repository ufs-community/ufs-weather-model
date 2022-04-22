module random_numbers

  implicit none

contains
!> Returns a random number between 0 and 1
!! See https://arxiv.org/abs/2004.06278. Not an exact reproduction of
!"squares" because Fortran
!! doesn't have a uint64 type, and not all compilers provide integers
!with > 64 bits...
real function random_01_CB(ctr, key)
  use iso_fortran_env, only : int64
  integer, intent(in)  :: ctr !< ctr should be incremented each time you call the function
  integer, intent(in)  :: key !< key is like a seed: use a different key for each random stream
  integer(kind=int64) :: x, y, z ! Follows "Squares" naming convention

  x = (ctr + 1) * (key + 65536) ! 65536 added because keys below that don't work.
  y = (ctr + 1) * (key + 65536)
  z = y + (key + 65536)
  x = x*x + y
  x = ior(ishft(x,32),ishft(x,-32))
  x = x*x + z
  x = ior(ishft(x,32),ishft(x,-32))
  x = x*x + y
  x = ior(ishft(x,32),ishft(x,-32))
  x = x*x + z
  random_01_CB = .5*(1. + .5*real(int(ishft(x,-32)))/real(2**30))

end function
end module random_numbers
