module parse_tracers

  integer, parameter :: NO_TRACER = -99

  public get_tracer_index, NO_TRACER

CONTAINS

  function get_tracer_index (tracer_names, name, me, master, debug)

    character(len=32), intent(in) :: tracer_names(:)
    character(len=*),  intent(in) :: name
    integer,           intent(in) :: me
    integer,           intent(in) :: master
    logical,           intent(in) :: debug
    !--- local variables
    integer :: get_tracer_index
    integer :: i

    get_tracer_index = NO_TRACER

    do i=1, size(tracer_names)
       if (trim(name) == trim(tracer_names(i))) then
           get_tracer_index = i
           exit
       endif
    enddo

    if (debug .and. (me == master)) then
      if (get_tracer_index == NO_TRACER) then
        print *,' PE ',me,' tracer with name '//trim(name)//' not found'
      else
        print *,' PE ',me,' tracer FOUND:',trim(name)
      endif
    endif

    return

  end function get_tracer_index

end module parse_tracers
