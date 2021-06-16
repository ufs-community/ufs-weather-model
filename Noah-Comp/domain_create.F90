module land_domain_mod

  use mpp_mod
  use mpp_domains_mod
  use fms_mod,             only: open_namelist_file, &
       check_nml_error, close_file, file_exist
  implicit none
  !type(domain2D) :: land_domain

  !! NEED TO GET from nml: npx,npy, layout
  !! halo = 0

  integer :: f_unit, ios, ierr

  ! namelist variables
  integer :: npx = 0, npy = 0, ntiles = 0
  integer :: layout(2) = (/0,0/)
  namelist /noah_nml/ npx, npy,layout, ntiles

  public domain_create
  
contains

  subroutine domain_create(land_domain)

    type(domain2d), intent(out) :: land_domain
    !type (land_data_type), intent(in) :: Land ! create this

    integer, allocatable :: pe_start(:), pe_end(:)
    integer :: n
    integer :: halo = 0 ! 0 for land


    call mpp_domains_init()

    call read_namelist_noah_nml()
    write(*,*) 'spit out noah nml: ', npx, npy,layout !tmp debug

    !--- define mosaic for domain
    !ntiles=6
    allocate(pe_start(ntiles))
    allocate(pe_end(ntiles))
    do n = 1, ntiles
       pe_start(n) = mpp_root_pe() + (n-1)*layout(1)*layout(2)
       pe_end(n)   = mpp_root_pe() +     n*layout(1)*layout(2)-1
    enddo

    call define_cubic_mosaic(land_domain, npx-1, npy-1, layout, pe_start, pe_end, halo)
    !write(*,*) 'some domain info: ', land_domain%pe, land_domain%ntiles  !tmp debug
    deallocate(pe_start)
    deallocate(pe_end)

  end subroutine domain_create


  subroutine read_namelist_noah_nml

    integer :: f_unit, ios, ierr

    f_unit = open_namelist_file('input.nml')
    ! Read Noah namelist                                                                                                                            
    read (f_unit,noah_nml,iostat=ios)
    ierr = check_nml_error(ios,'noah_nml')
    call close_file(f_unit)
    !#endif

  end subroutine read_namelist_noah_nml




  ! This subroutine is copied from FMS/test_fms/test_mpp_domains.F90
  ! and modified to make it simpler to use.
  ! domain_decomp in fv_mp_mod.F90 does something similar, but it does a
  ! few other unnecessary things (and requires more arguments).
  subroutine define_cubic_mosaic(domain, ni, nj, layout, pe_start, pe_end, halo)
    type(domain2d), intent(inout) :: domain
    integer,        intent(in)    :: ni, nj
    integer,        intent(in)    :: layout(:)
    integer,        intent(in)    :: pe_start(:), pe_end(:)
    integer,        intent(in)    :: halo
    !--- local variables
    integer                       :: global_indices(4,6), layout2D(2,6)
    integer, dimension(12)        :: istart1, iend1, jstart1, jend1, tile1
    integer, dimension(12)        :: istart2, iend2, jstart2, jend2, tile2
    integer                       :: ntiles, num_contact
    integer                       :: i

    ntiles = 6
    num_contact = 12
    if(size(pe_start(:)) .NE. 6 .OR. size(pe_end(:)) .NE. 6 ) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of pe_start and pe_end should be 6")
    if(size(layout) .NE. 2) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of layout should be 2")
    do i = 1, 6
       layout2D(:,i) = layout(:)
       global_indices(1,i) = 1
       global_indices(2,i) = ni
       global_indices(3,i) = 1
       global_indices(4,i) = nj
    enddo

    !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
    tile1(1) = 1;     tile2(1) = 2
    istart1(1) = ni;  iend1(1) = ni;  jstart1(1) = 1;      jend1(1) = nj
    istart2(1) = 1;   iend2(1) = 1;   jstart2(1) = 1;      jend2(1) = nj

    !--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
    tile1(2) = 1; tile2(2) = 3
    istart1(2) = 1;      iend1(2) = ni;  jstart1(2) = nj;  jend1(2) = nj
    istart2(2) = 1;      iend2(2) = 1;   jstart2(2) = nj;  jend2(2) = 1

    !--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
    tile1(3) = 1;     tile2(3) = 5
    istart1(3) = 1;   iend1(3) = 1;      jstart1(3) = 1;   jend1(3) = nj
    istart2(3) = ni;  iend2(3) = 1;      jstart2(3) = nj;  jend2(3) = nj

    !--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
    tile1(4) = 1; tile2(4) = 6
    istart1(4) = 1;      iend1(4) = ni;  jstart1(4) = 1;   jend1(4) = 1
    istart2(4) = 1;      iend2(4) = ni;  jstart2(4) = nj;  jend2(4) = nj

    !--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
    tile1(5) = 2;        tile2(5) = 3
    istart1(5) = 1;      iend1(5) = ni;  jstart1(5) = nj;  jend1(5) = nj
    istart2(5) = 1;      iend2(5) = ni;  jstart2(5) = 1;   jend2(5) = 1

    !--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
    tile1(6) = 2; tile2(6) = 4
    istart1(6) = ni;  iend1(6) = ni;  jstart1(6) = 1;      jend1(6) = nj
    istart2(6) = ni;  iend2(6) = 1;   jstart2(6) = 1;      jend2(6) = 1

    !--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
    tile1(7) = 2; tile2(7) = 6
    istart1(7) = 1;   iend1(7) = ni;  jstart1(7) = 1;   jend1(7) = 1
    istart2(7) = ni;  iend2(7) = ni;  jstart2(7) = nj;  jend2(7) = 1

    !--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
    tile1(8) = 3; tile2(8) = 4
    istart1(8) = ni;  iend1(8) = ni;  jstart1(8) = 1;      jend1(8) = nj
    istart2(8) = 1;   iend2(8) = 1;   jstart2(8) = 1;      jend2(8) = nj

    !--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
    tile1(9) = 3; tile2(9) = 5
    istart1(9) = 1;      iend1(9) = ni;  jstart1(9) = nj;  jend1(9) = nj
    istart2(9) = 1;      iend2(9) = 1;   jstart2(9) = nj;  jend2(9) = 1

    !--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
    tile1(10) = 4; tile2(10) = 5
    istart1(10) = 1;     iend1(10) = ni; jstart1(10) = nj; jend1(10) = nj
    istart2(10) = 1;     iend2(10) = ni; jstart2(10) = 1;  jend2(10) = 1

    !--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
    tile1(11) = 4; tile2(11) = 6
    istart1(11) = ni; iend1(11) = ni; jstart1(11) = 1;     jend1(11) = nj
    istart2(11) = ni; iend2(11) = 1;  jstart2(11) = 1;     jend2(11) = 1

    !--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
    tile1(12) = 5; tile2(12) = 6
    istart1(12) = ni; iend1(12) = ni; jstart1(12) = 1;     jend1(12) = nj
    istart2(12) = 1;  iend2(12) = 1;  jstart2(12) = 1;     jend2(12) = nj

    call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, &
         num_contact, tile1, tile2, istart1, iend1, jstart1, jend1, &
         istart2, iend2, jstart2, jend2, pe_start, pe_end, symmetry=.true., &
         whalo=halo, ehalo=halo, shalo=halo, nhalo=halo, &
         name='CA cubic mosaic')
  end subroutine define_cubic_mosaic

end module land_domain_mod
                    
