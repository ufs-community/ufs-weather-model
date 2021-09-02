!!! This should bprobably be cleaned up or removed

module proc_bounds

  ! proc and gridcell bounds for simple  mesh with regular decomp
  ! also, control and init type
  
  implicit none

  public procbounds_type, control_init_type
  
  private

  type procbounds_type
     integer :: gridbeg, gridend ! local gridcell range
     integer :: de               ! local de
     integer :: im               ! # gridcells on de
  end type procbounds_type

  type control_init_type
     ! modelled after FV3 GFS init and control types

     ! grid, domain, and blocking
     integer           :: npx, npy
     integer           :: ntiles
     integer           :: layout(2)
     character(len=64) :: grid
     integer           :: blocksize

     ! land model params
     integer :: ivegsrc
     integer :: isot

     ! run
     logical   :: first_time  ! flag for first time step

     ! MPI stuff
     integer :: me                                ! my MPI-rank
     integer :: master                            ! master MPI-rank
    ! !integer :: tile_num                          ! tile number for this MPI rank
    ! integer :: isc                               ! starting i-index for this MPI-domain
    ! integer :: jsc                               ! starting j-index for this MPI-domain
    ! integer :: nx                                ! number of points in i-dir for this MPI rank
    ! integer :: ny                                ! number of points in j-dir for this MPI rank
    ! integer, pointer :: blksz(:)                 ! for explicit data blocking
    ! character(len=64) :: fn_nml                  ! namelist filename

  contains
    procedure :: init  => control_initialize
  end type control_init_type
  
  type(procbounds_type), public :: procbounds

contains

!  subroutine control_initialize(Model,fn_nml,me,master,tile_num,isc, jsc, nx, ny,blksz)
  subroutine control_initialize(Model)

    use fms_mod,             only: check_nml_error, close_file, file_exist
    use mpp_mod,             only: mpp_pe, mpp_root_pe
#ifdef INTERNAL_FILE_NML
    use mpp_mod,             only: input_nml_file
#else
    use fms_mod,             only: open_namelist_file
#endif

    
    implicit none

    class(control_init_type)            :: Model
    ! integer,                intent(in) :: me                                ! my MPI-rank
    ! integer,                intent(in) :: master                            ! master MPI-rank
    ! integer,                intent(in) :: tile_num                          ! tile number for this MPI rank
    ! integer,                intent(in) :: isc                               ! starting i-index for this MPI-domain
    ! integer,                intent(in) :: jsc                               ! starting j-index for this MPI-domain
    ! integer,                intent(in) :: nx                                ! number of points in i-dir for this MPI rank
    ! integer,                intent(in) :: ny                                ! number of points in j-dir for this MPI rank
    ! integer,                intent(in) :: blksz(:)                           ! for explicit data blocking
    ! character(len=64),      intent(in) :: fn_nml                            ! namelist filename

    ! namelist variables
    ! ------------------------------------------
    ! grid, domain, and blocking
    integer           :: npx = 0, npy = 0, ntiles = 0
    integer           :: layout(2) = (/0,0/)
    character(len=64) :: grid      = 'none'
    integer           :: blocksize = -1
    ! land model params
    integer :: ivegsrc  = -1
    integer :: isot     = -1


    ! for namelist read
    integer :: unit, io, ierr
    namelist /noah_nml/ grid, npx, npy, layout, ntiles, &
         blocksize, ivegsrc, isot
    
    ! -------------------------------------------
    ! read in namelist

    if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
       read(input_nml_file, nml=noah_nml, iostat=io)
       ierr = check_nml_error(io, 'noah_nml')
#else
       unit = open_namelist_file ( )
       ierr=1
       do while (ierr /= 0)
          read(unit, nml=noah_nml, iostat=io)
          ierr = check_nml_error(io,'noah_nml')
       enddo
       call close_file(unit)
#endif
    endif

    Model%grid      = grid
    Model%blocksize = blocksize
    Model%npx       = npx
    Model%npy       = npy
    Model%layout    = layout
    Model%ntiles    = ntiles
    Model%ivegsrc   = ivegsrc
    Model%isot      = isot
    !--- MPI parameters
    Model%me        =  mpp_pe()
    Model%master    =  mpp_root_pe()
!     Model%fn_nml           = fn_nml


  end subroutine control_initialize

end module proc_bounds
