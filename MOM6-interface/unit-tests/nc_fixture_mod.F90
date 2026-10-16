!> @file nc_fixture_mod.F90
!> @brief NetCDF fixture to generate files with established behaviour of file creation and
!! completion. This feature creates files and then writes files with one of two observed
!! patterns:
!!
!! 1) in DATM mode, nlen flips from 0 at creation to 1 at completion; file size is ~constant
!! 2) in ATM mode, nlen>0 at creation and completion is determined by the file size increase
!!
!! For 1) we can't reproduce the observed behaviour using a variable that actually shares "time"
!! (writing any such variable immediately advances nlen) so this is a deliberate stand-in that
!! reproduces the observable nlen/fsize progression without replicating the real file's
!! internal variable structure. file_is_complete/get_file_state only look at nlen and fsize, so
!! that's sufficient for what's under test.
!!
!! Restart files test files provide use length of unlimited dim only, since restarts are written
!! complete at requested time by MOM6 directly (FMS is not used).
!!
!> @date 09-01-2026
module nc_fixture_mod

  use netcdf
  implicit none

  integer, parameter :: fixture_nx = 200         !< bulk "data" variable spatial size (ATM payload)
  integer, parameter :: fixture_pad_n = 100000   !< padding variable size -- models the DATM case's
                                                 !! already-written bulk payload
contains
  !> Create a fresh file with the schema defined but zero records written. nlen=0.
  !! Includes an extra "pad" variable that does NOT use the unlimited "time" dim
  !!
  !! @param[in]   fname                filename
  !! @param[out]  createsize_out       optional returned file size at creation
  subroutine create_schema(fname, createsize_out)
    character(len=*), intent(in)            :: fname
    integer,          intent(out), optional :: createsize_out

    integer :: ncid, dimid_t, dimid_x, dimid_pad, varid_t, varid_data, varid_pad, fsize

    call nf90_err(nf90_create(trim(fname), ior(NF90_CLOBBER, NF90_NETCDF4), ncid), "create")
    call nf90_err(nf90_def_dim(ncid, "time", NF90_UNLIMITED, dimid_t), "def_dim time")
    call nf90_err(nf90_def_dim(ncid, "x", fixture_nx, dimid_x), "def_dim x")
    call nf90_err(nf90_def_dim(ncid, "pad", fixture_pad_n, dimid_pad), "def_dim pad")
    call nf90_err(nf90_def_var(ncid, "time", NF90_DOUBLE, [dimid_t], varid_t), "def_var time")
    call nf90_err(nf90_def_var(ncid, "data", NF90_DOUBLE, [dimid_x, dimid_t], varid_data), "def_var data")
    call nf90_err(nf90_def_var(ncid, "pad", NF90_DOUBLE, [dimid_pad], varid_pad), "def_var pad")
    call nf90_err(nf90_enddef(ncid), "enddef")
    call nf90_err(nf90_sync(ncid), "sync (schema)")
    call nf90_err(nf90_close(ncid), "close (schema)")

    if (present(createsize_out)) then
      inquire(file=trim(fname), size=fsize)
      createsize_out = fsize
    end if
  end subroutine create_schema
  !> Extend the unlimited dim by writing record 1's time-coordinate value.
  !! nlen: 0 -> 1. Does NOT write the bulk data variable.
  !!
  !! @param[in]   fname            filename
  !! @param[in]   time_value       optional time value used
  subroutine write_record(fname, time_value)
    character(len=*), intent(in)           :: fname
    real(8),          intent(in), optional :: time_value

    integer :: ncid, varid_t
    real(8) :: tval

    tval = 1.0d0
    if (present(time_value)) tval = time_value

    call nf90_err(nf90_open(trim(fname), NF90_WRITE, ncid), "open for write_record")
    call nf90_err(nf90_inq_varid(ncid, "time", varid_t), "inq_varid time")
    call nf90_err(nf90_put_var(ncid, varid_t, [tval], start=[1]), "put_var time")
    call nf90_err(nf90_sync(ncid), "sync (record)")
    call nf90_err(nf90_close(ncid), "close (record)")
  end subroutine write_record
  !> Write the bulk "data" variable for record 1. Grows fsize. Requires
  !! write_record to have been called first (record 1 must already exist).
  !!
  !! @param[in]   fname            filename
  !! @param[in]   fill_value       optional fill value for data fill
  subroutine write_bulk_data(fname, fill_value)
    character(len=*), intent(in)           :: fname
    real(8),          intent(in), optional :: fill_value

    integer :: ncid, varid_data
    real(8) :: data_vals(fixture_nx)
    real(8) :: fval

    fval = 42.0d0
    if (present(fill_value)) fval = fill_value
    data_vals = fval

    call nf90_err(nf90_open(trim(fname), NF90_WRITE, ncid), "open for write_bulk_data")
    call nf90_err(nf90_inq_varid(ncid, "data", varid_data), "inq_varid data")
    call nf90_err(nf90_put_var(ncid, varid_data, data_vals, start=[1,1], count=[fixture_nx,1]), "put_var data")
    call nf90_err(nf90_sync(ncid), "sync (data)")
    call nf90_err(nf90_close(ncid), "close (data)")
  end subroutine write_bulk_data
  !> Write the "pad" variable -- deliberately not tied to the unlimited "time"
  !! dim.
  !!
  !! @param[in]   fname            filename
  !! @param[in]   fill_value       optional fill value for data fill
  subroutine write_padding(fname, fill_value)
    character(len=*), intent(in)           :: fname
    real(8),          intent(in), optional :: fill_value

    integer :: ncid, varid_pad
    real(8) :: pad_vals(fixture_pad_n)
    real(8) :: fval

    fval = 1.0d0
    if (present(fill_value)) fval = fill_value
    pad_vals = fval

    call nf90_err(nf90_open(trim(fname), NF90_WRITE, ncid), "open for write_padding")
    call nf90_err(nf90_inq_varid(ncid, "pad", varid_pad), "inq_varid pad")
    call nf90_err(nf90_put_var(ncid, varid_pad, pad_vals), "put_var pad")
    call nf90_err(nf90_sync(ncid), "sync (pad)")
    call nf90_err(nf90_close(ncid), "close (pad)")
  end subroutine write_padding
  ! -------------------------------------------------------------------------------
  ! Convenience wrappers for the four canonical states, built from the  primitives
  !--------------------------------------------------------------------------------
  !> Wrapper for DATM incomplete
  !!
  !! @param[in]   fname    filename
  subroutine make_datm_incomplete(fname)
    character(len=*), intent(in) :: fname
    call create_schema(fname)
    call write_padding(fname)
  end subroutine make_datm_incomplete
  !> Wrapper for DATM complete
  !!
  !! @param[in]   fname    filename
  subroutine make_datm_complete(fname)
    character(len=*), intent(in) :: fname
    call create_schema(fname)
    call write_padding(fname)
    call write_record(fname)
  end subroutine make_datm_complete
  !> Wrapper for ATM incomplete
  !!
  !! @param[in]   fname            filename
  !! @param[out]  createsize_out   file size at creation
  subroutine make_atm_incomplete(fname, createsize_out)
    character(len=*),  intent(in)  :: fname
    integer,           intent(out) :: createsize_out
    integer :: fsize
    call create_schema(fname)
    call write_record(fname)
    inquire(file=trim(fname), size=fsize)
    createsize_out = fsize
  end subroutine make_atm_incomplete
  !> Wrapper for ATM complete
  !!
  !! @param[in]   fname            filename
  !! @param[out]  createsize_out   file size at completion
  subroutine make_atm_complete(fname, createsize_out)
    character(len=*),  intent(in)  :: fname
    integer,           intent(out) :: createsize_out
    call make_atm_incomplete(fname, createsize_out)
    call write_bulk_data(fname)
  end subroutine make_atm_complete
  ! -------------------------------------------------------------------------------
  ! Restart-file fixtures.
  !--------------------------------------------------------------------------------
  !> Build the filename for restart using base name and  `part_index`
  !!
  !! @param[in]  base_fname   base name
  !! @param[in]  part_index   0 for the base file, 1,2,... for later parts
  function restart_part_fname(base_fname, part_index) result(fname)
    character(len=*), intent(in) :: base_fname
    integer,          intent(in) :: part_index
    character(len=256) :: fname

    if (part_index == 0) then
      fname = trim(base_fname)//'.nc'
    else
      write(fname,'(A,A,I0,A)') trim(base_fname), '_', part_index, '.nc'
    end if
  end function restart_part_fname
  !> Create a minimal restart-part fixture: just an unlimited "time" dim,
  !! optionally with record 1 written (nlen 0->1).
  !!
  !! @param[in]  fname      filename (see restart_part_fname)
  !! @param[in]  complete   .true. to also write record 1 (nlen=1);
  !!                        .false. leaves the file at nlen=0
  subroutine make_restart_fixture(fname, complete)
    character(len=*), intent(in) :: fname
    logical,          intent(in) :: complete

    integer :: ncid, dimid_t, varid_t

    call nf90_err(nf90_create(trim(fname), ior(NF90_CLOBBER, NF90_NETCDF4), ncid), "create (restart)")
    call nf90_err(nf90_def_dim(ncid, "time", NF90_UNLIMITED, dimid_t), "def_dim time (restart)")
    call nf90_err(nf90_def_var(ncid, "time", NF90_DOUBLE, [dimid_t], varid_t), "def_var time (restart)")
    call nf90_err(nf90_enddef(ncid), "enddef (restart)")
    call nf90_err(nf90_sync(ncid), "sync (restart schema)")
    call nf90_err(nf90_close(ncid), "close (restart schema)")

    if (complete) call write_record(fname)
  end subroutine make_restart_fixture
  !> Error return function for NetCDF
  !!
  !! @param[in]    ierr      error return value
  !! @param[in]    context   failure identifier
  subroutine nf90_err(ierr, context)
    integer,          intent(in) :: ierr
    character(len=*), intent(in) :: context
    if (ierr /= nf90_noerr) then
      write(0,'(A)') "FATAL (nc_fixture_mod): "//trim(context)//": "//trim(nf90_strerror(ierr))
      stop 99
    end if
  end subroutine nf90_err
end module nc_fixture_mod
