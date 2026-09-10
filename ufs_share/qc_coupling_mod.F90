module qc_coupling_mod
  use ESMF
  implicit none
  private

  ! Public API
  public :: init_qc_registry
  public :: apply_qc_check
  public :: get_limits        ! Exposed for unit tests and physics queries
  
  ! Configurable QC Actions 
  integer, parameter, public :: QC_ACTION_IGNORE = 0
  integer, parameter, public :: QC_ACTION_WARN   = 1
  integer, parameter, public :: QC_ACTION_CLAMP  = 2
  integer, parameter, public :: QC_ACTION_ABORT  = 3

  ! Internal registry structure
  type :: qc_rule
    character(len=ESMF_MAXSTR) :: var_name
    real(ESMF_KIND_R8) :: qc_min
    real(ESMF_KIND_R8) :: qc_max
  end type qc_rule

  ! Memory-resident dictionary of QC limits
  type(qc_rule), allocatable, save :: qc_registry(:)
  
  ! Thread-local guard
  logical, save :: is_initialized = .false.

  interface apply_qc_check
    module procedure apply_qc_check_2d_r8
    module procedure apply_qc_check_1d_r8
    module procedure apply_qc_check_2d_r4
    module procedure apply_qc_check_1d_r4
  end interface

contains

  ! ==============================================================================
  ! SUBROUTINE: init_qc_registry
  ! Parses a Fortran namelist file using a safe two-pass read to extract limits.
  ! ==============================================================================
  subroutine init_qc_registry(nml_file, rc)
    character(len=*), intent(in)  :: nml_file
    integer,          intent(out) :: rc
    
    ! Single scalar variables bound to the namelist
    character(len=ESMF_MAXSTR) :: var_name
    real(ESMF_KIND_R8) :: qc_min
    real(ESMF_KIND_R8) :: qc_max

    integer :: iunit, ios, i, valid_count
    character(len=ESMF_MAXSTR) :: log_msg

    namelist /qc_rule_nml/ var_name, qc_min, qc_max

    rc = ESMF_SUCCESS
    if (is_initialized) return

    open(newunit=iunit, file=trim(nml_file), status='old', iostat=ios)
    if (ios /= 0) then
       call ESMF_LogWrite("QC FATAL: Could not open "//trim(nml_file), ESMF_LOGMSG_ERROR, rc=rc)
       rc = ESMF_RC_FILE_OPEN
       return
    end if

    ! ---------------------------------------------------------
    ! PASS 1: Count the number of valid namelist blocks
    ! ---------------------------------------------------------
    valid_count = 0
    do
      read(iunit, nml=qc_rule_nml, iostat=ios)
      if (ios < 0) exit  ! EOF reached
      if (ios > 0) then
         call ESMF_LogWrite("QC FATAL: Syntax error in QC namelist pass 1.", ESMF_LOGMSG_ERROR, rc=rc)
         rc = ESMF_RC_FILE_READ
         close(iunit)
         return
      end if
      valid_count = valid_count + 1
    end do

    if (valid_count == 0) then
      close(iunit)
      is_initialized = .true.
      return
    end if

    ! ---------------------------------------------------------
    ! PASS 2: Allocate exact memory and populate
    ! ---------------------------------------------------------
    allocate(qc_registry(valid_count))
    rewind(iunit)

    do i = 1, valid_count
      ! Reset scalars to defaults to prevent bleed-through between blocks
      var_name = ''
      qc_min   = -huge(1.0_ESMF_KIND_R8)
      qc_max   =  huge(1.0_ESMF_KIND_R8)

      read(iunit, nml=qc_rule_nml, iostat=ios)
      if (ios /= 0) then
         call ESMF_LogWrite("QC FATAL: Read error in QC namelist pass 2.", ESMF_LOGMSG_ERROR, rc=rc)
         rc = ESMF_RC_FILE_READ
         close(iunit)
         return
      end if
      
      qc_registry(i)%var_name = trim(var_name)
      qc_registry(i)%qc_min   = qc_min
      qc_registry(i)%qc_max   = qc_max
    end do

    close(iunit)

    ! Alphabetize the cleanly allocated array for binary searching
    if (valid_count > 1) call quicksort_registry(1, valid_count)

    write(log_msg, '(A,I0,A)') "QC INIT: Successfully loaded ", valid_count, " rules."
    call ESMF_LogWrite(trim(log_msg), ESMF_LOGMSG_INFO, rc=rc)

    is_initialized = .true.
  end subroutine init_qc_registry

  ! ==============================================================================
  ! PRIVATE HELPER: quicksort_registry
  ! ==============================================================================
  recursive subroutine quicksort_registry(first, last)
    integer, intent(in) :: first, last
    type(qc_rule)       :: temp
    integer             :: i, j
    character(len=ESMF_MAXSTR) :: pivot

    if (first >= last) return

    pivot = qc_registry((first + last) / 2)%var_name
    i = first
    j = last

    do while (i <= j)
      do while (qc_registry(i)%var_name < pivot)
        i = i + 1
      end do
      do while (qc_registry(j)%var_name > pivot)
        j = j - 1
      end do
      
      if (i <= j) then
        temp = qc_registry(i)
        qc_registry(i) = qc_registry(j)
        qc_registry(j) = temp
        i = i + 1
        j = j - 1
      end if
    end do

    if (first < j) call quicksort_registry(first, j)
    if (i < last)  call quicksort_registry(i, last)
  end subroutine quicksort_registry

  ! ==============================================================================
  ! PUBLIC HELPER: get_limits (O(log N) Binary Search)
  ! ==============================================================================
  subroutine get_limits(var_name, limit_min, limit_max, is_known)
    character(len=*),   intent(in)  :: var_name
    real(ESMF_KIND_R8), intent(out) :: limit_min, limit_max
    logical,            intent(out) :: is_known
    
    integer           :: low, high, mid
    character(len=ESMF_MAXSTR) :: target_name

    is_known = .false.
    if (.not. allocated(qc_registry)) return

    target_name = var_name
    low = 1
    high = size(qc_registry)

    do while (low <= high)
      mid = (low + high) / 2
      
      if (qc_registry(mid)%var_name == target_name) then
        limit_min = qc_registry(mid)%qc_min
        limit_max = qc_registry(mid)%qc_max
        is_known  = .true.
        return
      else if (qc_registry(mid)%var_name < target_name) then
        low = mid + 1
      else
        high = mid - 1
      end if
    end do
  end subroutine get_limits

  ! ==============================================================================
  ! CPP INCLUSIONS FOR GENERIC INTERFACES
  ! ==============================================================================

#define ROUTINE_NAME apply_qc_check_2d_r8
#define VAR_KIND ESMF_KIND_R8
#define VAR_RANK :,:
#include "apply_qc_check_template.inc"
#undef ROUTINE_NAME
#undef VAR_KIND
#undef VAR_RANK

#define ROUTINE_NAME apply_qc_check_1d_r8
#define VAR_KIND ESMF_KIND_R8
#define VAR_RANK :
#include "apply_qc_check_template.inc"
#undef ROUTINE_NAME
#undef VAR_KIND
#undef VAR_RANK

#define ROUTINE_NAME apply_qc_check_2d_r4
#define VAR_KIND ESMF_KIND_R4
#define VAR_RANK :,:
#include "apply_qc_check_template.inc"
#undef ROUTINE_NAME
#undef VAR_KIND
#undef VAR_RANK

#define ROUTINE_NAME apply_qc_check_1d_r4
#define VAR_KIND ESMF_KIND_R4
#define VAR_RANK :
#include "apply_qc_check_template.inc"
#undef ROUTINE_NAME
#undef VAR_KIND
#undef VAR_RANK

end module qc_coupling_mod
