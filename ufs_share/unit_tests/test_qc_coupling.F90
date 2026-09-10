program test_qc_coupling
  use ESMF
  use qc_coupling_mod
  implicit none

  integer                       :: rc
  real(ESMF_KIND_R8)            :: min_val, max_val
  logical                       :: is_known
  real(ESMF_KIND_R8)            :: test_array(3)
  real(ESMF_KIND_R8), parameter :: TOLERANCE = 1.0e-5_ESMF_KIND_R8

  integer                       :: thread_errors, pet_id, pet_count
  real(ESMF_KIND_R8)            :: thread_array(3)
  type(ESMF_VM)                 :: vm

  print *, "--- Starting QC Coupling Unit Tests (Namelist Version) ---"

  call ESMF_Initialize(rc=rc)

  print *, "Testing Namelist Initialization..."
  call init_qc_registry('test_qc.nml', rc)
  if (rc /= ESMF_SUCCESS) then
     print *, "FAIL: Could not parse test_qc.nml"
     stop 1
  end if

  ! ---------------------------------------------------------
  ! TEST SUITE A: Namelist Parsing & Binary Search Verification
  ! ---------------------------------------------------------
  print *, "Testing Namelist parsing and usage..."

  call get_limits('test_temp', min_val, max_val, is_known)
  if (.not. is_known) stop 10
  if (abs(min_val - 150.5_ESMF_KIND_R8) > TOLERANCE) stop 11
  if (abs(max_val - 350.5_ESMF_KIND_R8) > TOLERANCE) stop 12
  print *, "  PASS: test_temp parsed correctly."

  call get_limits('test_multiline', min_val, max_val, is_known)
  if (.not. is_known) stop 20
  if (abs(min_val - (-10.0_ESMF_KIND_R8)) > TOLERANCE) stop 21
  if (abs(max_val - 10.0_ESMF_KIND_R8) > TOLERANCE) stop 22
  print *, "  PASS: test_multiline (spaced format) parsed correctly."

  call get_limits('missing_var', min_val, max_val, is_known)
  if (is_known) stop 30
  print *, "  PASS: missing_var handled correctly."

  print *, "Testing internal Quicksort routine..."
  
  ! If the array wasn't sorted perfectly, the Binary Search will fail to find these
  call get_limits('Z_var', min_val, max_val, is_known)
  if (.not. is_known) stop 31
  
  call get_limits('A_var_1', min_val, max_val, is_known)
  if (.not. is_known) stop 32
  
  call get_limits('A_var_2', min_val, max_val, is_known)
  if (.not. is_known) stop 33
  
  call get_limits('X_var', min_val, max_val, is_known)
  if (.not. is_known) stop 34

  print *, "  PASS: Quicksort correctly alphabetized reversed and similar strings."

  ! ---------------------------------------------------------
  ! TEST SUITE B: Array Clamping Execution
  ! ---------------------------------------------------------
  print *, "Testing QC Clamping..."
  
  test_array(1) = 200.0_ESMF_KIND_R8 ! Valid
  test_array(2) = 50.0_ESMF_KIND_R8  ! Too cold
  test_array(3) = 400.0_ESMF_KIND_R8 ! Too hot

  call apply_qc_check('test_temp', test_array, QC_ACTION_CLAMP, rc)

  if (abs(test_array(1) - 200.0_ESMF_KIND_R8) > TOLERANCE) stop 50
  if (abs(test_array(2) - 150.5_ESMF_KIND_R8) > TOLERANCE) stop 51
  if (abs(test_array(3) - 350.5_ESMF_KIND_R8) > TOLERANCE) stop 52

  print *, "  PASS: Arrays safely bounded to physical limits."
  print *, "SUCCESS: All QC Coupling Namelist tests passed!"

  ! ---------------------------------------------------------
  ! TEST SUITE C: OpenMP Thread Safety Verification
  ! ---------------------------------------------------------
  print *, "Testing OpenMP Thread Safety..."

  thread_errors = 0

  ! Spin up a parallel region.
  ! Each thread gets its own private copy of 'thread_array' and 'rc'.
  ! They all share 'thread_errors' and safely add to it using a reduction.
  !$omp parallel private(thread_array, rc) reduction(+:thread_errors)

  !$omp master
  print *, "  Parallel region initialized."
  !$omp end master

  ! Each thread sets up a faulty array
  thread_array(1) = 200.0_ESMF_KIND_R8 ! Valid
  thread_array(2) = 50.0_ESMF_KIND_R8  ! Too cold
  thread_array(3) = 400.0_ESMF_KIND_R8 ! Too hot

  ! Threads enter registry and logger
  call apply_qc_check('test_temp', thread_array, QC_ACTION_CLAMP, rc)

  ! 3. Each thread checks its own math
  if (abs(thread_array(1) - 200.0_ESMF_KIND_R8) > TOLERANCE) thread_errors = thread_errors + 1
  if (abs(thread_array(2) - 150.5_ESMF_KIND_R8) > TOLERANCE) thread_errors = thread_errors + 1
  if (abs(thread_array(3) - 350.5_ESMF_KIND_R8) > TOLERANCE) thread_errors = thread_errors + 1

  !$omp end parallel

  if (thread_errors > 0) then
     print *, "FAIL: Threading collisions detected. Errors: ", thread_errors
     stop 60
  end if

  ! ---------------------------------------------------------
  ! TEST SUITE D: MPI Parallel Verification
  ! ---------------------------------------------------------
  print *, "Testing MPI Parallel Execution..."

  ! Get the local MPI Rank (PET) from ESMF
  call ESMF_VMGetGlobal(vm, rc=rc)
  call ESMF_VMGet(vm, localPet=pet_id, petCount=pet_count, rc=rc)

  ! Each rank tests a different value based on its ID
  test_array(1) = 200.0_ESMF_KIND_R8 + pet_id

  call apply_qc_check('test_temp', test_array, QC_ACTION_CLAMP, rc)

  if (rc /= ESMF_SUCCESS) then
     print *, "FAIL: MPI Rank ", pet_id, " encountered an error."
     stop 70
  end if

  print *, "  PASS: MPI Rank ", pet_id, " (of ", pet_count, ") executed successfully."
  
  call ESMF_Finalize(rc=rc)
end program test_qc_coupling
