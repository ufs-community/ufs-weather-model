!>@file test_utils.F90
!>@brief Common assertion and test summary functions for MOM6 outputlog unit tests
!!
!> @date 08-12-2026

!> Module containing common assertion and test summary functions
module test_utils

  use ESMF, only : ESMF_SUCCESS

  implicit none

  private

  integer, parameter :: real_kind = selected_real_kind( 6) !< 4 byte real
  integer, parameter ::  dbl_kind = selected_real_kind(12) !< 8 byte real
  integer, parameter ::  int_kind = selected_int_kind ( 6) !< 4 byte integer
  integer, parameter :: int8_kind = selected_int_kind (13) !< 8 byte integer

  public :: testsummary, addresult, esmf_err, itoa
  public :: assert_equal

  !> Wrapper type for dynamically sized string messages
  type :: msg_type
     character(len=:), allocatable :: str    !< a message string
  end type msg_type

  !> Type for tracking test execution counts, passes, fails, and messages
  type testsummary
     integer :: count = 0         !< a test count
     integer :: npass = 0         !< passing test count
     integer :: nfail = 0         !< failing test count
     logical,        allocatable :: teststatus(:)       !< test status list
     type(msg_type), allocatable :: testmessage(:)      !< test message list
     type(msg_type), allocatable :: errmessage(:)       !< test error message list
   contains
     !> Bound procedure to initialize the summary
     procedure :: init => init_summary
  end type testsummary

  !> Generic interface for asserting equality between expected and actual values
  interface assert_equal
     module procedure assert_int_scalar
     module procedure assert_real_scalar
     module procedure assert_double_scalar
     module procedure assert_double_1d
     module procedure assert_logical
  end interface assert_equal

  !> Generic interface for adding a test result to the summary
  interface addresult
     module procedure add_test_result
  end interface addresult

contains

  !> Initialize the test summary object
  !!
  !! @param[inout] this       The testsummary object to initialize
  !! @param[in]    maxtests   The maximum number of tests to allocate space for
  subroutine init_summary(this, maxtests)

    class(testsummary), intent(inout) :: this
    integer,            intent(in)    :: maxtests

    allocate(this%teststatus(maxtests))
    allocate(this%testmessage(maxtests))
    allocate(this%errmessage(maxtests))

  end subroutine init_summary

  !> Add test results to summary
  !!
  !! @param[inout]     summary     test summary
  !! @param[in]        passed      test result
  !! @param[in]        message     test message
  !! @param[in]        errormsg    test error message
  subroutine add_test_result(summary, passed, message, errormsg)

    type(testsummary), intent(inout) :: summary
    logical,           intent(in)    :: passed
    character(len=*),  intent(in)    :: message
    character(len=*),  intent(in)    :: errormsg

    summary%count = summary%count + 1
    call grow_if_needed(summary)
    if (passed) then
       summary%npass = summary%npass + 1
    else
       summary%nfail = summary%nfail + 1
    end if

    summary%teststatus(summary%count) = passed
    summary%testmessage(summary%count)%str = message
    summary%errmessage(summary%count)%str = errormsg
  end subroutine add_test_result

  !> Increase the size of the testsummary type arrays if needed
  !!
  !! @param[inout]    summary    test summary type
  subroutine grow_if_needed(summary)
    type(testsummary), intent(inout) :: summary

    logical,        allocatable :: new_status(:)
    type(msg_type), allocatable :: new_testmsg(:), new_errmsg(:)

    integer :: old_size, new_size

    old_size = 0
    if (allocated(summary%teststatus)) old_size = size(summary%teststatus)
    if (summary%count <= old_size) return
    new_size = max(old_size*2, summary%count)
    allocate(new_status(new_size))
    if (old_size > 0) new_status(1:old_size) = summary%teststatus(1:old_size)
    call move_alloc(new_status, summary%teststatus)
    allocate(new_testmsg(new_size))
    if (old_size > 0) new_testmsg(1:old_size) = summary%testmessage(1:old_size)
    call move_alloc(new_testmsg, summary%testmessage)
    allocate(new_errmsg(new_size))
    if (old_size > 0) new_errmsg(1:old_size) = summary%errmessage(1:old_size)
    call move_alloc(new_errmsg, summary%errmessage)
  end subroutine grow_if_needed

  !> Assert equality between two logical values
  !!
  !! @param[in]  actual      The actual boolean result
  !! @param[in]  expected    The expected boolean result
  !! @param[in]  msg         Description of the test
  !! @param[out] rc          Return code (true if pass, false if fail)
  !! @param[out] returnmsg   Formatted result string
  subroutine assert_logical(actual, expected, msg, rc, returnmsg)

    logical,          intent(in)  :: actual, expected
    character(len=*), intent(in)  :: msg
    logical,          intent(out) :: rc
    character(len=*), intent(out) :: returnmsg

    rc = (actual .eqv. expected)
    if (rc) then
       returnmsg = "Pass: "// trim(msg)
    else
       write(returnmsg, '(2(a,L2))') "Fail: " // trim(msg) // " | Expected ", expected, ", got ", actual
    endif
  end subroutine assert_logical

  !> Assert equality between two integer scalars
  !!
  !! @param[in]  actual      The actual integer
  !! @param[in]  expected    The expected integer
  !! @param[in]  msg         Description of the test
  !! @param[out] rc          Return code (true if pass, false if fail)
  !! @param[out] returnmsg   Formatted result string
  subroutine assert_int_scalar(actual, expected, msg, rc, returnmsg)
    integer(int_kind), intent(in)  :: actual, expected
    character(len=*),  intent(in)  :: msg
    logical,           intent(out) :: rc
    character(len=*),  intent(out) :: returnmsg

    rc = (actual == expected)
    if (rc) then
       returnmsg = "Pass: " // trim(msg)
    else
       write(returnmsg, '(2(a,i0))') "Fail: " // trim(msg) // " | Expected ", expected, ", got ", actual
    end if
  end subroutine assert_int_scalar

  !> Assert equality between two real scalars within a tolerance
  !!
  !! @param[in]  actual      The actual real value
  !! @param[in]  expected    The expected real value
  !! @param[in]  tol         The tolerance for equality
  !! @param[in]  msg         Description of the test
  !! @param[out] rc          Return code (true if pass, false if fail)
  !! @param[out] returnmsg   Formatted result string
  subroutine assert_real_scalar(actual, expected, tol, msg, rc, returnmsg)
    real(real_kind),   intent(in)  :: actual, expected, tol
    character(len=*),  intent(in)  :: msg
    logical,           intent(out) :: rc
    character(len=*),  intent(out) :: returnmsg

    rc = (abs(actual - expected) <= tol)
    if (rc) then
       returnmsg = "Pass: " // trim(msg)
    else
       write(returnmsg, '(2(a,g15.8))') "Fail: " // trim(msg) // " | Expected ", expected, ", got ", actual
    end if
  end subroutine assert_real_scalar

  !> Assert equality between two double precision scalars within a tolerance
  !!
  !! @param[in]  actual      The actual double value
  !! @param[in]  expected    The expected double value
  !! @param[in]  tol         The tolerance for equality
  !! @param[in]  msg         Description of the test
  !! @param[out] rc          Return code (true if pass, false if fail)
  !! @param[out] returnmsg   Formatted result string
  subroutine assert_double_scalar(actual, expected, tol, msg, rc, returnmsg)
    real(dbl_kind),    intent(in)  :: actual, expected, tol
    character(len=*),  intent(in)  :: msg
    logical,           intent(out) :: rc
    character(len=*),  intent(out) :: returnmsg

    rc = (abs(actual - expected) <= tol)
    if (rc) then
       returnmsg = "Pass: " // trim(msg)
    else
       write(returnmsg, '(2(a,g20.13))') "Fail: " // trim(msg) // " | Expected ", expected, ", got ", actual
    end if
  end subroutine assert_double_scalar

  !> Assert equality between two double precision 1D arrays within a tolerance
  !!
  !! @param[in]  actual      The actual double array
  !! @param[in]  expected    The expected double array
  !! @param[in]  tol         The tolerance for equality
  !! @param[in]  msg         Description of the test
  !! @param[out] rc          Return code (true if pass, false if fail)
  !! @param[out] returnmsg   Formatted result string
  subroutine assert_double_1d(actual, expected, tol, msg, rc, returnmsg)
    real(dbl_kind),    intent(in)  :: actual(:), expected(:), tol
    character(len=*),  intent(in)  :: msg
    logical,           intent(out) :: rc
    character(len=*),  intent(out) :: returnmsg

    rc = all(abs(actual - expected) <= tol)
    if (rc) then
       returnmsg = "Pass: " // trim(msg)
    else
       returnmsg = "Fail: " // trim(msg) // " | At least one element mismatched."
    end if
  end subroutine assert_double_1d

  !> Convert an integer to a string
  !!
  !! @param[in] i  The integer to convert
  !! @return       A string representation of the integer
  function itoa(i) result(s)
    integer, intent(in) :: i
    character(len=4) :: s

    write(s,'(I0)') i
  end function itoa

  !> Check an ESMF return code and abort if an error is found
  !!
  !! @param[in] rc        The ESMF return code
  !! @param[in] subname   The name of the calling subroutine for logging
  !! @param[in] context   Context message for logging
  subroutine esmf_err(rc, subname, context)
    integer,          intent(in) :: rc
    character(len=*), intent(in) :: subname
    character(len=*), intent(in) :: context

    if (rc /= ESMF_SUCCESS) then
       write(0,'(A,I0)') "FATAL ("//trim(subname)//") : "//trim(context)//": rc=", rc
       stop 99
    end if
  end subroutine esmf_err

end module test_utils
