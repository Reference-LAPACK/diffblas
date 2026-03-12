! Test program for DNRM2 vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirs=4

program test_dnrm2_vector_forward
  implicit none
  integer, parameter :: nbdirs = 4

  real(8), external :: dnrm2
  external :: dnrm2_dv

  ! Test parameters
  integer :: n  ! Current size (set in loop)
  integer, parameter :: max_size = 100  ! Maximum array dimension (multi-size: 1,4,40,100)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, idir  ! Loop counters
  integer :: test_sizes(1), itest
  logical :: passed, all_passed
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for initialization

  integer :: nsize
  real(8), dimension(max_size) :: x
  integer :: incx_val

  ! Vector mode derivative variables (type-promoted)
  ! Scalars become arrays(nbdirs), arrays gain extra dimension
  real(8), dimension(nbdirs,max_size) :: x_dv
  ! Declare variables for storing original values
  real(8), dimension(max_size) :: x_orig
  real(8), dimension(nbdirs,max_size) :: x_dv_orig

  ! Function result variables
  real(8) :: dnrm2_result
  real(8), dimension(nbdirs) :: dnrm2_dv_result

  test_sizes = (/ 4 /)
  write(*,*) 'Testing DNRM2 (Vector Forward, multi-size: n = 4)'
  all_passed = .true.
  do itest = 1, 1
    n = test_sizes(itest)
    write(*,*) 'Testing DNRM2 (Vector Forward, n =', n, ')'

    call run_test_for_size(n, passed)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) then
    write(*,*) 'PASS: Vector forward mode - all sizes completed successfully'
  else
    write(*,*) 'FAIL: Vector forward mode - one or more sizes had derivative errors'
  end if

contains

  subroutine run_test_for_size(n, passed)
    implicit none
    integer, intent(in) :: n
    logical, intent(out) :: passed

    ! Initialize test parameters
    nsize = n
    incx_val = 1

    ! Initialize test data with random numbers
    ! Initialize random seed for reproducible results
    seed_array = 42
    call random_seed(put=seed_array)

    call random_number(x)
    x = x * 2.0 - 1.0  ! Scale to [-1,1]

    ! Initialize input derivatives to random values (exactly like scalar mode)
    do idir = 1, nbdirs
      call random_number(x_dv(idir,:))
      x_dv(idir,:) = x_dv(idir,:) * 2.0 - 1.0
    end do

    write(*,*) 'Testing DNRM2 (Vector Forward Mode)'
    ! Store original values before any function calls
    x_orig = x
    x_dv_orig = x_dv

    ! Call the vector mode differentiated function
    call dnrm2_dv(nsize, x, x_dv, incx_val, dnrm2_result, dnrm2_dv_result, nbdirs)

    ! Numerical differentiation check
    call check_derivatives_numerically(passed)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(passed)
    implicit none
    logical, intent(out) :: passed
    real(8), parameter :: h = 1.0e-7  ! Step size for finite differences
    real(8) :: relative_error, max_error
    real(8) :: abs_error, abs_reference, error_bound
    real(8) :: central_diff, ad_result
    integer :: i, j, idir
    logical :: has_large_errors
    real(8) :: dnrm2_forward, dnrm2_backward

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Number of directions:', nbdirs

    ! Test each derivative direction separately
    do idir = 1, nbdirs

      ! Forward perturbation: f(x + h * direction)
      x = x_orig + h * x_dv_orig(idir,:)
      dnrm2_forward = dnrm2(nsize, x, incx_val)

      ! Backward perturbation: f(x - h * direction)
      x = x_orig - h * x_dv_orig(idir,:)
      dnrm2_backward = dnrm2(nsize, x, incx_val)

      ! Central difference and AD comparison
      central_diff = (dnrm2_forward - dnrm2_backward) / (2.0e0 * h)
      ad_result = dnrm2_dv_result(idir)
      abs_error = abs(central_diff - ad_result)
      abs_reference = abs(ad_result)
      error_bound = 1.0e-5 + 1.0e-5 * abs_reference
      if (abs_error > error_bound) then
        has_large_errors = .true.
      end if
      relative_error = abs_error / max(abs_reference, 1.0e-10)
      max_error = max(max_error, relative_error)
    end do

    write(*,*) 'Maximum relative error across all directions:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in vector derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Vector derivatives are within tolerance (rtol + atol)'
    end if
  end subroutine check_derivatives_numerically

end program test_dnrm2_vector_forward