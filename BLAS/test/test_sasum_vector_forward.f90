! Test program for SASUM vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=4

program test_sasum_vector_forward
  implicit none
  integer, parameter :: nbdirs = 4

  real(4), external :: sasum
  external :: sasum_dv

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
  real(4), dimension(max_size) :: sx
  integer :: incx_val

  ! Vector mode derivative variables (type-promoted)
  ! Scalars become arrays(nbdirs), arrays gain extra dimension
  real(4), dimension(nbdirs,max_size) :: sx_dv
  ! Declare variables for storing original values
  real(4), dimension(max_size) :: sx_orig
  real(4), dimension(nbdirs,max_size) :: sx_dv_orig

  ! Function result variables
  real(4) :: sasum_result
  real(4), dimension(nbdirs) :: sasum_dv_result

  test_sizes = (/ 4 /)
  write(*,*) 'Testing SASUM (Vector Forward, multi-size: n = 4)'
  all_passed = .true.
  do itest = 1, 1
    n = test_sizes(itest)
    write(*,*) 'Testing SASUM (Vector Forward, n =', n, ')'

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

    call random_number(sx)
    sx = sx * 2.0 - 1.0  ! Scale to [-1,1]

    ! Initialize input derivatives to random values (exactly like scalar mode)
    do idir = 1, nbdirs
      call random_number(sx_dv(idir,:))
      sx_dv(idir,:) = sx_dv(idir,:) * 2.0 - 1.0
    end do

    write(*,*) 'Testing SASUM (Vector Forward Mode)'
    ! Store original values before any function calls
    sx_orig = sx
    sx_dv_orig = sx_dv

    ! Call the vector mode differentiated function
    call sasum_dv(nsize, sx, sx_dv, incx_val, sasum_result, sasum_dv_result, nbdirs)

    ! Numerical differentiation check
    call check_derivatives_numerically(passed)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(passed)
    implicit none
    logical, intent(out) :: passed
    real(4), parameter :: h = 1.0e-3  ! Step size for finite differences
    real(4) :: relative_error, max_error
    real(4) :: abs_error, abs_reference, error_bound
    real(4) :: central_diff, ad_result
    integer :: i, j, idir
    logical :: has_large_errors
    real(4) :: sasum_forward, sasum_backward

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Number of directions:', nbdirs

    ! Test each derivative direction separately
    do idir = 1, nbdirs

      ! Forward perturbation: f(x + h * direction)
      sx = sx_orig + h * sx_dv_orig(idir,:)
      sasum_forward = sasum(nsize, sx, incx_val)

      ! Backward perturbation: f(x - h * direction)
      sx = sx_orig - h * sx_dv_orig(idir,:)
      sasum_backward = sasum(nsize, sx, incx_val)

      ! Central difference and AD comparison
      central_diff = (sasum_forward - sasum_backward) / (2.0e0 * h)
      ad_result = sasum_dv_result(idir)
      abs_error = abs(central_diff - ad_result)
      abs_reference = abs(ad_result)
      error_bound = 2.0e-3 + 2.0e-3 * abs_reference
      if (abs_error > error_bound) then
        has_large_errors = .true.
      end if
      relative_error = abs_error / max(abs_reference, 1.0e-10)
      max_error = max(max_error, relative_error)
    end do

    write(*,*) 'Maximum relative error across all directions:', max_error
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in vector derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Vector derivatives are within tolerance (rtol + atol)'
    end if
  end subroutine check_derivatives_numerically

end program test_sasum_vector_forward