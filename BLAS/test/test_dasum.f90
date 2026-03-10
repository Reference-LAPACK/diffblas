! Test program for DASUM differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_dasum
  implicit none

  real(8), external :: dasum
  real(8), external :: dasum_d

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing DASUM (multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    call run_test_for_size(n_test, passed)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) then
    write(*,*) 'PASS: All sizes completed successfully'
  else
    write(*,*) 'FAIL: One or more sizes had derivative errors'
  end if

contains

  subroutine run_test_for_size(n, passed)
    implicit none
    integer, intent(in) :: n
    logical, intent(out) :: passed

    integer :: nsize
    real(8), dimension(n) :: dx
    integer :: incx

    ! Derivative variables
    real(8), dimension(n) :: dx_d
    real(8) :: dasum_d_result  ! Derivative of function result (avoid name clash with func_d)

    ! Array restoration and derivative storage
    real(8), dimension(n) :: dx_orig, dx_d_orig
    real(8) :: dasum_orig  ! Function result (no _d_orig - use _d_result)
    integer :: i, j

    nsize = n
    incx = 1

    call random_number(dx)
    dx = dx * 2.0d0 - 1.0d0  ! Scale to [-1,1]

    ! Initialize input derivatives
    call random_number(dx_d)
    dx_d = dx_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]

    ! Store _orig and _d_orig
    dx_d_orig = dx_d
    dx_orig = dx
    dasum_orig = dasum(nsize, dx, 1)

    write(*,*) 'Testing DASUM (n =', n, ')'

    ! Call the differentiated function
    dasum_d_result = dasum_d(nsize, dx, dx_d, 1, dasum_orig)

    write(*,*) 'Function calls completed successfully'

    ! Numerical differentiation check
    call check_derivatives_numerically(n, nsize, dx_orig, dasum_orig, dx_d_orig, dasum_d_result, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, nsize, dx_orig, dasum_orig, dx_d_orig, dasum_d_result, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nsize
    real(8), intent(in) :: dx_orig(n), dx_d_orig(n)
    real(8), intent(in) :: dasum_orig
    real(8), intent(in) :: dasum_d_result
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-6  ! Step size for finite differences
    real(8) :: relative_error, max_error
    real(8) :: abs_error, abs_reference, error_bound
    real(8) :: central_diff, ad_result
    logical :: has_large_errors
    real(8) :: dasum_forward, dasum_backward  ! Function result for FD check
    integer :: i, j
    real(8), dimension(n) :: dx

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    ! Forward perturbation: f(x + h)
    dx = dx_orig + h * dx_d_orig
    dasum_forward = dasum(nsize, dx, 1)

    ! Backward perturbation: f(x - h)
    dx = dx_orig - h * dx_d_orig
    dasum_backward = dasum(nsize, dx, 1)

    ! Compute central differences and compare with AD results
    central_diff = (dasum_forward - dasum_backward) / (2.0e0 * h)
    ad_result = dasum_d_result
    abs_error = abs(central_diff - ad_result)
    abs_reference = abs(ad_result)
    error_bound = 1.0e-5 + 1.0e-5 * abs_reference
    if (abs_error > error_bound) then
      has_large_errors = .true.
      relative_error = abs_error / max(abs_reference, 1.0e-10)
      write(*,*) 'Large error in function result DASUM:'
      write(*,*) '  Central diff: ', central_diff
      write(*,*) '  AD result:   ', ad_result
      write(*,*) '  Absolute error:', abs_error
      write(*,*) '  Error bound:', error_bound
      write(*,*) '  Relative error:', relative_error
    end if
    relative_error = abs_error / max(abs_reference, 1.0e-10)
    max_error = max(max_error, relative_error)

    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_dasum