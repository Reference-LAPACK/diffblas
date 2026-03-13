! Test program for DNRM2 differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_dnrm2
  implicit none

  real(8), external :: dnrm2
  real(8), external :: dnrm2_d

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing DNRM2 (multi-size: n = 4)'
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
    real(8), dimension(n) :: x
    integer :: incx

    ! Derivative variables
    real(8) :: dnrm2_d_result  ! Derivative of function result (avoid name clash with func_d)
    real(8), dimension(n) :: x_d

    ! Array restoration and derivative storage
    real(8) :: dnrm2_orig  ! Function result (no _d_orig - use _d_result)
    real(8), dimension(n) :: x_orig, x_d_orig
    integer :: i, j

    nsize = n
    incx = 1

    call random_number(x)
    x = x * 2.0d0 - 1.0d0  ! Scale to [-1,1]

    ! Initialize input derivatives
    call random_number(x_d)
    x_d = x_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]

    ! Store _orig and _d_orig
    x_d_orig = x_d
    dnrm2_orig = dnrm2(nsize, x, 1)
    x_orig = x

    write(*,*) 'Testing DNRM2 (n =', n, ')'

    ! Call the differentiated function
    dnrm2_d_result = dnrm2_d(nsize, x, x_d, 1, dnrm2_orig)

    write(*,*) 'Function calls completed successfully'

    ! Numerical differentiation check
    call check_derivatives_numerically(n, nsize, x_orig, dnrm2_orig, x_d_orig, dnrm2_d_result, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, nsize, x_orig, dnrm2_orig, x_d_orig, dnrm2_d_result, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nsize
    real(8), intent(in) :: x_orig(n), x_d_orig(n)
    real(8), intent(in) :: dnrm2_orig
    real(8), intent(in) :: dnrm2_d_result
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-6  ! Step size for finite differences
    real(8) :: relative_error, max_error
    real(8) :: abs_error, abs_reference, error_bound
    real(8) :: central_diff, ad_result
    logical :: has_large_errors
    real(8) :: dnrm2_forward, dnrm2_backward  ! Function result for FD check
    integer :: i, j
    real(8), dimension(n) :: x

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    ! Forward perturbation: f(x + h)
    x = x_orig + h * x_d_orig
    dnrm2_forward = dnrm2(nsize, x, 1)

    ! Backward perturbation: f(x - h)
    x = x_orig - h * x_d_orig
    dnrm2_backward = dnrm2(nsize, x, 1)

    ! Compute central differences and compare with AD results
    central_diff = (dnrm2_forward - dnrm2_backward) / (2.0e0 * h)
    ad_result = dnrm2_d_result
    abs_error = abs(central_diff - ad_result)
    abs_reference = abs(ad_result)
    error_bound = 1.0e-5 + 1.0e-5 * abs_reference
    if (abs_error > error_bound) then
      has_large_errors = .true.
      relative_error = abs_error / max(abs_reference, 1.0e-10)
      write(*,*) 'Large error in function result DNRM2:'
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
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_dnrm2