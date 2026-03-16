! Test program for SNRM2 differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_snrm2
  implicit none

  real(4), external :: snrm2
  real(4), external :: snrm2_d

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(3)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing SNRM2 (multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 3
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
    real(4), dimension(n) :: x
    integer :: incx

    ! Derivative variables
    real(4) :: snrm2_d_result  ! Derivative of function result (avoid name clash with func_d)
    real(4), dimension(n) :: x_d

    ! Array restoration and derivative storage
    real(4) :: snrm2_orig  ! Function result (no _d_orig - use _d_result)
    real(4), dimension(n) :: x_orig, x_d_orig
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
    snrm2_orig = snrm2(nsize, x, 1)
    x_orig = x

    write(*,*) 'Testing SNRM2 (n =', n, ')'

    ! Call the differentiated function
    snrm2_d_result = snrm2_d(nsize, x, x_d, 1, snrm2_orig)
    x_d = x_d_orig

    write(*,*) 'Function calls completed successfully'

    ! Numerical differentiation check
    call check_derivatives_numerically(n, nsize, x_orig, snrm2_orig, x_d_orig, snrm2_d_result, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, nsize, x_orig, snrm2_orig, x_d_orig, snrm2_d_result, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nsize
    real(4), intent(in) :: x_orig(n), x_d_orig(n)
    real(4), intent(in) :: snrm2_orig
    real(4), intent(in) :: snrm2_d_result
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3  ! Step size for finite differences
    real(4) :: relative_error, max_error
    real(4) :: abs_error, abs_reference, error_bound
    real(4) :: central_diff, ad_result
    logical :: has_large_errors
    real(4) :: snrm2_forward, snrm2_backward  ! Function result for FD check
    integer :: i, j
    real(4), dimension(n) :: x

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    ! Forward perturbation: f(x + h)
    x = x_orig + h * x_d_orig
    snrm2_forward = snrm2(nsize, x, 1)

    ! Backward perturbation: f(x - h)
    x = x_orig - h * x_d_orig
    snrm2_backward = snrm2(nsize, x, 1)

    ! Compute central differences and compare with AD results
    central_diff = (snrm2_forward - snrm2_backward) / (2.0e0 * h)
    ad_result = snrm2_d_result
    abs_error = abs(central_diff - ad_result)
    abs_reference = abs(ad_result)
    error_bound = 2.0e-3 + 2.0e-3 * abs_reference
    if (abs_error > error_bound) then
      has_large_errors = .true.
      relative_error = abs_error / max(abs_reference, 1.0e-10)
      write(*,*) 'Large error in function result SNRM2:'
      write(*,*) '  Central diff: ', central_diff
      write(*,*) '  AD result:   ', ad_result
      write(*,*) '  Absolute error:', abs_error
      write(*,*) '  Error bound:', error_bound
      write(*,*) '  Relative error:', relative_error
    end if
    relative_error = abs_error / max(abs_reference, 1.0e-10)
    max_error = max(max_error, relative_error)

    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_snrm2