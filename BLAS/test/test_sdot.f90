! Test program for SDOT differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_sdot
  implicit none

  real(4), external :: sdot
  real(4), external :: sdot_d

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing SDOT (multi-size: n = 4)'
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
    real(4), dimension(n) :: sx
    integer :: incx
    real(4), dimension(n) :: sy
    integer :: incy

    ! Derivative variables
    real(4), dimension(n) :: sx_d
    real(4), dimension(n) :: sy_d
    real(4) :: sdot_d_result  ! Derivative of function result (avoid name clash with func_d)

    ! Array restoration and derivative storage
    real(4), dimension(n) :: sx_orig, sx_d_orig
    real(4), dimension(n) :: sy_orig, sy_d_orig
    real(4) :: sdot_orig  ! Function result (no _d_orig - use _d_result)
    integer :: i, j

    nsize = n
    incx = 1
    incy = 1

    call random_number(sx)
    sx = sx * 2.0d0 - 1.0d0  ! Scale to [-1,1]
    call random_number(sy)
    sy = sy * 2.0d0 - 1.0d0  ! Scale to [-1,1]

    ! Initialize input derivatives
    call random_number(sx_d)
    sx_d = sx_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
    call random_number(sy_d)
    sy_d = sy_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]

    ! Store _orig and _d_orig
    sx_d_orig = sx_d
    sy_d_orig = sy_d
    sx_orig = sx
    sy_orig = sy
    sdot_orig = sdot(nsize, sx, 1, sy, 1)

    write(*,*) 'Testing SDOT (n =', n, ')'

    ! Call the differentiated function
    sdot_d_result = sdot_d(nsize, sx, sx_d, 1, sy, sy_d, 1, sdot_orig)

    write(*,*) 'Function calls completed successfully'

    ! Numerical differentiation check
    call check_derivatives_numerically(n, nsize, sx_orig, sy_orig, sdot_orig, sx_d_orig, sy_d_orig, sdot_d_result, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, nsize, sx_orig, sy_orig, sdot_orig, sx_d_orig, sy_d_orig, sdot_d_result, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nsize
    real(4), intent(in) :: sx_orig(n), sx_d_orig(n)
    real(4), intent(in) :: sy_orig(n), sy_d_orig(n)
    real(4), intent(in) :: sdot_orig
    real(4), intent(in) :: sdot_d_result
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3  ! Step size for finite differences
    real(4) :: relative_error, max_error
    real(4) :: abs_error, abs_reference, error_bound
    real(4) :: central_diff, ad_result
    logical :: has_large_errors
    real(4) :: sdot_forward, sdot_backward  ! Function result for FD check
    integer :: i, j
    real(4), dimension(n) :: sx
    real(4), dimension(n) :: sy

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    ! Forward perturbation: f(x + h)
    sx = sx_orig + h * sx_d_orig
    sy = sy_orig + h * sy_d_orig
    sdot_forward = sdot(nsize, sx, 1, sy, 1)

    ! Backward perturbation: f(x - h)
    sx = sx_orig - h * sx_d_orig
    sy = sy_orig - h * sy_d_orig
    sdot_backward = sdot(nsize, sx, 1, sy, 1)

    ! Compute central differences and compare with AD results
    central_diff = (sdot_forward - sdot_backward) / (2.0e0 * h)
    ad_result = sdot_d_result
    abs_error = abs(central_diff - ad_result)
    abs_reference = abs(ad_result)
    error_bound = 2.0e-3 + 2.0e-3 * abs_reference
    if (abs_error > error_bound) then
      has_large_errors = .true.
      relative_error = abs_error / max(abs_reference, 1.0e-10)
      write(*,*) 'Large error in function result SDOT:'
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
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_sdot