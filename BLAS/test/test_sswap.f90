! Test program for SSWAP differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_sswap
  implicit none

  external :: sswap
  external :: sswap_d

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(3)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing SSWAP (multi-size: n = 4)'
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
    real(4), dimension(n) :: sx
    integer :: incx
    real(4), dimension(n) :: sy
    integer :: incy

    ! Derivative variables
    real(4), dimension(n) :: sx_d
    real(4), dimension(n) :: sy_d

    ! Array restoration and derivative storage
    real(4), dimension(n) :: sx_orig, sx_d_orig
    real(4), dimension(n) :: sy_orig, sy_d_orig
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

    write(*,*) 'Testing SSWAP (n =', n, ')'
    sx_orig = sx
    sy_orig = sy

    ! Call the differentiated function
    call sswap_d(nsize, sx, sx_d, 1, sy, sy_d, 1)

    write(*,*) 'Function calls completed successfully'

    ! Numerical differentiation check
    call check_derivatives_numerically(n, nsize, sx_orig, sy_orig, sx_d_orig, sy_d_orig, sx_d, sy_d, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, nsize, sx_orig, sy_orig, sx_d_orig, sy_d_orig, sx_d, sy_d, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nsize
    real(4), intent(in) :: sx_orig(n), sx_d_orig(n)
    real(4), intent(in) :: sy_orig(n), sy_d_orig(n)
    real(4), intent(in) :: sx_d(n)
    real(4), intent(in) :: sy_d(n)
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3  ! Step size for finite differences
    real(4) :: relative_error, max_error
    real(4) :: abs_error, abs_reference, error_bound
    real(4) :: central_diff, ad_result
    logical :: has_large_errors
    real(4), dimension(n) :: sx_forward, sx_backward
    real(4), dimension(n) :: sy_forward, sy_backward
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
    call sswap(nsize, sx, 1, sy, 1)
    sx_forward = sx
    sy_forward = sy

    ! Backward perturbation: f(x - h)
    sx = sx_orig - h * sx_d_orig
    sy = sy_orig - h * sy_d_orig
    call sswap(nsize, sx, 1, sy, 1)
    sx_backward = sx
    sy_backward = sy

    ! Compute central differences and compare with AD results
    do i = 1, n
        central_diff = (sx_forward(i) - sx_backward(i)) / (2.0e0 * h)
        ad_result = sx_d(i)
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 2.0e-3 + 2.0e-3 * abs_reference
        if (abs_error > error_bound) then
          has_large_errors = .true.
          relative_error = abs_error / max(abs_reference, 1.0e-10)
          write(*,*) 'Large error in output SX(', i, '):'
          write(*,*) '  Central diff: ', central_diff
          write(*,*) '  AD result:   ', ad_result
          write(*,*) '  Absolute error:', abs_error
          write(*,*) '  Error bound:', error_bound
          write(*,*) '  Relative error:', relative_error
        end if
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        max_error = max(max_error, relative_error)
    end do
    do i = 1, n
        central_diff = (sy_forward(i) - sy_backward(i)) / (2.0e0 * h)
        ad_result = sy_d(i)
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 2.0e-3 + 2.0e-3 * abs_reference
        if (abs_error > error_bound) then
          has_large_errors = .true.
          relative_error = abs_error / max(abs_reference, 1.0e-10)
          write(*,*) 'Large error in output SY(', i, '):'
          write(*,*) '  Central diff: ', central_diff
          write(*,*) '  AD result:   ', ad_result
          write(*,*) '  Absolute error:', abs_error
          write(*,*) '  Error bound:', error_bound
          write(*,*) '  Relative error:', relative_error
        end if
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        max_error = max(max_error, relative_error)
    end do

    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_sswap