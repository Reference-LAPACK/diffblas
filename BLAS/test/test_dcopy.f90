! Test program for DCOPY differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_dcopy
  implicit none

  external :: dcopy
  external :: dcopy_d

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing DCOPY (multi-size: n = 4)'
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
    real(8), dimension(n) :: dy
    integer :: incy

    ! Derivative variables
    real(8), dimension(n) :: dy_d
    real(8), dimension(n) :: dx_d

    ! Array restoration and derivative storage
    real(8), dimension(n) :: dy_orig, dy_d_orig
    real(8), dimension(n) :: dx_orig, dx_d_orig
    integer :: i, j

    nsize = n
    incx = 1
    incy = 1

    call random_number(dx)
    dx = dx * 2.0d0 - 1.0d0  ! Scale to [-1,1]
    call random_number(dy)
    dy = dy * 2.0d0 - 1.0d0  ! Scale to [-1,1]

    ! Initialize input derivatives
    call random_number(dy_d)
    dy_d = dy_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
    call random_number(dx_d)
    dx_d = dx_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]

    ! Store _orig and _d_orig
    dy_d_orig = dy_d
    dx_d_orig = dx_d
    dy_orig = dy
    dx_orig = dx

    write(*,*) 'Testing DCOPY (n =', n, ')'

    ! Set ISIZE globals required by differentiated routine
    call set_ISIZE1OFDy(n)


    ! Call the differentiated function
    call dcopy_d(nsize, dx, dx_d, 1, dy, dy_d, 1)

    ! Reset ISIZE globals to uninitialized (-1)
    call set_ISIZE1OFDy(-1)

    write(*,*) 'Function calls completed successfully'

    ! Numerical differentiation check
    call check_derivatives_numerically(n, nsize, dx_orig, dy_orig, dx_d_orig, dy_d_orig, dy_d, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, nsize, dx_orig, dy_orig, dx_d_orig, dy_d_orig, dy_d, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nsize
    real(8), intent(in) :: dx_orig(n), dx_d_orig(n)
    real(8), intent(in) :: dy_orig(n), dy_d_orig(n)
    real(8), intent(in) :: dy_d(n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-6  ! Step size for finite differences
    real(8) :: relative_error, max_error
    real(8) :: abs_error, abs_reference, error_bound
    real(8) :: central_diff, ad_result
    logical :: has_large_errors
    real(8), dimension(n) :: dy_forward, dy_backward
    integer :: i, j
    real(8), dimension(n) :: dx
    real(8), dimension(n) :: dy

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    ! Forward perturbation: f(x + h)
    dx = dx_orig + h * dx_d_orig
    dy = dy_orig + h * dy_d_orig
    call dcopy(nsize, dx, 1, dy, 1)
    dy_forward = dy

    ! Backward perturbation: f(x - h)
    dx = dx_orig - h * dx_d_orig
    dy = dy_orig - h * dy_d_orig
    call dcopy(nsize, dx, 1, dy, 1)
    dy_backward = dy

    ! Compute central differences and compare with AD results
    do i = 1, n
        central_diff = (dy_forward(i) - dy_backward(i)) / (2.0e0 * h)
        ad_result = dy_d(i)
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 1.0e-5 + 1.0e-5 * abs_reference
        if (abs_error > error_bound) then
          has_large_errors = .true.
          relative_error = abs_error / max(abs_reference, 1.0e-10)
          write(*,*) 'Large error in output DY(', i, '):'
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
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_dcopy