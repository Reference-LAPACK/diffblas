! Test program for DDOT vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirs=n
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_ddot_vector_reverse
  implicit none

  real(8), external :: ddot
  external :: ddot_bv

  integer :: nbdirs
  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing DDOT (Vector Reverse, multi-size: n =', test_sizes(1), ')'
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    nbdirs = test_sizes(i)
    call run_test_for_size(n_test, passed, nbdirs)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) then
    write(*,*) 'PASS: Vector reverse mode - all sizes completed successfully'
  else
    write(*,*) 'FAIL: Vector reverse mode - one or more sizes had derivative errors'
  end if

contains

  subroutine run_test_for_size(n, passed, nbdirs)
    implicit none
    integer, intent(in) :: n
    logical, intent(out) :: passed
    integer, intent(in) :: nbdirs

    integer :: nsize, incx_val, incy_val
    real(8), dimension(n) :: x, y
    real(8), dimension(nbdirs,n) :: xb, yb
    real(8), dimension(nbdirs) :: result_b, result_b_seed
    real(8), dimension(n) :: x_orig, y_orig
    integer :: k, i
    real(4) :: temp_real, temp_imag

    nsize = n
    incx_val = 1
    incy_val = 1

    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    call random_number(y)
    y = y * 2.0d0 - 1.0d0

    x_orig = x
    y_orig = y

    do k = 1, nbdirs
      call random_number(temp_real)
      result_b(k) = temp_real * 2.0d0 - 1.0d0
    end do
    result_b_seed = result_b

    xb = 0.0d0
    yb = 0.0d0

    write(*,*) 'Testing DDOT (Vector Reverse, n =', n, ')'

    call set_ISIZE1OFDx(n)
    call set_ISIZE1OFDy(n)

    call ddot_bv(nsize, x, xb, incx_val, y, yb, incy_val, result_b, nbdirs)

    call set_ISIZE1OFDx(-1)
    call set_ISIZE1OFDy(-1)

    call check_vjp_numerically(n, nbdirs, nsize, incx_val, incy_val, x_orig, y_orig, result_b_seed, xb, yb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nbdirs, nsize, incx_val, incy_val, x_orig, y_orig, result_b_seed, xb, yb, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    integer, intent(in) :: nsize, incx_val, incy_val
    real(8), intent(in) :: x_orig(n), y_orig(n)
    real(8), intent(in) :: result_b_seed(nbdirs)
    real(8), intent(in) :: xb(nbdirs,n), yb(nbdirs,n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    real(8), dimension(n) :: x_dir, y_dir
    real(8) :: result_forward, result_backward, result_central_diff
    real(8), dimension(n) :: x, y
    integer :: i, k
    real(4) :: temp_real, temp_imag
    logical :: has_large_errors

    max_error = 0.0d0
    has_large_errors = .false.

    do k = 1, nbdirs
      call random_number(x_dir)
      x_dir = x_dir * 2.0d0 - 1.0d0
      call random_number(y_dir)
      y_dir = y_dir * 2.0d0 - 1.0d0
      x = x_orig + h * x_dir
      y = y_orig + h * y_dir
      result_forward = ddot(nsize, x, incx_val, y, incy_val)
      x = x_orig - h * x_dir
      y = y_orig - h * y_dir
      result_backward = ddot(nsize, x, incx_val, y, incy_val)
      result_central_diff = (result_forward - result_backward) / (2.0d0 * h)
      vjp_fd = result_b_seed(k) * result_central_diff
      vjp_ad = 0.0d0
      do i = 1, n
        vjp_ad = vjp_ad + x_dir(i) * xb(k,i)
        vjp_ad = vjp_ad + y_dir(i) * yb(k,i)
      end do
      abs_error = abs(vjp_fd - vjp_ad)
      abs_reference = abs(vjp_ad)
      error_bound = 1.0e-5 + 1.0e-5 * abs_reference
      if (abs_error > error_bound) has_large_errors = .true.
      if (abs_reference > 1.0e-10) then
        relative_error = abs_error / abs_reference
      else
        relative_error = abs_error
      end if
      if (relative_error > max_error) max_error = relative_error
    end do

    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: VJP errors outside tolerance'
    else
      write(*,*) 'PASS: VJP within tolerance'
    end if

  end subroutine check_vjp_numerically

end program test_ddot_vector_reverse