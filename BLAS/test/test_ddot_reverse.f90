! Test program for DDOT reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_ddot_reverse
  implicit none

  real(8), external :: ddot
  external :: ddot_b

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing DDOT (multi-size: n = 4)'
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
    integer :: incx_val
    real(8), dimension(n) :: dy
    integer :: incy_val
    real(8), dimension(n) :: dxb
    real(8), dimension(n) :: dyb
    real(8) :: ddotb, ddotb_orig
    real(8), dimension(n) :: dx_orig
    real(8), dimension(n) :: dy_orig
    integer :: i, j

    nsize = n
    incx_val = 1
    incy_val = 1

    call random_number(dx)
    dx = dx * 2.0 - 1.0
    call random_number(dy)
    dy = dy * 2.0 - 1.0

    dx_orig = dx
    dy_orig = dy


    call random_number(ddotb)
    ddotb = ddotb * 2.0 - 1.0
    ddotb_orig = ddotb

    dxb = 0.0
    dyb = 0.0

    write(*,*) 'Testing DDOT (n =', n, ')'

    call set_ISIZE1OFDx(n)
    call set_ISIZE1OFDy(n)

    call ddot_b(nsize, dx, dxb, incx_val, dy, dyb, incy_val, ddotb)

    call set_ISIZE1OFDx(-1)
    call set_ISIZE1OFDy(-1)

    call check_vjp_numerically(n, nsize, incx_val, incy_val, dx_orig, dy_orig, dxb, dyb, ddotb_orig, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nsize, incx_val, incy_val, dx_orig, dy_orig, dxb, dyb, ddotb_orig, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nsize
    integer, intent(in) :: incx_val
    integer, intent(in) :: incy_val
    real(8), intent(in) :: dx_orig(n)
    real(8), intent(in) :: dy_orig(n)
    real(8), intent(in) :: dxb(n)
    real(8), intent(in) :: dyb(n)
    real(8), intent(in) :: ddotb_orig
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j, n_products
    real(8), dimension(n) :: temp_products

    real(8), dimension(n) :: dx_dir
    real(8), dimension(n) :: dy_dir

    real(8) :: ddot_plus, ddot_minus

    real(8), dimension(n) :: dx
    real(8), dimension(n) :: dy

    max_error = 0.0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    call random_number(dx_dir)
    dx_dir = dx_dir * 2.0 - 1.0
    call random_number(dy_dir)
    dy_dir = dy_dir * 2.0 - 1.0

    dx = dx_orig + h * dx_dir
    dy = dy_orig + h * dy_dir
    ddot_plus = ddot(nsize, dx, incx_val, dy, incy_val)

    dx = dx_orig - h * dx_dir
    dy = dy_orig - h * dy_dir
    ddot_minus = ddot(nsize, dx, incx_val, dy, incy_val)


    vjp_fd = ddotb_orig * (ddot_plus - ddot_minus) / (2.0 * h)

    vjp_ad = 0.0
    n_products = n
    do i = 1, n
      temp_products(i) = dx_dir(i) * dxb(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    n_products = n
    do i = 1, n
      temp_products(i) = dy_dir(i) * dyb(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
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
    max_error = relative_error

    write(*,*) ''
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_vjp_numerically

  subroutine sort_array(arr, n)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(inout) :: arr
    integer :: i, j, min_idx
    real(8) :: temp
    do i = 1, n-1
      min_idx = i
      do j = i+1, n
        if (abs(arr(j)) < abs(arr(min_idx))) min_idx = j
      end do
      if (min_idx /= i) then
        temp = arr(i)
        arr(i) = arr(min_idx)
        arr(min_idx) = temp
      end if
    end do
  end subroutine sort_array

end program test_ddot_reverse