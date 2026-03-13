! Test program for DSWAP reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_dswap_reverse
  implicit none

  external :: dswap
  external :: dswap_b

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing DSWAP (multi-size: n = 4)'
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
    real(8), dimension(n) :: dx_orig
    real(8), dimension(n) :: dy_orig
    real(8), dimension(n) :: dxb_orig
    real(8), dimension(n) :: dyb_orig
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

    call random_number(dxb)
    dxb = dxb * 2.0 - 1.0
    call random_number(dyb)
    dyb = dyb * 2.0 - 1.0
    dxb_orig = dxb
    dyb_orig = dyb


    write(*,*) 'Testing DSWAP (n =', n, ')'

    call dswap_b(nsize, dx, dxb, incx_val, dy, dyb, incy_val)

    call check_vjp_numerically(n, nsize, incx_val, incy_val, dx_orig, dy_orig, dxb_orig, dyb_orig, dxb, dyb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nsize, incx_val, incy_val, dx_orig, dy_orig, dxb_orig, dyb_orig, dxb, dyb, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nsize
    integer, intent(in) :: incx_val
    integer, intent(in) :: incy_val
    real(8), intent(in) :: dx_orig(n)
    real(8), intent(in) :: dy_orig(n)
    real(8), intent(in) :: dxb_orig(n)
    real(8), intent(in) :: dyb_orig(n)
    real(8), intent(in) :: dxb(n)
    real(8), intent(in) :: dyb(n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j, n_products
    real(8), dimension(n) :: temp_products

    real(8), dimension(n) :: dx_dir
    real(8), dimension(n) :: dy_dir

    real(8), dimension(n) :: dy_plus, dy_minus, dy_central_diff
    real(8), dimension(n) :: dx_plus, dx_minus, dx_central_diff

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
    call dswap(nsize, dx, incx_val, dy, incy_val)
    dy_plus = dy
    dx_plus = dx

    dx = dx_orig - h * dx_dir
    dy = dy_orig - h * dy_dir
    call dswap(nsize, dx, incx_val, dy, incy_val)
    dy_minus = dy
    dx_minus = dx

    dy_central_diff = (dy_plus - dy_minus) / (2.0 * h)
    dx_central_diff = (dx_plus - dx_minus) / (2.0 * h)

    vjp_fd = 0.0
    n_products = n
    do i = 1, n
      temp_products(i) = dyb_orig(i) * dy_central_diff(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do
    n_products = n
    do i = 1, n
      temp_products(i) = dxb_orig(i) * dx_central_diff(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do

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
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
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

end program test_dswap_reverse