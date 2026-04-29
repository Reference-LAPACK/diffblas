! Test program for SSCAL reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_sscal_reverse
  implicit none

  external :: sscal
  external :: sscal_b

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(3)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing SSCAL (multi-size: n = 4)'
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
    real(4) :: sa
    real(4), dimension(n) :: sx
    integer :: incx_val
    real(4) :: sab
    real(4), dimension(n) :: sxb
    real(4) :: sa_orig
    real(4), dimension(n) :: sx_orig
    real(4), dimension(n) :: sxb_orig
    integer :: i, j

    nsize = n
    incx_val = 1

    call random_number(sa)
    sa = sa * 2.0 - 1.0
    call random_number(sx)
    sx = sx * 2.0 - 1.0

    sa_orig = sa
    sx_orig = sx

    call random_number(sxb)
    sxb = sxb * 2.0 - 1.0
    sxb_orig = sxb

    sab = 0.0

    write(*,*) 'Testing SSCAL (n =', n, ')'

    call sscal_b(nsize, sa, sab, sx, sxb, incx_val)

    call check_vjp_numerically(n, nsize, incx_val, sa_orig, sx_orig, sxb_orig, sab, sxb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nsize, incx_val, sa_orig, sx_orig, sxb_orig, sab, sxb, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nsize
    integer, intent(in) :: incx_val
    real(4), intent(in) :: sa_orig
    real(4), intent(in) :: sx_orig(n)
    real(4), intent(in) :: sxb_orig(n)
    real(4), intent(in) :: sab
    real(4), intent(in) :: sxb(n)
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j, n_products
    real(4), dimension(n) :: temp_products

    real(4) :: sa_dir
    real(4), dimension(n) :: sx_dir

    real(4), dimension(n) :: sx_plus, sx_minus, sx_central_diff

    real(4) :: sa
    real(4), dimension(n) :: sx

    max_error = 0.0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    call random_number(sa_dir)
    sa_dir = sa_dir * 2.0 - 1.0
    call random_number(sx_dir)
    sx_dir = sx_dir * 2.0 - 1.0

    sa = sa_orig + h * sa_dir
    sx = sx_orig + h * sx_dir
    call sscal(nsize, sa, sx, incx_val)
    sx_plus = sx

    sa = sa_orig - h * sa_dir
    sx = sx_orig - h * sx_dir
    call sscal(nsize, sa, sx, incx_val)
    sx_minus = sx

    sx_central_diff = (sx_plus - sx_minus) / (2.0 * h)

    vjp_fd = 0.0
    n_products = n
    do i = 1, n
      temp_products(i) = sxb_orig(i) * sx_central_diff(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do

    vjp_ad = 0.0
    vjp_ad = vjp_ad + sa_dir * sab
    n_products = n
    do i = 1, n
      temp_products(i) = sx_dir(i) * sxb(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do

    abs_error = abs(vjp_fd - vjp_ad)
    abs_reference = abs(vjp_ad)
    error_bound = 2.0e-3 + 2.0e-3 * abs_reference
    if (abs_error > error_bound) has_large_errors = .true.
    if (abs_reference > 1.0e-10) then
      relative_error = abs_error / abs_reference
    else
      relative_error = abs_error
    end if
    max_error = relative_error
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
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
    real(4), dimension(n), intent(inout) :: arr
    integer :: i, j, min_idx
    real(4) :: temp
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

end program test_sscal_reverse