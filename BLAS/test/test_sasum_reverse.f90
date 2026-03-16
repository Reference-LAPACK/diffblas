! Test program for SASUM reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_sasum_reverse
  implicit none

  real(4), external :: sasum
  external :: sasum_b

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(3)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing SASUM (multi-size: n = 4)'
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
    integer :: incx_val
    real(4), dimension(n) :: sxb
    real(4) :: sasumb, sasumb_orig
    real(4), dimension(n) :: sx_orig
    integer :: i, j

    nsize = n
    incx_val = 1

    call random_number(sx)
    sx = sx * 2.0 - 1.0

    sx_orig = sx


    call random_number(sasumb)
    sasumb = sasumb * 2.0 - 1.0
    sasumb_orig = sasumb

    sxb = 0.0

    write(*,*) 'Testing SASUM (n =', n, ')'

    call set_ISIZE1OFSx(n)

    call sasum_b(nsize, sx, sxb, incx_val, sasumb)

    call set_ISIZE1OFSx(-1)

    call check_vjp_numerically(n, nsize, incx_val, sx_orig, sxb, sasumb_orig, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nsize, incx_val, sx_orig, sxb, sasumb_orig, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nsize
    integer, intent(in) :: incx_val
    real(4), intent(in) :: sx_orig(n)
    real(4), intent(in) :: sxb(n)
    real(4), intent(in) :: sasumb_orig
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j, n_products
    real(4), dimension(n) :: temp_products

    real(4), dimension(n) :: sx_dir

    real(4) :: sasum_plus, sasum_minus

    real(4), dimension(n) :: sx

    max_error = 0.0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    call random_number(sx_dir)
    sx_dir = sx_dir * 2.0 - 1.0

    sx = sx_orig + h * sx_dir
    sasum_plus = sasum(nsize, sx, incx_val)

    sx = sx_orig - h * sx_dir
    sasum_minus = sasum(nsize, sx, incx_val)


    vjp_fd = sasumb_orig * (sasum_plus - sasum_minus) / (2.0 * h)

    vjp_ad = 0.0
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

end program test_sasum_reverse