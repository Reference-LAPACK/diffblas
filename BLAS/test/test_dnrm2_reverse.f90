! Test program for DNRM2 reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_dnrm2_reverse
  implicit none

  real(8), external :: dnrm2
  external :: dnrm2_b

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
    integer :: incx_val
    real(8), dimension(n) :: xb
    real(8) :: dnrm2b, dnrm2b_orig
    real(8), dimension(n) :: x_orig
    integer :: i, j

    nsize = n
    incx_val = 1

    call random_number(x)
    x = x * 2.0 - 1.0

    x_orig = x


    call random_number(dnrm2b)
    dnrm2b = dnrm2b * 2.0 - 1.0
    dnrm2b_orig = dnrm2b

    xb = 0.0

    write(*,*) 'Testing DNRM2 (n =', n, ')'

    call dnrm2_b(nsize, x, xb, incx_val, dnrm2b)

    call check_vjp_numerically(n, nsize, incx_val, x_orig, xb, dnrm2b_orig, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nsize, incx_val, x_orig, xb, dnrm2b_orig, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nsize
    integer, intent(in) :: incx_val
    real(8), intent(in) :: x_orig(n)
    real(8), intent(in) :: xb(n)
    real(8), intent(in) :: dnrm2b_orig
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j, n_products
    real(8), dimension(n) :: temp_products

    real(8), dimension(n) :: x_dir

    real(8) :: dnrm2_plus, dnrm2_minus

    real(8), dimension(n) :: x

    max_error = 0.0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    call random_number(x_dir)
    x_dir = x_dir * 2.0 - 1.0

    x = x_orig + h * x_dir
    dnrm2_plus = dnrm2(nsize, x, incx_val)

    x = x_orig - h * x_dir
    dnrm2_minus = dnrm2(nsize, x, incx_val)


    vjp_fd = dnrm2b_orig * (dnrm2_plus - dnrm2_minus) / (2.0 * h)

    vjp_ad = 0.0
    n_products = n
    do i = 1, n
      temp_products(i) = x_dir(i) * xb(i)
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

end program test_dnrm2_reverse