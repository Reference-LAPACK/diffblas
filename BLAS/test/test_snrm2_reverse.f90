! Test program for SNRM2 reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_snrm2_reverse
  implicit none

  real(4), external :: snrm2
  external :: snrm2_b

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing SNRM2 (multi-size: n = 4)'
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
    real(4), dimension(n) :: x
    integer :: incx_val
    real(4), dimension(n) :: xb
    real(4) :: snrm2b, snrm2b_orig
    real(4), dimension(n) :: x_orig
    integer :: i, j

    nsize = n
    incx_val = 1

    call random_number(x)
    x = x * 2.0 - 1.0

    x_orig = x


    call random_number(snrm2b)
    snrm2b = snrm2b * 2.0 - 1.0
    snrm2b_orig = snrm2b

    xb = 0.0

    write(*,*) 'Testing SNRM2 (n =', n, ')'

    call snrm2_b(nsize, x, xb, incx_val, snrm2b)

    call check_vjp_numerically(n, nsize, incx_val, x_orig, xb, snrm2b_orig, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nsize, incx_val, x_orig, xb, snrm2b_orig, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nsize
    integer, intent(in) :: incx_val
    real(4), intent(in) :: x_orig(n)
    real(4), intent(in) :: xb(n)
    real(4), intent(in) :: snrm2b_orig
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j, n_products
    real(4), dimension(n) :: temp_products

    real(4), dimension(n) :: x_dir

    real(4) :: snrm2_plus, snrm2_minus

    real(4), dimension(n) :: x

    max_error = 0.0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    call random_number(x_dir)
    x_dir = x_dir * 2.0 - 1.0

    x = x_orig + h * x_dir
    snrm2_plus = snrm2(nsize, x, incx_val)

    x = x_orig - h * x_dir
    snrm2_minus = snrm2(nsize, x, incx_val)


    vjp_fd = snrm2b_orig * (snrm2_plus - snrm2_minus) / (2.0 * h)

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
    error_bound = 1.0e-3 + 1.0e-3 * abs_reference
    if (abs_error > error_bound) has_large_errors = .true.
    if (abs_reference > 1.0e-10) then
      relative_error = abs_error / abs_reference
    else
      relative_error = abs_error
    end if
    max_error = relative_error
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
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

end program test_snrm2_reverse