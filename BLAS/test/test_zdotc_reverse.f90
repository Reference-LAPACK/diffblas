! Test program for ZDOTC reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_zdotc_reverse
  implicit none

  complex(8), external :: zdotc
  external :: zdotc_b

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing ZDOTC (multi-size: n = 4)'
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
    complex(8), dimension(n) :: zx
    integer :: incx_val
    complex(8), dimension(n) :: zy
    integer :: incy_val
    complex(8), dimension(n) :: zxb
    complex(8), dimension(n) :: zyb
    complex(8) :: zdotcb, zdotcb_orig
    complex(8), dimension(n) :: zx_orig
    complex(8), dimension(n) :: zy_orig
    real(4) :: temp_re, temp_im
    integer :: i, j

    nsize = n
    incx_val = 1
    incy_val = 1

    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      zx(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      zy(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do

    zx_orig = zx
    zy_orig = zy


    call random_number(temp_re)
    call random_number(temp_im)
    zdotcb = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    zdotcb_orig = zdotcb

    zxb = 0.0
    zyb = 0.0

    write(*,*) 'Testing ZDOTC (n =', n, ')'

    call set_ISIZE1OFZx(n)
    call set_ISIZE1OFZy(n)

    call zdotc_b(nsize, zx, zxb, incx_val, zy, zyb, incy_val, zdotcb)

    call set_ISIZE1OFZx(-1)
    call set_ISIZE1OFZy(-1)

    call check_vjp_numerically(n, nsize, incx_val, incy_val, zx_orig, zy_orig, zxb, zyb, zdotcb_orig, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nsize, incx_val, incy_val, zx_orig, zy_orig, zxb, zyb, zdotcb_orig, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nsize
    integer, intent(in) :: incx_val
    integer, intent(in) :: incy_val
    complex(8), intent(in) :: zx_orig(n)
    complex(8), intent(in) :: zy_orig(n)
    complex(8), intent(in) :: zxb(n)
    complex(8), intent(in) :: zyb(n)
    complex(8), intent(in) :: zdotcb_orig
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j, n_products
    real(8), dimension(n) :: temp_products
    real(4) :: temp_re, temp_im

    complex(8), dimension(n) :: zx_dir
    complex(8), dimension(n) :: zy_dir

    complex(8) :: zdotc_plus, zdotc_minus

    complex(8), dimension(n) :: zx
    complex(8), dimension(n) :: zy

    max_error = 0.0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      zx_dir(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      zy_dir(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do

    zx = zx_orig + cmplx(h, 0.0) * zx_dir
    zy = zy_orig + cmplx(h, 0.0) * zy_dir
    zdotc_plus = zdotc(nsize, zx, incx_val, zy, incy_val)

    zx = zx_orig - cmplx(h, 0.0) * zx_dir
    zy = zy_orig - cmplx(h, 0.0) * zy_dir
    zdotc_minus = zdotc(nsize, zx, incx_val, zy, incy_val)


    vjp_fd = real(conjg(zdotcb_orig) * (zdotc_plus - zdotc_minus) / (2.0 * h))

    vjp_ad = 0.0
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(zx_dir(i)) * zxb(i))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(zy_dir(i)) * zyb(i))
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

end program test_zdotc_reverse