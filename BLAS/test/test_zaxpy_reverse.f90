! Test program for ZAXPY reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_zaxpy_reverse
  implicit none

  external :: zaxpy
  external :: zaxpy_b

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(3)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing ZAXPY (multi-size: n = 4)'
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
    complex(8) :: za
    complex(8), dimension(n) :: zx
    integer :: incx_val
    complex(8), dimension(n) :: zy
    integer :: incy_val
    complex(8) :: zab
    complex(8), dimension(n) :: zxb
    complex(8), dimension(n) :: zyb
    complex(8) :: za_orig
    complex(8), dimension(n) :: zx_orig
    complex(8), dimension(n) :: zy_orig
    complex(8), dimension(n) :: zyb_orig
    real(4) :: temp_re, temp_im
    integer :: i, j

    nsize = n
    incx_val = 1
    incy_val = 1

    call random_number(temp_re)
    call random_number(temp_im)
    za = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
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

    za_orig = za
    zx_orig = zx
    zy_orig = zy

    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      zyb(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do
    zyb_orig = zyb

    zab = 0.0
    zxb = 0.0

    write(*,*) 'Testing ZAXPY (n =', n, ')'

    call set_ISIZE1OFZx(n)

    call zaxpy_b(nsize, za, zab, zx, zxb, incx_val, zy, zyb, incy_val)

    call set_ISIZE1OFZx(-1)

    call check_vjp_numerically(n, nsize, incx_val, incy_val, za_orig, zx_orig, zy_orig, zyb_orig, zab, zxb, zyb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nsize, incx_val, incy_val, za_orig, zx_orig, zy_orig, zyb_orig, zab, zxb, zyb, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nsize
    integer, intent(in) :: incx_val
    integer, intent(in) :: incy_val
    complex(8), intent(in) :: za_orig
    complex(8), intent(in) :: zx_orig(n)
    complex(8), intent(in) :: zy_orig(n)
    complex(8), intent(in) :: zyb_orig(n)
    complex(8), intent(in) :: zab
    complex(8), intent(in) :: zxb(n)
    complex(8), intent(in) :: zyb(n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j, n_products
    real(8), dimension(n) :: temp_products
    real(4) :: temp_re, temp_im

    complex(8) :: za_dir
    complex(8), dimension(n) :: zx_dir
    complex(8), dimension(n) :: zy_dir

    complex(8), dimension(n) :: zy_plus, zy_minus, zy_central_diff

    complex(8) :: za
    complex(8), dimension(n) :: zx
    complex(8), dimension(n) :: zy

    max_error = 0.0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    call random_number(temp_re)
    call random_number(temp_im)
    za_dir = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
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

    za = za_orig + cmplx(h, 0.0) * za_dir
    zx = zx_orig + cmplx(h, 0.0) * zx_dir
    zy = zy_orig + cmplx(h, 0.0) * zy_dir
    call zaxpy(nsize, za, zx, incx_val, zy, incy_val)
    zy_plus = zy

    za = za_orig - cmplx(h, 0.0) * za_dir
    zx = zx_orig - cmplx(h, 0.0) * zx_dir
    zy = zy_orig - cmplx(h, 0.0) * zy_dir
    call zaxpy(nsize, za, zx, incx_val, zy, incy_val)
    zy_minus = zy

    zy_central_diff = (zy_plus - zy_minus) / (2.0 * h)

    vjp_fd = 0.0
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(zyb_orig(i)) * zy_central_diff(i))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do

    vjp_ad = 0.0
    vjp_ad = vjp_ad + real(conjg(za_dir) * zab)
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

end program test_zaxpy_reverse