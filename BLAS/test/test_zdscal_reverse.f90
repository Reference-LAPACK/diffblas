! Test program for ZDSCAL reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_zdscal_reverse
  implicit none

  external :: zdscal
  external :: zdscal_b

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing ZDSCAL (multi-size: n = 4)'
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
    real(8) :: da
    complex(8), dimension(n) :: zx
    integer :: incx_val
    real(8) :: dab
    complex(8), dimension(n) :: zxb
    real(8) :: da_orig
    complex(8), dimension(n) :: zx_orig
    complex(8), dimension(n) :: zxb_orig
    real(4) :: temp_re, temp_im
    integer :: i, j

    nsize = n
    incx_val = 1

    call random_number(da)
    da = da * 2.0 - 1.0
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      zx(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do

    da_orig = da
    zx_orig = zx

    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      zxb(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do
    zxb_orig = zxb

    dab = 0.0

    write(*,*) 'Testing ZDSCAL (n =', n, ')'

    call zdscal_b(nsize, da, dab, zx, zxb, incx_val)

    call check_vjp_numerically(n, nsize, incx_val, da_orig, zx_orig, zxb_orig, dab, zxb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nsize, incx_val, da_orig, zx_orig, zxb_orig, dab, zxb, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nsize
    integer, intent(in) :: incx_val
    real(8), intent(in) :: da_orig
    complex(8), intent(in) :: zx_orig(n)
    complex(8), intent(in) :: zxb_orig(n)
    real(8), intent(in) :: dab
    complex(8), intent(in) :: zxb(n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j, n_products
    real(8), dimension(n) :: temp_products
    real(4) :: temp_re, temp_im

    real(8) :: da_dir
    complex(8), dimension(n) :: zx_dir

    complex(8), dimension(n) :: zx_plus, zx_minus, zx_central_diff

    real(8) :: da
    complex(8), dimension(n) :: zx

    max_error = 0.0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    call random_number(da_dir)
    da_dir = da_dir * 2.0 - 1.0
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      zx_dir(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do

    da = da_orig + h * da_dir
    zx = zx_orig + cmplx(h, 0.0) * zx_dir
    call zdscal(nsize, da, zx, incx_val)
    zx_plus = zx

    da = da_orig - h * da_dir
    zx = zx_orig - cmplx(h, 0.0) * zx_dir
    call zdscal(nsize, da, zx, incx_val)
    zx_minus = zx

    zx_central_diff = (zx_plus - zx_minus) / (2.0 * h)

    vjp_fd = 0.0
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(zxb_orig(i)) * zx_central_diff(i))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do

    vjp_ad = 0.0
    vjp_ad = vjp_ad + da_dir * dab
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(zx_dir(i)) * zxb(i))
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

end program test_zdscal_reverse