! Test program for CDOTC reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_cdotc_reverse
  implicit none

  complex(4), external :: cdotc
  external :: cdotc_b

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing CDOTC (multi-size: n = 4)'
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
    complex(4), dimension(n) :: cx
    integer :: incx_val
    complex(4), dimension(n) :: cy
    integer :: incy_val
    complex(4), dimension(n) :: cxb
    complex(4), dimension(n) :: cyb
    complex(4) :: cdotcb, cdotcb_orig
    complex(4), dimension(n) :: cx_orig
    complex(4), dimension(n) :: cy_orig
    real(4) :: temp_re, temp_im
    integer :: i, j

    nsize = n
    incx_val = 1
    incy_val = 1

    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      cx(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      cy(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do

    cx_orig = cx
    cy_orig = cy


    call random_number(temp_re)
    call random_number(temp_im)
    cdotcb = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    cdotcb_orig = cdotcb

    cxb = 0.0
    cyb = 0.0

    write(*,*) 'Testing CDOTC (n =', n, ')'

    call set_ISIZE1OFCx(n)
    call set_ISIZE1OFCy(n)

    call cdotc_b(nsize, cx, cxb, incx_val, cy, cyb, incy_val, cdotcb)

    call set_ISIZE1OFCx(-1)
    call set_ISIZE1OFCy(-1)

    call check_vjp_numerically(n, nsize, incx_val, incy_val, cx_orig, cy_orig, cxb, cyb, cdotcb_orig, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nsize, incx_val, incy_val, cx_orig, cy_orig, cxb, cyb, cdotcb_orig, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nsize
    integer, intent(in) :: incx_val
    integer, intent(in) :: incy_val
    complex(4), intent(in) :: cx_orig(n)
    complex(4), intent(in) :: cy_orig(n)
    complex(4), intent(in) :: cxb(n)
    complex(4), intent(in) :: cyb(n)
    complex(4), intent(in) :: cdotcb_orig
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j, n_products
    real(4), dimension(n) :: temp_products
    real(4) :: temp_re, temp_im

    complex(4), dimension(n) :: cx_dir
    complex(4), dimension(n) :: cy_dir

    complex(4) :: cdotc_plus, cdotc_minus

    complex(4), dimension(n) :: cx
    complex(4), dimension(n) :: cy

    max_error = 0.0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      cx_dir(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      cy_dir(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do

    cx = cx_orig + cmplx(h, 0.0) * cx_dir
    cy = cy_orig + cmplx(h, 0.0) * cy_dir
    cdotc_plus = cdotc(nsize, cx, incx_val, cy, incy_val)

    cx = cx_orig - cmplx(h, 0.0) * cx_dir
    cy = cy_orig - cmplx(h, 0.0) * cy_dir
    cdotc_minus = cdotc(nsize, cx, incx_val, cy, incy_val)


    vjp_fd = real(conjg(cdotcb_orig) * (cdotc_plus - cdotc_minus) / (2.0 * h))

    vjp_ad = 0.0
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(cx_dir(i)) * cxb(i))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(cy_dir(i)) * cyb(i))
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

end program test_cdotc_reverse