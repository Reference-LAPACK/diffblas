! Test program for CGERC reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_cgerc_reverse
  implicit none

  external :: cgerc
  external :: cgerc_b

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(3)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing CGERC (multi-size: n = 4)'
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

    integer :: msize
    integer :: nsize
    complex(4) :: alpha
    complex(4), dimension(n) :: x
    integer :: incx_val
    complex(4), dimension(n) :: y
    integer :: incy_val
    complex(4), dimension(n,n) :: a
    integer :: lda_val
    complex(4) :: alphab
    complex(4), dimension(n) :: xb
    complex(4), dimension(n) :: yb
    complex(4), dimension(n,n) :: ab
    complex(4) :: alpha_orig
    complex(4), dimension(n) :: x_orig
    complex(4), dimension(n) :: y_orig
    complex(4), dimension(n,n) :: a_orig
    complex(4), dimension(n,n) :: ab_orig
    real(4) :: temp_re, temp_im
    integer :: i, j

    nsize = n
    msize = n
    incx_val = 1
    incy_val = 1
    lda_val = n

    call random_number(temp_re)
    call random_number(temp_im)
    alpha = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      x(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      y(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do
    do j = 1, n
      do i = 1, n
        call random_number(temp_re)
        call random_number(temp_im)
        a(i,j) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
      end do
    end do

    alpha_orig = alpha
    x_orig = x
    y_orig = y
    a_orig = a

    call random_number(temp_re)
    call random_number(temp_im)
    ab = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    ab_orig = ab

    alphab = 0.0
    xb = 0.0
    yb = 0.0

    write(*,*) 'Testing CGERC (n =', n, ')'

    call set_ISIZE1OFX(n)
    call set_ISIZE1OFY(n)

    call cgerc_b(msize, nsize, alpha, alphab, x, xb, incx_val, y, yb, incy_val, a, ab, lda_val)

    call set_ISIZE1OFX(-1)
    call set_ISIZE1OFY(-1)

    call check_vjp_numerically(n, msize, nsize, incx_val, incy_val, lda_val, alpha_orig, x_orig, y_orig, a_orig, ab_orig, alphab, xb, yb, ab, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, msize, nsize, incx_val, incy_val, lda_val, alpha_orig, x_orig, y_orig, a_orig, ab_orig, alphab, xb, yb, ab, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: msize
    integer, intent(in) :: nsize
    integer, intent(in) :: incx_val
    integer, intent(in) :: incy_val
    integer, intent(in) :: lda_val
    complex(4), intent(in) :: alpha_orig
    complex(4), intent(in) :: x_orig(n)
    complex(4), intent(in) :: y_orig(n)
    complex(4), intent(in) :: a_orig(n,n)
    complex(4), intent(in) :: ab_orig(n,n)
    complex(4), intent(in) :: alphab
    complex(4), intent(in) :: xb(n)
    complex(4), intent(in) :: yb(n)
    complex(4), intent(in) :: ab(n,n)
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j, n_products
    real(4), dimension(n) :: temp_products
    real(4) :: temp_re, temp_im

    complex(4) :: alpha_dir
    complex(4), dimension(n) :: x_dir
    complex(4), dimension(n) :: y_dir
    complex(4), dimension(n,n) :: a_dir

    complex(4), dimension(n,n) :: a_plus, a_minus, a_central_diff

    complex(4) :: alpha
    complex(4), dimension(n) :: x
    complex(4), dimension(n) :: y
    complex(4), dimension(n,n) :: a

    max_error = 0.0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    call random_number(temp_re)
    call random_number(temp_im)
    alpha_dir = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      x_dir(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      y_dir(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do
    do j = 1, n
      do i = 1, n
        call random_number(temp_re)
        call random_number(temp_im)
        a_dir(i,j) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
      end do
    end do

    alpha = alpha_orig + cmplx(h, 0.0) * alpha_dir
    x = x_orig + cmplx(h, 0.0) * x_dir
    y = y_orig + cmplx(h, 0.0) * y_dir
    a = a_orig + cmplx(h, 0.0) * a_dir
    call cgerc(msize, nsize, alpha, x, incx_val, y, incy_val, a, lda_val)
    a_plus = a

    alpha = alpha_orig - cmplx(h, 0.0) * alpha_dir
    x = x_orig - cmplx(h, 0.0) * x_dir
    y = y_orig - cmplx(h, 0.0) * y_dir
    a = a_orig - cmplx(h, 0.0) * a_dir
    call cgerc(msize, nsize, alpha, x, incx_val, y, incy_val, a, lda_val)
    a_minus = a

    a_central_diff = (a_plus - a_minus) / (2.0 * h)

    vjp_fd = 0.0
    do j = 1, n
      do i = 1, n
        vjp_fd = vjp_fd + real(conjg(ab_orig(i,j)) * a_central_diff(i,j))
      end do
    end do

    vjp_ad = 0.0
    vjp_ad = vjp_ad + real(conjg(alpha_dir) * alphab)
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(x_dir(i)) * xb(i))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(y_dir(i)) * yb(i))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    do j = 1, n
      do i = 1, n
        vjp_ad = vjp_ad + real(conjg(a_dir(i,j)) * ab(i,j))
      end do
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

end program test_cgerc_reverse