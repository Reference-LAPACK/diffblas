! Test program for ZGEMV reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_zgemv_reverse
  implicit none

  external :: zgemv
  external :: zgemv_b

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing ZGEMV (multi-size: n = 4)'
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

    character :: trans
    integer :: msize
    integer :: nsize
    complex(8) :: alpha
    complex(8), dimension(n,n) :: a
    integer :: lda_val
    complex(8), dimension(n) :: x
    integer :: incx_val
    complex(8) :: beta
    complex(8), dimension(n) :: y
    integer :: incy_val
    complex(8) :: alphab
    complex(8), dimension(n,n) :: ab
    complex(8), dimension(n) :: xb
    complex(8) :: betab
    complex(8), dimension(n) :: yb
    complex(8) :: alpha_orig
    complex(8), dimension(n,n) :: a_orig
    complex(8), dimension(n) :: x_orig
    complex(8) :: beta_orig
    complex(8), dimension(n) :: y_orig
    complex(8), dimension(n) :: yb_orig
    real(4) :: temp_re, temp_im
    integer :: i, j

    nsize = n
    msize = n
    lda_val = n
    incx_val = 1
    incy_val = 1
    trans = 'N'

    call random_number(temp_re)
    call random_number(temp_im)
    alpha = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    do j = 1, n
      do i = 1, n
        call random_number(temp_re)
        call random_number(temp_im)
        a(i,j) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
      end do
    end do
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      x(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do
    call random_number(temp_re)
    call random_number(temp_im)
    beta = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      y(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do

    alpha_orig = alpha
    a_orig = a
    x_orig = x
    beta_orig = beta
    y_orig = y

    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      yb(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do
    yb_orig = yb

    alphab = 0.0
    ab = 0.0
    xb = 0.0
    betab = 0.0

    write(*,*) 'Testing ZGEMV (n =', n, ')'

    call set_ISIZE1OFX(n)
    call set_ISIZE2OFA(n)

    call zgemv_b(trans, msize, nsize, alpha, alphab, a, ab, lda_val, x, xb, incx_val, beta, betab, y, yb, incy_val)

    call set_ISIZE1OFX(-1)
    call set_ISIZE2OFA(-1)

    call check_vjp_numerically(n, trans, msize, nsize, lda_val, incx_val, incy_val, alpha_orig, a_orig, x_orig, beta_orig, y_orig, yb_orig, alphab, ab, xb, betab, yb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, trans, msize, nsize, lda_val, incx_val, incy_val, alpha_orig, a_orig, x_orig, beta_orig, y_orig, yb_orig, alphab, ab, xb, betab, yb, passed)
    implicit none
    integer, intent(in) :: n
    character, intent(in) :: trans
    integer, intent(in) :: msize
    integer, intent(in) :: nsize
    integer, intent(in) :: lda_val
    integer, intent(in) :: incx_val
    integer, intent(in) :: incy_val
    complex(8), intent(in) :: alpha_orig
    complex(8), intent(in) :: a_orig(n,n)
    complex(8), intent(in) :: x_orig(n)
    complex(8), intent(in) :: beta_orig
    complex(8), intent(in) :: y_orig(n)
    complex(8), intent(in) :: yb_orig(n)
    complex(8), intent(in) :: alphab
    complex(8), intent(in) :: ab(n,n)
    complex(8), intent(in) :: xb(n)
    complex(8), intent(in) :: betab
    complex(8), intent(in) :: yb(n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j, n_products
    real(8), dimension(n) :: temp_products
    real(4) :: temp_re, temp_im

    complex(8) :: alpha_dir
    complex(8), dimension(n,n) :: a_dir
    complex(8), dimension(n) :: x_dir
    complex(8) :: beta_dir
    complex(8), dimension(n) :: y_dir

    complex(8), dimension(n) :: y_plus, y_minus, y_central_diff

    complex(8) :: alpha
    complex(8), dimension(n,n) :: a
    complex(8), dimension(n) :: x
    complex(8) :: beta
    complex(8), dimension(n) :: y

    max_error = 0.0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    call random_number(temp_re)
    call random_number(temp_im)
    alpha_dir = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    do j = 1, n
      do i = 1, n
        call random_number(temp_re)
        call random_number(temp_im)
        a_dir(i,j) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
      end do
    end do
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      x_dir(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do
    call random_number(temp_re)
    call random_number(temp_im)
    beta_dir = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      y_dir(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do

    alpha = alpha_orig + cmplx(h, 0.0) * alpha_dir
    a = a_orig + cmplx(h, 0.0) * a_dir
    x = x_orig + cmplx(h, 0.0) * x_dir
    beta = beta_orig + cmplx(h, 0.0) * beta_dir
    y = y_orig + cmplx(h, 0.0) * y_dir
    call zgemv(trans, msize, nsize, alpha, a, lda_val, x, incx_val, beta, y, incy_val)
    y_plus = y

    alpha = alpha_orig - cmplx(h, 0.0) * alpha_dir
    a = a_orig - cmplx(h, 0.0) * a_dir
    x = x_orig - cmplx(h, 0.0) * x_dir
    beta = beta_orig - cmplx(h, 0.0) * beta_dir
    y = y_orig - cmplx(h, 0.0) * y_dir
    call zgemv(trans, msize, nsize, alpha, a, lda_val, x, incx_val, beta, y, incy_val)
    y_minus = y

    y_central_diff = (y_plus - y_minus) / (2.0 * h)

    vjp_fd = 0.0
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(yb_orig(i)) * y_central_diff(i))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do

    vjp_ad = 0.0
    vjp_ad = vjp_ad + real(conjg(alpha_dir) * alphab)
    do j = 1, n
      do i = 1, n
        vjp_ad = vjp_ad + real(conjg(a_dir(i,j)) * ab(i,j))
      end do
    end do
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(x_dir(i)) * xb(i))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    vjp_ad = vjp_ad + real(conjg(beta_dir) * betab)
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(y_dir(i)) * yb(i))
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

end program test_zgemv_reverse