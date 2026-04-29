! Test program for ZGEMV vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirs=n
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_zgemv_vector_reverse
  implicit none

  external :: zgemv
  external :: zgemv_bv

  integer :: nbdirs
  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(3)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing ZGEMV (Vector Reverse, multi-size: n =', test_sizes(1), ')'
  all_passed = .true.
  do i = 1, 3
    n_test = test_sizes(i)
    nbdirs = test_sizes(i)
    call run_test_for_size(n_test, passed, nbdirs)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) then
    write(*,*) 'PASS: All sizes completed successfully'
  else
    write(*,*) 'FAIL: One or more sizes had derivative errors'
  end if

contains

  subroutine run_test_for_size(n, passed, nbdirs)
    implicit none
    integer, intent(in) :: n
    logical, intent(out) :: passed
    integer, intent(in) :: nbdirs

    character :: trans
    integer :: msize, nsize, lda_val, incx_val, incy_val
    complex(8) :: alpha, beta
    complex(8), dimension(n,n) :: a
    complex(8), dimension(n) :: x, y
    complex(8), dimension(nbdirs) :: alphab, betab
    complex(8), dimension(nbdirs,n,n) :: ab
    complex(8), dimension(nbdirs,n) :: xb, yb
    complex(8) :: alpha_orig, beta_orig
    complex(8), dimension(n,n) :: a_orig
    complex(8), dimension(n) :: x_orig, y_orig
    complex(8), dimension(nbdirs,n) :: yb_orig
    integer :: k, ii, jj
    real(4) :: temp_real, temp_imag

    trans = 'N'
    msize = n
    nsize = n
    lda_val = n
    incx_val = 1
    incy_val = 1

    call random_number(temp_real)
    call random_number(temp_imag)
    alpha = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(alpha))
    do jj = 1, n
      do ii = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        a(ii,jj) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(a))
      end do
    end do
    do ii = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      x(ii) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(x))
    end do
    call random_number(temp_real)
    call random_number(temp_imag)
    beta = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(beta))
    do ii = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      y(ii) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(y))
    end do

    alpha_orig = alpha
    a_orig = a
    x_orig = x
    beta_orig = beta
    y_orig = y

    do k = 1, nbdirs
      do ii = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        yb(k,ii) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(yb))
      end do
    end do
    yb_orig = yb

    alphab = 0.0d0
    ab = 0.0d0
    xb = 0.0d0
    betab = 0.0d0

    write(*,*) 'Testing ZGEMV (Vector Reverse, n =', n, ')'

    call set_ISIZE2OFA(n)
    call set_ISIZE1OFX(n)

    call zgemv_bv(trans, msize, nsize, alpha, alphab, a, ab, lda_val, x, xb, incx_val, beta, betab, y, yb, incy_val, nbdirs)

    call set_ISIZE2OFA(-1)
    call set_ISIZE1OFX(-1)

    call check_vjp_numerically(n, nbdirs, trans, msize, nsize, lda_val, incx_val, incy_val, alpha_orig, a_orig, x_orig, beta_orig, y_orig, yb_orig, alphab, ab, xb, betab, yb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nbdirs, trans, msize, nsize, lda_val, incx_val, incy_val, alpha_orig, a_orig, x_orig, beta_orig, y_orig, yb_orig, alphab, ab, xb, betab, yb, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    character, intent(in) :: trans
    integer, intent(in) :: msize, nsize, lda_val, incx_val, incy_val
    complex(8), intent(in) :: alpha_orig, beta_orig
    complex(8), intent(in) :: a_orig(n,n)
    complex(8), intent(in) :: x_orig(n), y_orig(n)
    complex(8), intent(in) :: yb_orig(nbdirs,n)
    complex(8), intent(in) :: alphab(nbdirs), betab(nbdirs)
    complex(8), intent(in) :: ab(nbdirs,n,n), xb(nbdirs,n), yb(nbdirs,n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    complex(8) :: alpha_dir, beta_dir
    complex(8), dimension(n,n) :: a_dir
    complex(8), dimension(n) :: x_dir, y_dir
    complex(8) :: alpha, beta
    complex(8), dimension(n,n) :: a
    complex(8), dimension(n) :: x, y, y_plus, y_minus, y_central_diff
    complex(8), dimension(n) :: temp_products
    integer :: n_products, i, k, ii, jj
    real(4) :: temp_real, temp_imag
    logical :: has_large_errors

    max_error = 0.0d0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    do k = 1, nbdirs
      call random_number(temp_real)
      call random_number(temp_imag)
      alpha_dir = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(alpha_dir))
      do jj = 1, n
        do ii = 1, n
          call random_number(temp_real)
          call random_number(temp_imag)
          a_dir(ii,jj) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(a_dir))
        end do
      end do
      do ii = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        x_dir(ii) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(x_dir))
      end do
      call random_number(temp_real)
      call random_number(temp_imag)
      beta_dir = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(beta_dir))
      do ii = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        y_dir(ii) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(y_dir))
      end do
      alpha = alpha_orig + h * alpha_dir
      a = a_orig + h * a_dir
      x = x_orig + h * x_dir
      beta = beta_orig + h * beta_dir
      y = y_orig + h * y_dir
      call zgemv(trans, msize, nsize, alpha, a, lda_val, x, incx_val, beta, y, incy_val)
      y_plus = y
      alpha = alpha_orig - h * alpha_dir
      a = a_orig - h * a_dir
      x = x_orig - h * x_dir
      beta = beta_orig - h * beta_dir
      y = y_orig - h * y_dir
      call zgemv(trans, msize, nsize, alpha, a, lda_val, x, incx_val, beta, y, incy_val)
      y_minus = y
      y_central_diff = (y_plus - y_minus) / (2.0d0 * h)
      vjp_fd = 0.0d0
      n_products = 0
      do i = 1, n
        n_products = n_products + 1
        temp_products(n_products) = conjg(yb_orig(k,i)) * y_central_diff(i)
        vjp_fd = vjp_fd + real(temp_products(n_products))
      end do
      vjp_ad = 0.0d0
      vjp_ad = vjp_ad + real(conjg(beta_dir) * betab(k))
      vjp_ad = vjp_ad + real(conjg(alpha_dir) * alphab(k))
      do jj = 1, n
        do ii = 1, n
          vjp_ad = vjp_ad + real(conjg(a_dir(ii,jj)) * ab(k,ii,jj))
        end do
      end do
      do ii = 1, n
        vjp_ad = vjp_ad + real(conjg(x_dir(ii)) * xb(k,ii))
        vjp_ad = vjp_ad + real(conjg(y_dir(ii)) * yb(k,ii))
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
      if (relative_error > max_error) max_error = relative_error
    end do
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_vjp_numerically

end program test_zgemv_vector_reverse