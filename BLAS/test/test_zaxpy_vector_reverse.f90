! Test program for ZAXPY vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirs=n
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_zaxpy_vector_reverse
  implicit none

  external :: zaxpy
  external :: zaxpy_bv

  integer :: nbdirs
  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing ZAXPY (Vector Reverse, multi-size: n =', test_sizes(1), ')'
  all_passed = .true.
  do i = 1, 1
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

    integer :: nsize, incx_val, incy_val
    complex(8) :: alpha
    complex(8), dimension(n) :: x, y
    complex(8), dimension(nbdirs) :: alphab
    complex(8), dimension(nbdirs,n) :: xb, yb
    complex(8) :: alpha_orig
    complex(8), dimension(n) :: x_orig, y_orig
    complex(8), dimension(nbdirs,n) :: yb_orig
    integer :: k, i
    real(4) :: temp_real, temp_imag

    nsize = n
    incx_val = 1
    incy_val = 1

    call random_number(temp_real)
    call random_number(temp_imag)
    alpha = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(alpha))
    do i = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      x(i) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(x))
      call random_number(temp_real)
      call random_number(temp_imag)
      y(i) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(y))
    end do

    alpha_orig = alpha
    x_orig = x
    y_orig = y

    do k = 1, nbdirs
      do i = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        yb(k,i) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(yb))
      end do
    end do
    yb_orig = yb

    alphab = 0.0d0
    xb = 0.0d0

    write(*,*) 'Testing ZAXPY (Vector Reverse, n =', n, ')'

    ! Set ISIZE globals required by AXPY bv routine (dimension 1 of vectors).
    call set_ISIZE1OFZx(n)

    call zaxpy_bv(nsize, alpha, alphab, x, xb, incx_val, y, yb, incy_val, nbdirs)

    call set_ISIZE1OFZx(-1)

    call check_vjp_numerically(n, nbdirs, nsize, incx_val, incy_val, alpha_orig, x_orig, y_orig, yb_orig, alphab, xb, yb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nbdirs, nsize, incx_val, incy_val, alpha_orig, x_orig, y_orig, yb_orig, alphab, xb, yb, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    integer, intent(in) :: nsize, incx_val, incy_val
    complex(8), intent(in) :: alpha_orig
    complex(8), intent(in) :: x_orig(n), y_orig(n)
    complex(8), intent(in) :: yb_orig(nbdirs,n)
    complex(8), intent(in) :: alphab(nbdirs)
    complex(8), intent(in) :: xb(nbdirs,n), yb(nbdirs,n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    complex(8) :: alpha_dir
    complex(8), dimension(n) :: x_dir, y_dir
    complex(8) :: alpha
    complex(8), dimension(n) :: x, y, y_plus, y_minus, y_central_diff
    complex(8), dimension(n) :: temp_products
    integer :: n_products, i, k
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
      do i = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        x_dir(i) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(x_dir))
        call random_number(temp_real)
        call random_number(temp_imag)
        y_dir(i) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(y_dir))
      end do
      alpha = alpha_orig + h * alpha_dir
      x = x_orig + h * x_dir
      y = y_orig + h * y_dir
      call zaxpy(nsize, alpha, x, incx_val, y, incy_val)
      y_plus = y
      alpha = alpha_orig - h * alpha_dir
      x = x_orig - h * x_dir
      y = y_orig - h * y_dir
      call zaxpy(nsize, alpha, x, incx_val, y, incy_val)
      y_minus = y
      y_central_diff = (y_plus - y_minus) / (2.0d0 * h)
      vjp_fd = 0.0d0
      n_products = 0
      do i = 1, n
        temp_products(i) = conjg(yb_orig(k,i)) * y_central_diff(i)
        vjp_fd = vjp_fd + real(temp_products(i))
      end do
      vjp_ad = 0.0d0
      vjp_ad = vjp_ad + real(conjg(alpha_dir) * alphab(k))
      do i = 1, n
        vjp_ad = vjp_ad + real(conjg(x_dir(i)) * xb(k,i))
        vjp_ad = vjp_ad + real(conjg(y_dir(i)) * yb(k,i))
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

end program test_zaxpy_vector_reverse