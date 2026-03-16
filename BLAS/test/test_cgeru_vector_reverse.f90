! Test program for CGERU vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=n
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_cgeru_vector_reverse
  implicit none

  external :: cgeru
  external :: cgeru_bv

  integer :: nbdirs
  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(3)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing CGERU (Vector Reverse, multi-size: n =', test_sizes(1), ')'
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

    integer :: msize, nsize, lda_val, incx_val, incy_val
    complex(4) :: alpha
    complex(4), dimension(n) :: x, y
    complex(4), dimension(n,n) :: a
    complex(4), dimension(nbdirs) :: alphab
    complex(4), dimension(nbdirs,n) :: xb, yb
    complex(4), dimension(nbdirs,n,n) :: ab
    complex(4) :: alpha_orig
    complex(4), dimension(n,n) :: a_orig
    complex(4), dimension(n) :: x_orig, y_orig
    complex(4), dimension(nbdirs,n,n) :: ab_orig
    integer :: k, ii, jj
    real(4) :: temp_real, temp_imag

    msize = n
    nsize = n
    lda_val = n
    incx_val = 1
    incy_val = 1

    call random_number(temp_real)
    call random_number(temp_imag)
    alpha = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(alpha))
    do ii = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      x(ii) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(x))
      y(ii) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(y))
    end do
    do jj = 1, n
      do ii = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        a(ii,jj) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(a))
      end do
    end do
    do k = 1, nbdirs
      do jj = 1, n
        do ii = 1, n
          call random_number(temp_real)
          call random_number(temp_imag)
          ab(k,ii,jj) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(ab))
        end do
      end do
    end do

    alpha_orig = alpha
    a_orig = a
    x_orig = x
    y_orig = y
    ab_orig = ab

    alphab = 0.0d0
    xb = 0.0d0
    yb = 0.0d0

    write(*,*) 'Testing CGERU (Vector Reverse, n =', n, ')'

    ! Set ISIZE globals required by GER bv routine (dimension 1 of vectors).
    call set_ISIZE1OFX(n)
    call set_ISIZE1OFY(n)

    call cgeru_bv(msize, nsize, alpha, alphab, x, xb, incx_val, y, yb, incy_val, a, ab, lda_val, nbdirs)

    call set_ISIZE1OFX(-1)
    call set_ISIZE1OFY(-1)

    call check_vjp_numerically(n, nbdirs, msize, nsize, lda_val, incx_val, incy_val, alpha_orig, x_orig, y_orig, a_orig, ab_orig, alphab, xb, yb, ab, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nbdirs, msize, nsize, lda_val, incx_val, incy_val, alpha_orig, x_orig, y_orig, a_orig, ab_orig, alphab, xb, yb, ab, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    integer, intent(in) :: msize, nsize, lda_val, incx_val, incy_val
    complex(4), intent(in) :: alpha_orig
    complex(4), intent(in) :: x_orig(n), y_orig(n)
    complex(4), intent(in) :: a_orig(n,n)
    complex(4), intent(in) :: ab_orig(nbdirs,n,n)
    complex(4), intent(in) :: alphab(nbdirs)
    complex(4), intent(in) :: xb(nbdirs,n), yb(nbdirs,n)
    complex(4), intent(in) :: ab(nbdirs,n,n)
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    complex(4) :: alpha_dir
    complex(4), dimension(n) :: x_dir, y_dir
    complex(4), dimension(n,n) :: a_dir
    complex(4) :: alpha
    complex(4), dimension(n) :: x, y
    complex(4), dimension(n,n) :: a, a_plus, a_minus, a_central_diff
    integer :: i, j, k, ii, jj
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
      do ii = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        x_dir(ii) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(x_dir))
        y_dir(ii) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(y_dir))
      end do
      do jj = 1, n
        do ii = 1, n
          call random_number(temp_real)
          call random_number(temp_imag)
          a_dir(ii,jj) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(a_dir))
        end do
      end do
      alpha = alpha_orig + h * alpha_dir
      x = x_orig + h * x_dir
      y = y_orig + h * y_dir
      a = a_orig + h * a_dir
      call cgeru(msize, nsize, alpha, x, incx_val, y, incy_val, a, lda_val)
      a_plus = a
      alpha = alpha_orig - h * alpha_dir
      x = x_orig - h * x_dir
      y = y_orig - h * y_dir
      a = a_orig - h * a_dir
      call cgeru(msize, nsize, alpha, x, incx_val, y, incy_val, a, lda_val)
      a_minus = a
      a_central_diff = (a_plus - a_minus) / (2.0d0 * h)
      vjp_fd = 0.0d0
      do jj = 1, n
        do ii = 1, n
          vjp_fd = vjp_fd + real(conjg(ab_orig(k,ii,jj)) * a_central_diff(ii,jj))
        end do
      end do
      vjp_ad = 0.0d0
      vjp_ad = vjp_ad + real(conjg(alpha_dir) * alphab(k))
      do ii = 1, n
        vjp_ad = vjp_ad + real(conjg(x_dir(ii)) * xb(k,ii))
        vjp_ad = vjp_ad + real(conjg(y_dir(ii)) * yb(k,ii))
      end do
      do jj = 1, n
        do ii = 1, n
          vjp_ad = vjp_ad + real(conjg(a_dir(ii,jj)) * ab(k,ii,jj))
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
      if (relative_error > max_error) max_error = relative_error
    end do

    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_vjp_numerically

end program test_cgeru_vector_reverse