! Test program for ZGERC vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirs=n
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_zgerc_vector_forward
  implicit none

  external :: zgerc
  external :: zgerc_dv

  integer :: nbdirs
  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing ZGERC (Vector Forward, multi-size: n = 4)'
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

    integer :: msize, nsize, lda_val, incx_val, incy_val
    complex(8) :: alpha
    complex(8), dimension(n) :: x, y
    complex(8), dimension(n,n) :: a
    complex(8), dimension(nbdirs) :: alpha_dv, alpha_dv_orig
    complex(8), dimension(nbdirs,n) :: x_dv, y_dv
    complex(8), dimension(nbdirs,n,n) :: a_dv
    complex(8) :: alpha_orig
    complex(8), dimension(n,n) :: a_orig
    complex(8), dimension(nbdirs,n,n) :: a_dv_orig
    complex(8), dimension(n) :: x_orig, y_orig
    complex(8), dimension(nbdirs,n) :: x_dv_orig, y_dv_orig
    integer :: idir, ii, jj
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
    end do
    do ii = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      y(ii) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(y))
    end do
    do jj = 1, n
      do ii = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        a(ii,jj) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(a))
      end do
    end do
    do idir = 1, nbdirs
      call random_number(temp_real)
      call random_number(temp_imag)
      alpha_dv(idir) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(alpha_dv))
    end do
    do idir = 1, nbdirs
      do ii = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        x_dv(idir,ii) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(x_dv))
        call random_number(temp_real)
        call random_number(temp_imag)
        y_dv(idir,ii) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(y_dv))
      end do
    end do
    do idir = 1, nbdirs
      do jj = 1, n
        do ii = 1, n
          call random_number(temp_real)
          call random_number(temp_imag)
          a_dv(idir,ii,jj) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(a_dv))
        end do
      end do
    end do

    alpha_orig = alpha
    alpha_dv_orig = alpha_dv
    a_orig = a
    a_dv_orig = a_dv
    x_orig = x
    x_dv_orig = x_dv
    y_orig = y
    y_dv_orig = y_dv

    write(*,*) 'Testing ZGERC (Vector Forward, n =', n, ')'

    call zgerc_dv(msize, nsize, alpha, alpha_dv, x, x_dv, incx_val, y, y_dv, incy_val, a, a_dv, lda_val, nbdirs)

    write(*,*) 'Function calls completed successfully'

    call check_derivatives_numerically(n, nbdirs, msize, nsize, lda_val, incx_val, incy_val, alpha_orig, alpha_dv_orig, x_orig, x_dv_orig, y_orig, y_dv_orig, a_orig, a_dv_orig, a_dv, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, nbdirs, msize, nsize, lda_val, incx_val, incy_val, alpha_orig, alpha_dv_orig, x_orig, x_dv_orig, y_orig, y_dv_orig, a_orig, a_dv_orig, a_dv, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    integer, intent(in) :: msize, nsize, lda_val, incx_val, incy_val
    complex(8), intent(in) :: alpha_orig
    complex(8), intent(in) :: alpha_dv_orig(nbdirs)
    complex(8), intent(in) :: x_orig(n), x_dv_orig(nbdirs,n)
    complex(8), intent(in) :: y_orig(n), y_dv_orig(nbdirs,n)
    complex(8), intent(in) :: a_orig(n,n), a_dv_orig(nbdirs,n,n)
    complex(8), intent(in) :: a_dv(nbdirs,n,n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: relative_error, max_error, abs_error, abs_reference, error_bound
    complex(8) :: central_diff, ad_result
    logical :: has_large_errors
    complex(8), dimension(n,n) :: a_forward, a_backward
    integer :: i, j, idir
    complex(8) :: alpha
    complex(8), dimension(n) :: x, y
    complex(8), dimension(n,n) :: a

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    do idir = 1, nbdirs
      alpha = alpha_orig + h * alpha_dv_orig(idir)
      x = x_orig + h * x_dv_orig(idir,:)
      y = y_orig + h * y_dv_orig(idir,:)
      a = a_orig + h * a_dv_orig(idir,:,:)
      call zgerc(msize, nsize, alpha, x, incx_val, y, incy_val, a, lda_val)
      a_forward = a
      alpha = alpha_orig - h * alpha_dv_orig(idir)
      x = x_orig - h * x_dv_orig(idir,:)
      y = y_orig - h * y_dv_orig(idir,:)
      a = a_orig - h * a_dv_orig(idir,:,:)
      call zgerc(msize, nsize, alpha, x, incx_val, y, incy_val, a, lda_val)
      a_backward = a
      do j = 1, min(4, n)
        do i = 1, min(4, n)
          central_diff = (a_forward(i,j) - a_backward(i,j)) / (2.0e0 * h)
          ad_result = a_dv(idir,i,j)
          abs_error = abs(central_diff - ad_result)
          abs_reference = abs(ad_result)
          error_bound = 1.0e-5 + 1.0e-5 * abs_reference
          if (abs_error > error_bound) has_large_errors = .true.
          relative_error = abs_error / max(abs_reference, 1.0e-10)
          max_error = max(max_error, relative_error)
        end do
      end do
    end do

    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_zgerc_vector_forward