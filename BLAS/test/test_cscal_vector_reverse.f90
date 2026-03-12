! Test program for CSCAL vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=n
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_cscal_vector_reverse
  implicit none

  external :: cscal
  external :: cscal_bv

  integer :: nbdirs
  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing CSCAL (Vector Reverse, multi-size: n =', test_sizes(1), ')'
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    nbdirs = test_sizes(i)
    call run_test_for_size(n_test, passed, nbdirs)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) then
    write(*,*) 'PASS: Vector reverse mode - all sizes completed successfully'
  else
    write(*,*) 'FAIL: Vector reverse mode - one or more sizes had derivative errors'
  end if

contains

  subroutine run_test_for_size(n, passed, nbdirs)
    implicit none
    integer, intent(in) :: n
    logical, intent(out) :: passed
    integer, intent(in) :: nbdirs

    integer :: nsize, incx_val
    complex(4) :: alpha
    complex(4), dimension(n) :: x
    complex(4), dimension(nbdirs) :: alphab
    complex(4), dimension(nbdirs,n) :: xb
    complex(4) :: alpha_orig
    complex(4), dimension(n) :: x_orig
    complex(4), dimension(nbdirs,n) :: xb_orig
    integer :: k, i
    real(4) :: temp_real, temp_imag

    nsize = n
    incx_val = 1

    call random_number(temp_real)
    call random_number(temp_imag)
    alpha = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(alpha))
    do i = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      x(i) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(x))
    end do

    alpha_orig = alpha
    x_orig = x

    do k = 1, nbdirs
      do i = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        xb(k,i) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(xb))
      end do
    end do
    xb_orig = xb

    alphab = 0.0d0

    write(*,*) 'Testing CSCAL (Vector Reverse, n =', n, ')'

    call cscal_bv(nsize, alpha, alphab, x, xb, incx_val, nbdirs)

    call check_vjp_numerically(n, nbdirs, nsize, incx_val, alpha_orig, x_orig, xb_orig, alphab, xb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nbdirs, nsize, incx_val, alpha_orig, x_orig, xb_orig, alphab, xb, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    integer, intent(in) :: nsize, incx_val
    complex(4), intent(in) :: alpha_orig
    complex(4), intent(in) :: x_orig(n)
    complex(4), intent(in) :: xb_orig(nbdirs,n)
    complex(4), intent(in) :: alphab(nbdirs)
    complex(4), intent(in) :: xb(nbdirs,n)
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    complex(4) :: alpha_dir
    complex(4), dimension(n) :: x_dir
    complex(4) :: alpha
    complex(4), dimension(n) :: x, x_plus, x_minus, x_central_diff
    complex(4), dimension(n) :: temp_products
    integer :: i, k
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
      end do
      alpha = alpha_orig + h * alpha_dir
      x = x_orig + h * x_dir
      call cscal(nsize, alpha, x, incx_val)
      x_plus = x
      alpha = alpha_orig - h * alpha_dir
      x = x_orig - h * x_dir
      call cscal(nsize, alpha, x, incx_val)
      x_minus = x
      x_central_diff = (x_plus - x_minus) / (2.0d0 * h)
      vjp_fd = 0.0d0
      do i = 1, n
        temp_products(i) = conjg(xb_orig(k,i)) * x_central_diff(i)
        vjp_fd = vjp_fd + real(temp_products(i))
      end do
      vjp_ad = 0.0d0
      vjp_ad = vjp_ad + real(conjg(alpha_dir) * alphab(k))
      do i = 1, n
        vjp_ad = vjp_ad + real(conjg(x_dir(i)) * xb(k,i))
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

    write(*,*) ''
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_vjp_numerically

end program test_cscal_vector_reverse