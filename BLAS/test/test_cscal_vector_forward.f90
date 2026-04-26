! Test program for CSCAL vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=n
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_cscal_vector_forward
  implicit none

  external :: cscal
  external :: cscal_dv

  integer :: nbdirs
  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(3)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing CSCAL (Vector Forward, multi-size: n = 4)'
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

    integer :: nsize, incx_val
    complex(4) :: alpha
    complex(4), dimension(n) :: x
    complex(4), dimension(nbdirs) :: alpha_dv
    complex(4), dimension(nbdirs,n) :: x_dv
    complex(4) :: alpha_orig
    complex(4), dimension(n) :: x_orig
    complex(4), dimension(nbdirs) :: alpha_dv_orig
    complex(4), dimension(nbdirs,n) :: x_dv_orig
    integer :: idir, i
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

    do idir = 1, nbdirs
      call random_number(temp_real)
      call random_number(temp_imag)
      alpha_dv(idir) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(alpha_dv))
    end do
    do idir = 1, nbdirs
      do i = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        x_dv(idir,i) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(x_dv))
      end do
    end do

    alpha_orig = alpha
    alpha_dv_orig = alpha_dv
    x_orig = x
    x_dv_orig = x_dv

    write(*,*) 'Testing CSCAL (Vector Forward, n =', n, ')'

    call cscal_dv(nsize, alpha, alpha_dv, x, x_dv, incx_val, nbdirs)

    write(*,*) 'Function calls completed successfully'

    call check_derivatives_numerically(n, nbdirs, nsize, incx_val, alpha_orig, alpha_dv_orig, x_orig, x_dv_orig, x_dv, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, nbdirs, nsize, incx_val, alpha_orig, alpha_dv_orig, x_orig, x_dv_orig, x_dv, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    integer, intent(in) :: nsize, incx_val
    complex(4), intent(in) :: alpha_orig
    complex(4), intent(in) :: alpha_dv_orig(nbdirs)
    complex(4), intent(in) :: x_orig(n), x_dv_orig(nbdirs,n)
    complex(4), intent(in) :: x_dv(nbdirs,n)
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3
    real(4) :: relative_error, max_error, abs_error, abs_reference, error_bound
    complex(4) :: central_diff, ad_result
    logical :: has_large_errors
    complex(4), dimension(n) :: x_forward, x_backward
    integer :: i, idir
    complex(4) :: alpha
    complex(4), dimension(n) :: x

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    do idir = 1, nbdirs
      alpha = alpha_orig + h * alpha_dv_orig(idir)
      x = x_orig + h * x_dv_orig(idir,:)
      call cscal(nsize, alpha, x, incx_val)
      x_forward = x
      alpha = alpha_orig - h * alpha_dv_orig(idir)
      x = x_orig - h * x_dv_orig(idir,:)
      call cscal(nsize, alpha, x, incx_val)
      x_backward = x
      do i = 1, min(4, n)
        central_diff = (x_forward(i) - x_backward(i)) / (2.0e0 * h)
        ad_result = x_dv(idir,i)
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 1.0e-3 + 1.0e-3 * abs_reference
        if (abs_error > error_bound) has_large_errors = .true.
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        max_error = max(max_error, relative_error)
      end do
    end do

    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_cscal_vector_forward