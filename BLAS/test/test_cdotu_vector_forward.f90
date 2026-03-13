! Test program for CDOTU vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=n
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_cdotu_vector_forward
  implicit none

  complex(4), external :: cdotu
  external :: cdotu_dv

  integer :: nbdirs
  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing CDOTU (Vector Forward, multi-size: n = 4)'
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
    complex(4), dimension(n) :: x, y
    complex(4), dimension(nbdirs,n) :: x_dv, y_dv
    complex(4) :: result_val
    complex(4), dimension(nbdirs) :: result_dv
    complex(4), dimension(n) :: x_orig, y_orig
    complex(4), dimension(nbdirs,n) :: x_dv_orig, y_dv_orig
    integer :: idir, i
    real(4) :: temp_real, temp_imag

    nsize = n
    incx_val = 1
    incy_val = 1

    do i = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      x(i) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(x))
      call random_number(temp_real)
      call random_number(temp_imag)
      y(i) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(y))
    end do
    do idir = 1, nbdirs
      do i = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        x_dv(idir,i) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(x_dv))
        call random_number(temp_real)
        call random_number(temp_imag)
        y_dv(idir,i) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(y_dv))
      end do
    end do

    x_orig = x
    y_orig = y
    x_dv_orig = x_dv
    y_dv_orig = y_dv

    result_val = cdotu(nsize, x, incx_val, y, incy_val)

    write(*,*) 'Testing CDOTU (Vector Forward, n =', n, ')'

    call cdotu_dv(nsize, x, x_dv, incx_val, y, y_dv, incy_val, result_val, result_dv, nbdirs)

    write(*,*) 'Function calls completed successfully'

    call check_derivatives_numerically(n, nbdirs, nsize, incx_val, incy_val, x_orig, x_dv_orig, y_orig, y_dv_orig, result_dv, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, nbdirs, nsize, incx_val, incy_val, x_orig, x_dv_orig, y_orig, y_dv_orig, result_dv, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    integer, intent(in) :: nsize, incx_val, incy_val
    complex(4), intent(in) :: x_orig(n), x_dv_orig(nbdirs,n)
    complex(4), intent(in) :: y_orig(n), y_dv_orig(nbdirs,n)
    complex(4), intent(in) :: result_dv(nbdirs)
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3
    real(4) :: relative_error, max_error, abs_error, abs_reference, error_bound
    complex(4) :: central_diff, ad_result, result_forward, result_backward
    logical :: has_large_errors
    integer :: idir
    complex(4), dimension(n) :: x, y

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    do idir = 1, nbdirs
      x = x_orig + h * x_dv_orig(idir,:)
      y = y_orig + h * y_dv_orig(idir,:)
      result_forward = cdotu(nsize, x, incx_val, y, incy_val)
      x = x_orig - h * x_dv_orig(idir,:)
      y = y_orig - h * y_dv_orig(idir,:)
      result_backward = cdotu(nsize, x, incx_val, y, incy_val)
      central_diff = (result_forward - result_backward) / (2.0e0 * h)
      ad_result = result_dv(idir)
      abs_error = abs(central_diff - ad_result)
      abs_reference = abs(ad_result)
      error_bound = 1.0e-3 + 1.0e-3 * abs_reference
      if (abs_error > error_bound) has_large_errors = .true.
      relative_error = abs_error / max(abs_reference, 1.0e-10)
      max_error = max(max_error, relative_error)
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

end program test_cdotu_vector_forward