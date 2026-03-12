! Test program for DCOPY vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirs=n
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_dcopy_vector_forward
  implicit none

  external :: dcopy
  external :: dcopy_dv

  integer :: nbdirs
  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing DCOPY (Vector Forward, multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    nbdirs = test_sizes(i)
    call run_test_for_size(n_test, passed, nbdirs)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) then
    write(*,*) 'PASS: Vector forward mode - all sizes completed successfully'
  else
    write(*,*) 'FAIL: Vector forward mode - one or more sizes had derivative errors'
  end if

contains

  subroutine run_test_for_size(n, passed, nbdirs)
    implicit none
    integer, intent(in) :: n
    logical, intent(out) :: passed
    integer, intent(in) :: nbdirs

    integer :: nsize, incx_val, incy_val
    real(8), dimension(n) :: x, y
    real(8), dimension(nbdirs,n) :: x_dv, y_dv
    real(8), dimension(n) :: x_orig, y_orig
    real(8), dimension(nbdirs,n) :: x_dv_orig, y_dv_orig
    integer :: idir, i
    real(4) :: temp_real, temp_imag

    nsize = n
    incx_val = 1
    incy_val = 1

    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    call random_number(y)
    y = y * 2.0d0 - 1.0d0

    do idir = 1, nbdirs
      call random_number(x_dv(idir,:))
      x_dv(idir,:) = x_dv(idir,:) * 2.0d0 - 1.0d0
    end do
    do idir = 1, nbdirs
      call random_number(y_dv(idir,:))
      y_dv(idir,:) = y_dv(idir,:) * 2.0d0 - 1.0d0
    end do

    x_orig = x
    x_dv_orig = x_dv
    y_orig = y
    y_dv_orig = y_dv

    write(*,*) 'Testing DCOPY (Vector Forward, n =', n, ')'

    call set_ISIZE1OFDy(n)

    call dcopy_dv(nsize, x, x_dv, incx_val, y, y_dv, incy_val, nbdirs)

    call set_ISIZE1OFDy(-1)

    write(*,*) 'Function calls completed successfully'

    call check_derivatives_numerically(n, nbdirs, nsize, incx_val, incy_val, x_orig, x_dv_orig, y_orig, y_dv_orig, y_dv, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, nbdirs, nsize, incx_val, incy_val, x_orig, x_dv_orig, y_orig, y_dv_orig, y_dv, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    integer, intent(in) :: nsize, incx_val, incy_val
    real(8), intent(in) :: x_orig(n), x_dv_orig(nbdirs,n)
    real(8), intent(in) :: y_orig(n), y_dv_orig(nbdirs,n)
    real(8), intent(in) :: y_dv(nbdirs,n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: relative_error, max_error, abs_error, abs_reference, error_bound
    real(8) :: central_diff, ad_result
    logical :: has_large_errors
    real(8), dimension(n) :: y_forward, y_backward
    integer :: i, idir
    real(8), dimension(n) :: x, y

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    do idir = 1, nbdirs
      x = x_orig + h * x_dv_orig(idir,:)
      y = y_orig
      call dcopy(nsize, x, incx_val, y, incy_val)
      y_forward = y
      x = x_orig - h * x_dv_orig(idir,:)
      y = y_orig
      call dcopy(nsize, x, incx_val, y, incy_val)
      y_backward = y
      do i = 1, min(4, n)
        central_diff = (y_forward(i) - y_backward(i)) / (2.0e0 * h)
        ad_result = y_dv(idir,i)
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 1.0e-5 + 1.0e-5 * abs_reference
        if (abs_error > error_bound) has_large_errors = .true.
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        max_error = max(max_error, relative_error)
      end do
    end do

    write(*,*) 'Maximum relative error across all directions:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in vector derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Vector derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_dcopy_vector_forward