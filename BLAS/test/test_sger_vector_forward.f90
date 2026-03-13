! Test program for SGER vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=n
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_sger_vector_forward
  implicit none

  external :: sger
  external :: sger_dv

  integer :: nbdirs
  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing SGER (Vector Forward, multi-size: n = 4)'
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
    real(4) :: alpha
    real(4), dimension(n) :: x, y
    real(4), dimension(n,n) :: a
    real(4), dimension(nbdirs) :: alpha_dv, alpha_dv_orig
    real(4), dimension(nbdirs,n) :: x_dv, y_dv
    real(4), dimension(nbdirs,n,n) :: a_dv
    real(4) :: alpha_orig
    real(4), dimension(n,n) :: a_orig
    real(4), dimension(nbdirs,n,n) :: a_dv_orig
    real(4), dimension(n) :: x_orig, y_orig
    real(4), dimension(nbdirs,n) :: x_dv_orig, y_dv_orig
    integer :: idir, ii, jj
    real(4) :: temp_real, temp_imag

    msize = n
    nsize = n
    lda_val = n
    incx_val = 1
    incy_val = 1

    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    call random_number(y)
    y = y * 2.0d0 - 1.0d0
    call random_number(a)
    a = a * 2.0d0 - 1.0d0

    do idir = 1, nbdirs
      call random_number(temp_real)
      alpha_dv(idir) = temp_real * 2.0d0 - 1.0d0
    end do
    do idir = 1, nbdirs
      call random_number(x_dv(idir,:))
      x_dv(idir,:) = x_dv(idir,:) * 2.0d0 - 1.0d0
      call random_number(y_dv(idir,:))
      y_dv(idir,:) = y_dv(idir,:) * 2.0d0 - 1.0d0
      call random_number(a_dv(idir,:,:))
      a_dv(idir,:,:) = a_dv(idir,:,:) * 2.0d0 - 1.0d0
    end do

    alpha_orig = alpha
    alpha_dv_orig = alpha_dv
    a_orig = a
    a_dv_orig = a_dv
    x_orig = x
    x_dv_orig = x_dv
    y_orig = y
    y_dv_orig = y_dv

    write(*,*) 'Testing SGER (Vector Forward, n =', n, ')'

    call sger_dv(msize, nsize, alpha, alpha_dv, x, x_dv, incx_val, y, y_dv, incy_val, a, a_dv, lda_val, nbdirs)

    write(*,*) 'Function calls completed successfully'

    call check_derivatives_numerically(n, nbdirs, msize, nsize, lda_val, incx_val, incy_val, alpha_orig, alpha_dv_orig, x_orig, x_dv_orig, y_orig, y_dv_orig, a_orig, a_dv_orig, a_dv, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, nbdirs, msize, nsize, lda_val, incx_val, incy_val, alpha_orig, alpha_dv_orig, x_orig, x_dv_orig, y_orig, y_dv_orig, a_orig, a_dv_orig, a_dv, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    integer, intent(in) :: msize, nsize, lda_val, incx_val, incy_val
    real(4), intent(in) :: alpha_orig
    real(4), intent(in) :: alpha_dv_orig(nbdirs)
    real(4), intent(in) :: x_orig(n), x_dv_orig(nbdirs,n)
    real(4), intent(in) :: y_orig(n), y_dv_orig(nbdirs,n)
    real(4), intent(in) :: a_orig(n,n), a_dv_orig(nbdirs,n,n)
    real(4), intent(in) :: a_dv(nbdirs,n,n)
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3
    real(4) :: relative_error, max_error, abs_error, abs_reference, error_bound
    real(4) :: central_diff, ad_result
    logical :: has_large_errors
    real(4), dimension(n,n) :: a_forward, a_backward
    integer :: i, j, idir
    real(4) :: alpha
    real(4), dimension(n) :: x, y
    real(4), dimension(n,n) :: a

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    do idir = 1, nbdirs
      alpha = alpha_orig + h * alpha_dv_orig(idir)
      x = x_orig + h * x_dv_orig(idir,:)
      y = y_orig + h * y_dv_orig(idir,:)
      a = a_orig + h * a_dv_orig(idir,:,:)
      call sger(msize, nsize, alpha, x, incx_val, y, incy_val, a, lda_val)
      a_forward = a
      alpha = alpha_orig - h * alpha_dv_orig(idir)
      x = x_orig - h * x_dv_orig(idir,:)
      y = y_orig - h * y_dv_orig(idir,:)
      a = a_orig - h * a_dv_orig(idir,:,:)
      call sger(msize, nsize, alpha, x, incx_val, y, incy_val, a, lda_val)
      a_backward = a
      do j = 1, min(4, n)
        do i = 1, min(4, n)
          central_diff = (a_forward(i,j) - a_backward(i,j)) / (2.0e0 * h)
          ad_result = a_dv(idir,i,j)
          abs_error = abs(central_diff - ad_result)
          abs_reference = abs(ad_result)
          error_bound = 1.0e-3 + 1.0e-3 * abs_reference
          if (abs_error > error_bound) has_large_errors = .true.
          relative_error = abs_error / max(abs_reference, 1.0e-10)
          max_error = max(max_error, relative_error)
        end do
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

end program test_sger_vector_forward