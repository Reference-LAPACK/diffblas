! Test program for SGEMV vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=n
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_sgemv_vector_forward
  implicit none

  external :: sgemv
  external :: sgemv_dv

  integer :: nbdirs
  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing SGEMV (Vector Forward, multi-size: n = 4)'
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

    character :: trans
    integer :: msize, nsize, lda_val, incx_val, incy_val
    real(4) :: alpha, beta
    real(4), dimension(n,n) :: a
    real(4), dimension(n) :: x, y
    real(4), dimension(nbdirs) :: alpha_dv, beta_dv
    real(4), dimension(nbdirs,n,n) :: a_dv
    real(4), dimension(nbdirs,n) :: x_dv, y_dv
    real(4) :: alpha_orig, beta_orig
    real(4), dimension(nbdirs) :: alpha_dv_orig, beta_dv_orig
    real(4), dimension(n,n) :: a_orig
    real(4), dimension(nbdirs,n,n) :: a_dv_orig
    real(4), dimension(n) :: x_orig, y_orig
    real(4), dimension(nbdirs,n) :: x_dv_orig, y_dv_orig
    integer :: idir, ii, jj
    real(4) :: temp_real, temp_imag

    trans = 'N'
    msize = n
    nsize = n
    lda_val = n
    incx_val = 1
    incy_val = 1

    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(a)
    a = a * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    call random_number(beta)
    beta = beta * 2.0d0 - 1.0d0
    call random_number(y)
    y = y * 2.0d0 - 1.0d0

    do idir = 1, nbdirs
      call random_number(temp_real)
      alpha_dv(idir) = temp_real * 2.0d0 - 1.0d0
    end do
    do idir = 1, nbdirs
      call random_number(a_dv(idir,:,:))
      a_dv(idir,:,:) = a_dv(idir,:,:) * 2.0d0 - 1.0d0
    end do
    do idir = 1, nbdirs
      call random_number(x_dv(idir,:))
      x_dv(idir,:) = x_dv(idir,:) * 2.0d0 - 1.0d0
    end do
    do idir = 1, nbdirs
      call random_number(temp_real)
      beta_dv(idir) = temp_real * 2.0d0 - 1.0d0
    end do
    do idir = 1, nbdirs
      call random_number(y_dv(idir,:))
      y_dv(idir,:) = y_dv(idir,:) * 2.0d0 - 1.0d0
    end do

    alpha_orig = alpha
    alpha_dv_orig = alpha_dv
    a_orig = a
    a_dv_orig = a_dv
    x_orig = x
    x_dv_orig = x_dv
    beta_orig = beta
    beta_dv_orig = beta_dv
    y_orig = y
    y_dv_orig = y_dv

    write(*,*) 'Testing SGEMV (Vector Forward, n =', n, ')'

    call sgemv_dv(trans, msize, nsize, alpha, alpha_dv, a, a_dv, lda_val, x, x_dv, incx_val, beta, beta_dv, y, y_dv, incy_val, nbdirs)

    write(*,*) 'Function calls completed successfully'

    call check_derivatives_numerically(n, nbdirs, trans, msize, nsize, lda_val, incx_val, incy_val, alpha_orig, alpha_dv_orig, a_orig, a_dv_orig, x_orig, x_dv_orig, beta_orig, beta_dv_orig, y_orig, y_dv_orig, y_dv, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, nbdirs, trans, msize, nsize, lda_val, incx_val, incy_val, alpha_orig, alpha_dv_orig, a_orig, a_dv_orig, x_orig, x_dv_orig, beta_orig, beta_dv_orig, y_orig, y_dv_orig, y_dv, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    character, intent(in) :: trans
    integer, intent(in) :: msize, nsize, lda_val, incx_val, incy_val
    real(4), intent(in) :: alpha_orig, beta_orig
    real(4), intent(in) :: alpha_dv_orig(nbdirs), beta_dv_orig(nbdirs)
    real(4), intent(in) :: a_orig(n,n), a_dv_orig(nbdirs,n,n)
    real(4), intent(in) :: x_orig(n), x_dv_orig(nbdirs,n)
    real(4), intent(in) :: y_orig(n), y_dv_orig(nbdirs,n)
    real(4), intent(in) :: y_dv(nbdirs,n)
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3
    real(4) :: relative_error, max_error, abs_error, abs_reference, error_bound
    real(4) :: central_diff, ad_result
    logical :: has_large_errors
    real(4), dimension(n) :: y_forward, y_backward
    integer :: i, idir
    real(4) :: alpha, beta
    real(4), dimension(n,n) :: a
    real(4), dimension(n) :: x, y

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    do idir = 1, nbdirs
      alpha = alpha_orig + h * alpha_dv_orig(idir)
      a = a_orig + h * a_dv_orig(idir,:,:)
      x = x_orig + h * x_dv_orig(idir,:)
      beta = beta_orig + h * beta_dv_orig(idir)
      y = y_orig + h * y_dv_orig(idir,:)
      call sgemv(trans, msize, nsize, alpha, a, lda_val, x, incx_val, beta, y, incy_val)
      y_forward = y
      alpha = alpha_orig - h * alpha_dv_orig(idir)
      a = a_orig - h * a_dv_orig(idir,:,:)
      x = x_orig - h * x_dv_orig(idir,:)
      beta = beta_orig - h * beta_dv_orig(idir)
      y = y_orig - h * y_dv_orig(idir,:)
      call sgemv(trans, msize, nsize, alpha, a, lda_val, x, incx_val, beta, y, incy_val)
      y_backward = y
      do i = 1, min(4, n)
        central_diff = (y_forward(i) - y_backward(i)) / (2.0e0 * h)
        ad_result = y_dv(idir,i)
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 1.0e-3 + 1.0e-3 * abs_reference
        if (abs_error > error_bound) has_large_errors = .true.
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        max_error = max(max_error, relative_error)
      end do
    end do

    write(*,*) 'Maximum relative error across all directions:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in vector derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Vector derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_sgemv_vector_forward