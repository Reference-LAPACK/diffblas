! Test program for SSCAL vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=n
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_sscal_vector_reverse
  implicit none

  external :: sscal
  external :: sscal_bv

  integer :: nbdirs
  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(3)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing SSCAL (Vector Reverse, multi-size: n =', test_sizes(1), ')'
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
    real(4) :: alpha
    real(4), dimension(n) :: x
    real(4), dimension(nbdirs) :: alphab
    real(4), dimension(nbdirs,n) :: xb
    real(4) :: alpha_orig
    real(4), dimension(n) :: x_orig
    real(4), dimension(nbdirs,n) :: xb_orig
    integer :: k, i
    real(4) :: temp_real, temp_imag

    nsize = n
    incx_val = 1

    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0

    alpha_orig = alpha
    x_orig = x

    do k = 1, nbdirs
      call random_number(xb(k,:))
      xb(k,:) = xb(k,:) * 2.0d0 - 1.0d0
    end do
    xb_orig = xb

    alphab = 0.0d0

    write(*,*) 'Testing SSCAL (Vector Reverse, n =', n, ')'

    call sscal_bv(nsize, alpha, alphab, x, xb, incx_val, nbdirs)

    call check_vjp_numerically(n, nbdirs, nsize, incx_val, alpha_orig, x_orig, xb_orig, alphab, xb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nbdirs, nsize, incx_val, alpha_orig, x_orig, xb_orig, alphab, xb, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    integer, intent(in) :: nsize, incx_val
    real(4), intent(in) :: alpha_orig
    real(4), intent(in) :: x_orig(n)
    real(4), intent(in) :: xb_orig(nbdirs,n)
    real(4), intent(in) :: alphab(nbdirs)
    real(4), intent(in) :: xb(nbdirs,n)
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    real(4) :: alpha_dir
    real(4), dimension(n) :: x_dir
    real(4) :: alpha
    real(4), dimension(n) :: x, x_plus, x_minus, x_central_diff
    real(4), dimension(n) :: temp_products
    integer :: i, k
    real(4) :: temp_real, temp_imag
    logical :: has_large_errors

    max_error = 0.0d0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    do k = 1, nbdirs
      call random_number(alpha_dir)
      alpha_dir = alpha_dir * 2.0d0 - 1.0d0
      call random_number(x_dir)
      x_dir = x_dir * 2.0d0 - 1.0d0
      alpha = alpha_orig + h * alpha_dir
      x = x_orig + h * x_dir
      call sscal(nsize, alpha, x, incx_val)
      x_plus = x
      alpha = alpha_orig - h * alpha_dir
      x = x_orig - h * x_dir
      call sscal(nsize, alpha, x, incx_val)
      x_minus = x
      x_central_diff = (x_plus - x_minus) / (2.0d0 * h)
      vjp_fd = 0.0d0
      do i = 1, n
        temp_products(i) = xb_orig(k,i) * x_central_diff(i)
        vjp_fd = vjp_fd + temp_products(i)
      end do
      vjp_ad = 0.0d0
      vjp_ad = vjp_ad + alpha_dir * alphab(k)
      do i = 1, n
        vjp_ad = vjp_ad + x_dir(i) * xb(k,i)
      end do
      abs_error = abs(vjp_fd - vjp_ad)
      abs_reference = abs(vjp_ad)
      error_bound = 2.0e-3 + 2.0e-3 * abs_reference
      if (abs_error > error_bound) has_large_errors = .true.
      if (abs_reference > 1.0e-10) then
        relative_error = abs_error / abs_reference
      else
        relative_error = abs_error
      end if
      if (relative_error > max_error) max_error = relative_error
    end do
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_vjp_numerically

end program test_sscal_vector_reverse