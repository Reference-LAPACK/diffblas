! Test program for SSPMV differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size outlined run_test_for_size(n) - SPMV (symmetric packed matrix-vector)

program test_sspmv
  implicit none
  external :: sspmv
  external :: sspmv_d
  integer :: n_test, seed_array(33), test_sizes(3), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing SSPMV (multi-size: n =', test_sizes(1), ')'
  all_passed = .true.
  do i = 1, 3
    n_test = test_sizes(i)
    call run_test_for_size(n_test, passed)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) then
    write(*,*) 'PASS: All sizes completed successfully'
  else
    write(*,*) 'FAIL: One or more sizes had derivative errors'
  end if
contains
  subroutine run_test_for_size(n, passed)
    implicit none
    integer, intent(in) :: n
    logical, intent(out) :: passed
    character :: uplo
    integer :: nsize, incx_val, incy_val, npack
    real(4) :: alpha, alpha_d, beta, beta_d
    real(4), dimension(n) :: x, x_d, y, y_d, y_d_seed, y_orig, y_plus, y_minus
    real(4), dimension(:), allocatable :: ap, ap_d, ap_t, ap_orig
    real(4) :: alpha_t, beta_t
    real(4), dimension(n) :: x_t
    real(4) :: h
    parameter (h = 1.0e-3)
    real(4) :: abs_error, abs_ref, err_bound, max_err
    integer :: ii
    uplo = 'U'
    nsize = n
    incx_val = 1
    incy_val = 1
    npack = (n * (n + 1)) / 2
    allocate(ap(npack), ap_d(npack), ap_t(npack), ap_orig(npack))
    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(beta)
    beta = beta * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    call random_number(y)
    y = y * 2.0d0 - 1.0d0
    call random_number(ap)
    ap = ap * 2.0d0 - 1.0d0
    call random_number(alpha_d)
    alpha_d = alpha_d * 2.0d0 - 1.0d0
    call random_number(beta_d)
    beta_d = beta_d * 2.0d0 - 1.0d0
    call random_number(x_d)
    x_d = x_d * 2.0d0 - 1.0d0
    call random_number(y_d)
    y_d = y_d * 2.0d0 - 1.0d0
    call random_number(ap_d)
    ap_d = ap_d * 2.0d0 - 1.0d0
    ap_orig = ap
    y_orig = y
    y_d_seed = y_d
    write(*,*) 'Testing SSPMV (n =', n, ')'
    call sspmv_d(uplo, nsize, alpha, alpha_d, ap, ap_d, x, x_d, incx_val, beta, beta_d, y, y_d, incy_val)
    ! FD check: perturb all inputs and inout y by directions (y_d_seed for inout y); use ap_orig for base
    alpha_t = alpha + h * alpha_d
    beta_t = beta + h * beta_d
    x_t = x + h * x_d
    y_plus = y_orig + h * y_d_seed
    ap_t = ap_orig + h * ap_d
    call sspmv(uplo, nsize, alpha_t, ap_t, x_t, incx_val, beta_t, y_plus, incy_val)
    alpha_t = alpha - h * alpha_d
    beta_t = beta - h * beta_d
    x_t = x - h * x_d
    y_minus = y_orig - h * y_d_seed
    ap_t = ap_orig - h * ap_d
    call sspmv(uplo, nsize, alpha_t, ap_t, x_t, incx_val, beta_t, y_minus, incy_val)
    max_err = 0.0d0
    do ii = 1, n
      abs_error = abs((y_plus(ii) - y_minus(ii)) / (2.0d0 * h) - y_d(ii))
      if (abs_error > max_err) max_err = abs_error
    end do
    abs_ref = maxval(abs(y_d)) + 1.0d0
    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Maximum relative error:', max_err / abs_ref
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
    passed = (max_err <= 2.0e-3 * abs_ref)
    if (.not. passed) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
    deallocate(ap, ap_d, ap_t, ap_orig)
  end subroutine run_test_for_size
end program test_sspmv