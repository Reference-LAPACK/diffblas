! Test program for DSPR differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size outlined run_test_for_size(n) - SPR/SPR2 packed (declarations in subroutines)

program test_dspr
  implicit none
  external :: dspr
  external :: dspr_d
  integer :: n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing DSPR (multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 1
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
    real(8) :: alpha, alpha_d
    real(8), dimension(n) :: x, x_d
    real(8), allocatable :: ap(:), ap_d(:), ap_d_seed(:), ap_orig(:)
    integer :: ii
    uplo = 'U'
    nsize = n
    incx_val = 1
    incy_val = 1
    npack = (n * (n + 1)) / 2
    allocate(ap(npack), ap_d(npack), ap_d_seed(npack), ap_orig(npack))
    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    call random_number(ap)
    ap = ap * 2.0d0 - 1.0d0
    call random_number(alpha_d)
    alpha_d = alpha_d * 2.0d0 - 1.0d0
    call random_number(x_d)
    x_d = x_d * 2.0d0 - 1.0d0
    call random_number(ap_d)
    ap_d = ap_d * 2.0d0 - 1.0d0
    ap_orig = ap
    ap_d_seed = ap_d
    write(*,*) 'Testing DSPR (n =', n, ')'
    call dspr_d(uplo, nsize, alpha, alpha_d, x, x_d, incx_val, ap, ap_d)
    call check_derivatives_numerically(n, npack, uplo, nsize, incx_val, alpha, alpha_d, x, x_d, ap_orig, ap_d, ap_d_seed, passed)
    deallocate(ap, ap_d, ap_d_seed, ap_orig)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, npack, uplo, nsize, incx_val, alpha, alpha_d, x, x_d, ap_orig, ap_d, ap_d_seed, passed)
    implicit none
    integer, intent(in) :: n, npack
    character, intent(in) :: uplo
    integer, intent(in) :: nsize, incx_val
    real(8), intent(in) :: alpha, alpha_d
    real(8), intent(in) :: x(n), x_d(n)
    real(8), intent(in) :: ap_orig(npack), ap_d(npack), ap_d_seed(npack)
    logical, intent(out) :: passed
    real(8), parameter :: h = 1.0e-7
    real(8) :: abs_error, abs_ref, err_bound, max_error, relative_error
    real(8), dimension(npack) :: ap_fwd, ap_bwd, ap_t
    real(8) :: alpha_t
    real(8), dimension(n) :: x_t
    integer :: ii
    logical :: has_err
    has_err = .false.
    alpha_t = alpha + h * alpha_d
    x_t = x + h * x_d
    ap_t = ap_orig + h * ap_d_seed
    call dspr(uplo, nsize, alpha_t, x_t, incx_val, ap_t)
    ap_fwd = ap_t
    alpha_t = alpha - h * alpha_d
    x_t = x - h * x_d
    ap_t = ap_orig - h * ap_d_seed
    call dspr(uplo, nsize, alpha_t, x_t, incx_val, ap_t)
    ap_bwd = ap_t
    has_err = .false.
    max_error = 0.0e0
    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    do ii = 1, npack
      abs_error = abs((ap_fwd(ii) - ap_bwd(ii)) / (2.0e0 * h) - ap_d(ii))
      abs_ref = abs(ap_d(ii))
      err_bound = 1.0e-5 + 1.0e-5 * abs_ref
      if (abs_error > max_error) max_error = abs_error
      if (abs_error > err_bound) has_err = .true.
    end do
    relative_error = 0.0e0
    abs_ref = maxval(abs(ap_d)) + 1.0e0
    if (abs_ref > 1.0e-10) relative_error = max_error / abs_ref
    write(*,*) 'Maximum relative error:', relative_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_err
    if (has_err) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
  end subroutine check_derivatives_numerically
end program test_dspr