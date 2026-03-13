! Test program for DTPMV differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size outlined run_test_for_size(n) - TPMV/TPSV packed triangular

program test_dtpmv
  implicit none
  external :: dtpmv
  external :: dtpmv_d
  integer :: n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing DTPMV (multi-size: n = 4)'
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
    character :: uplo, trans, diag
    integer :: nsize, incx_val, npack
    real(8), allocatable :: ap(:), ap_d(:), x(:), x_d(:)
    real(8), allocatable :: ap_t(:), x_t(:), x_plus(:), x_minus(:)
    real(8), allocatable :: ap_d_seed(:), x_d_seed(:)
    real(8), allocatable :: ap_orig(:), x_orig(:)
    integer :: ii
    uplo = 'U'
    trans = 'N'
    diag = 'N'
    nsize = n
    incx_val = 1
    npack = (n * (n + 1)) / 2
    allocate(ap(npack), ap_d(npack), x(n), x_d(n))
    allocate(ap_t(npack), x_t(n), x_plus(n), x_minus(n))
    allocate(ap_d_seed(npack), x_d_seed(n))
    allocate(ap_orig(npack), x_orig(n))
    call random_number(ap)
    ap = ap * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    call random_number(ap_d)
    ap_d = ap_d * 2.0d0 - 1.0d0
    call random_number(x_d)
    x_d = x_d * 2.0d0 - 1.0d0
    ap_orig = ap
    x_orig = x
    ap_d_seed = ap_d
    x_d_seed = x_d
    call dtpmv_d(uplo, trans, diag, nsize, ap, ap_d, x, x_d, incx_val)
    write(*,*) 'Testing DTPMV (n =', n, ')'
    write(*,*) 'Function calls completed successfully'
    call check_derivatives_numerically(n, npack, uplo, trans, diag, nsize, incx_val, ap_orig, ap_d_seed, x_orig, x_d_seed, x_d, passed)
    deallocate(ap, ap_d, x, x_d, ap_t, x_t, x_plus, x_minus, ap_d_seed, x_d_seed, ap_orig, x_orig)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, npack, uplo, trans, diag, nsize, incx_val, ap, ap_d_seed, x, x_d_seed, x_d, passed)
    implicit none
    integer, intent(in) :: n, npack, nsize, incx_val
    character, intent(in) :: uplo, trans, diag
    real(8), intent(in) :: ap(npack), ap_d_seed(npack), x(n), x_d_seed(n), x_d(n)
    logical, intent(out) :: passed
    real(8), parameter :: h = 1.0e-7
    real(8) :: ap_t(npack), x_t(n), x_plus(n), x_minus(n)
    real(8) :: central_diff, ad_result
    logical :: has_err
    integer :: ii
    real(8) :: abs_error, abs_ref, err_bound, relative_error, max_error
    has_err = .false.
    max_error = 0.0d0
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    ap_t = ap + h * ap_d_seed
    x_t = x + h * x_d_seed
    call dtpmv(uplo, trans, diag, nsize, ap_t, x_t, incx_val)
    x_plus = x_t
    ap_t = ap - h * ap_d_seed
    x_t = x - h * x_d_seed
    call dtpmv(uplo, trans, diag, nsize, ap_t, x_t, incx_val)
    x_minus = x_t
    do ii = 1, min(2, n)
      central_diff = (x_plus(ii) - x_minus(ii)) / (2.0d0 * h)
      ad_result = x_d(ii)
      abs_error = abs(central_diff - ad_result)
      abs_ref = abs(ad_result)
      err_bound = 1.0e-5 + 1.0e-5 * abs_ref
      if (abs_error > err_bound) then
        has_err = .true.
        relative_error = abs_error / max(abs_ref, 1.0e-10)
        write(*,*) 'Large error in output X(', ii, '):'
        write(*,*) '  Central diff: ', central_diff
        write(*,*) '  AD result:   ', ad_result
        write(*,*) '  Absolute error:', abs_error
        write(*,*) '  Error bound:', err_bound
        write(*,*) '  Relative error:', relative_error
      end if
      relative_error = abs_error / max(abs_ref, 1.0e-10)
      max_error = max(max_error, relative_error)
    end do
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_err
    if (has_err) write(*,*) 'FAIL: Derivatives are outside tolerance'
    if (.not. has_err) write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
  end subroutine check_derivatives_numerically
end program test_dtpmv