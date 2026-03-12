! Test program for CTPMV differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size outlined run_test_for_size(n) - TPMV/TPSV packed triangular

program test_ctpmv
  implicit none
  external :: ctpmv
  external :: ctpmv_d
  integer :: n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing CTPMV (multi-size: n = 4)'
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
    complex(4), allocatable :: ap(:), ap_d(:), x(:), x_d(:)
    complex(4), allocatable :: ap_t(:), x_t(:), x_plus(:), x_minus(:)
    complex(4), allocatable :: ap_d_seed(:), x_d_seed(:)
    complex(4), allocatable :: ap_orig(:), x_orig(:)
    integer :: ii
    real(4) :: tr, ti
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
    do ii = 1, npack
      call random_number(tr)
      call random_number(ti)
      ap(ii) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(ap))
    end do
    do ii = 1, n
      call random_number(tr)
      call random_number(ti)
      x(ii) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(x))
    end do
    do ii = 1, npack
      call random_number(tr)
      call random_number(ti)
      ap_d(ii) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(ap_d))
    end do
    do ii = 1, n
      call random_number(tr)
      call random_number(ti)
      x_d(ii) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(x_d))
    end do
    ap_orig = ap
    x_orig = x
    ap_d_seed = ap_d
    x_d_seed = x_d
    call ctpmv_d(uplo, trans, diag, nsize, ap, ap_d, x, x_d, incx_val)
    call check_derivatives_numerically(n, npack, uplo, trans, diag, nsize, incx_val, ap_orig, ap_d_seed, x_orig, x_d_seed, x_d, passed)
    deallocate(ap, ap_d, x, x_d, ap_t, x_t, x_plus, x_minus, ap_d_seed, x_d_seed, ap_orig, x_orig)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, npack, uplo, trans, diag, nsize, incx_val, ap, ap_d_seed, x, x_d_seed, x_d, passed)
    implicit none
    integer, intent(in) :: n, npack, nsize, incx_val
    character, intent(in) :: uplo, trans, diag
    complex(4), intent(in) :: ap(npack), ap_d_seed(npack), x(n), x_d_seed(n), x_d(n)
    logical, intent(out) :: passed
    real(4), parameter :: h = 1.0e-3
    complex(4) :: ap_t(npack), x_t(n), x_plus(n), x_minus(n)
    complex(4) :: central_diff, ad_result
    logical :: has_err
    integer :: ii
    real(4) :: abs_error, abs_ref, err_bound, relative_error, max_error
    has_err = .false.
    max_error = 0.0e0
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    ap_t = ap + h * ap_d_seed
    x_t = x + h * x_d_seed
    call ctpmv(uplo, trans, diag, nsize, ap_t, x_t, incx_val)
    x_plus = x_t
    ap_t = ap - h * ap_d_seed
    x_t = x - h * x_d_seed
    call ctpmv(uplo, trans, diag, nsize, ap_t, x_t, incx_val)
    x_minus = x_t
    do ii = 1, min(2, n)
      central_diff = (x_plus(ii) - x_minus(ii)) / (2.0e0 * h)
      ad_result = x_d(ii)
      abs_error = abs(central_diff - ad_result)
      abs_ref = abs(ad_result)
      err_bound = 1.0e-3 + 1.0e-3 * abs_ref
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
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    passed = .not. has_err
    if (has_err) write(*,*) 'FAIL: TPMV/TPSV scalar derivatives'
    if (.not. has_err) write(*,*) 'PASS: TPMV/TPSV scalar derivatives'
  end subroutine check_derivatives_numerically
end program test_ctpmv