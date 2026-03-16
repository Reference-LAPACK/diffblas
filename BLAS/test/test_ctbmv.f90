! Test program for CTBMV differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size outlined run_test_for_size(n) - BLAS2 band (declarations in subroutines)

program test_ctbmv
  implicit none
  external :: ctbmv
  external :: ctbmv_d
  integer :: n_test, seed_array(33), test_sizes(3), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing CTBMV (multi-size: n = 4)'
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
    character :: uplo, trans, diag
    integer :: nsize, ksize, lda_val, incx_val, incy_val
    complex(4) :: alpha, alpha_d, alpha_orig, alpha_d_seed
    complex(4), dimension(:,:), allocatable :: a, a_d, a_orig, a_d_seed
    complex(4), dimension(:), allocatable :: x, x_d, x_orig, x_d_seed
    integer :: band_row, j
    real(4) :: temp_real, temp_imag
    ksize = max(0, n - 1)
    nsize = n
    lda_val = ksize + 1
    incx_val = 1
    incy_val = 1
    uplo = 'U'
    trans = 'N'
    diag = 'N'
    allocate(a(lda_val, n), a_d(lda_val, n), a_orig(lda_val, n), a_d_seed(lda_val, n))
    allocate(x(n), x_d(n), x_orig(n), x_d_seed(n))
    ! Initialize a as triangular band matrix (upper band storage)
    do j = 1, n
    do band_row = max(1, ksize+2-j), ksize+1
    call random_number(temp_real)
    call random_number(temp_imag)
    a(band_row, j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
    end do
    ! Keep direction consistent with triangular band: only band entries used
    do j = 1, n
    do band_row = max(1, ksize+2-j), ksize+1
    call random_number(temp_real)
    call random_number(temp_imag)
    a_d(band_row, j) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
    end do
    end do
    call random_number(temp_real)
    call random_number(temp_imag)
    alpha = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(alpha))
    call random_number(temp_real)
    call random_number(temp_imag)
    alpha_d = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(alpha_d))
    do j = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      x(j) = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(x))
      x_d(j) = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(x_d))
    end do
    write(*,*) 'Testing CTBMV (n =', n, ')'
    a_orig = a
    a_d_seed = a_d
    x_orig = x
    x_d_seed = x_d
    alpha_orig = alpha
    alpha_d_seed = alpha_d
    call ctbmv_d(uplo, trans, diag, nsize, ksize, a, a_d, lda_val, x, x_d, incx_val)
    ! Reset input derivative vars from seeds; output derivative (x_d or y_d) keeps AD result
    a_d = a_d_seed
    alpha_d = alpha_d_seed
    write(*,*) 'Function calls completed successfully'
    call check_derivatives_numerically_band(n, lda_val, ksize, uplo, trans, diag, nsize, incx_val, a_orig, a_d_seed, x_orig, x_d_seed, x_d, passed)
    deallocate(a, a_d, a_orig, a_d_seed, x, x_d, x_orig, x_d_seed)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically_band(n, lda_val, ksize, uplo, trans, diag, nsize, incx_val, a_orig, a_d_seed, x_orig, x_d_seed, x_d_out, passed)
    implicit none
    integer, intent(in) :: n, lda_val, ksize, nsize, incx_val
    character, intent(in) :: uplo, trans, diag
    complex(4), intent(in) :: a_orig(lda_val, n), a_d_seed(lda_val, n), x_orig(n), x_d_seed(n), x_d_out(n)
    logical, intent(out) :: passed
    real(4), parameter :: h = 1.0e-3
    real(4) :: abs_error, abs_ref, err_bound, max_error, relative_error
    complex(4), dimension(n) :: x_fwd, x_bwd, x_t
    complex(4), dimension(lda_val, n) :: a_t
    integer :: ii, j, band_row
    logical :: has_err
    has_err = .false.
    max_error = 0.0e0
    a_t = a_orig
    do j = 1, n
      do band_row = max(1, ksize+2-j), ksize+1
        a_t(band_row, j) = a_orig(band_row, j) + h * a_d_seed(band_row, j)
      end do
    end do
    x_t = x_orig + h * x_d_seed
    call ctbmv(uplo, trans, diag, nsize, ksize, a_t, lda_val, x_t, incx_val)
    x_fwd = x_t
    a_t = a_orig
    do j = 1, n
      do band_row = max(1, ksize+2-j), ksize+1
        a_t(band_row, j) = a_orig(band_row, j) - h * a_d_seed(band_row, j)
      end do
    end do
    x_t = x_orig - h * x_d_seed
    call ctbmv(uplo, trans, diag, nsize, ksize, a_t, lda_val, x_t, incx_val)
    x_bwd = x_t
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    do ii = 1, n
      abs_error = abs((x_fwd(ii) - x_bwd(ii)) / (2.0e0 * h) - x_d_out(ii))
      abs_ref = abs(x_d_out(ii))
      err_bound = 1.0e-3 + 1.0e-3 * abs_ref
      if (abs_error > err_bound) has_err = .true.
      relative_error = abs_error / max(abs_ref, 1.0e-10)
      if (relative_error > max_error) max_error = relative_error
    end do
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    passed = .not. has_err
    if (has_err) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
  end subroutine check_derivatives_numerically_band
end program test_ctbmv