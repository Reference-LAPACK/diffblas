! Test program for DTBMV differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size outlined run_test_for_size(n) - BLAS2 band (declarations in subroutines)

program test_dtbmv
  implicit none
  external :: dtbmv
  external :: dtbmv_d
  integer :: n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing DTBMV (multi-size: n = 4)'
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
    integer :: nsize, ksize, lda_val, incx_val, incy_val
    real(8) :: alpha, alpha_d, alpha_orig, alpha_d_seed
    real(8), dimension(:,:), allocatable :: a, a_d, a_orig, a_d_seed
    real(8), dimension(:), allocatable :: x, x_d, x_orig, x_d_seed
    integer :: band_row, j
    real(4) :: temp_real
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
    ! A(band_row, j) = full(i,j) with band_row = ksize+1+i-j, i = max(1,j-ksize)..j
    do j = 1, n
    do band_row = max(1, ksize+2-j), ksize+1
    call random_number(temp_real)
    a(band_row, j) = temp_real * 2.0 - 1.0  ! Scale to [-1,1]
    end do
    end do
    ! Keep direction consistent with triangular band: only band entries used
    do j = 1, n
    do band_row = max(1, ksize+2-j), ksize+1
    call random_number(temp_real)
    a_d(band_row, j) = temp_real * 2.0 - 1.0
    end do
    end do
    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(alpha_d)
    alpha_d = alpha_d * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    call random_number(x_d)
    x_d = x_d * 2.0d0 - 1.0d0
    write(*,*) 'Testing DTBMV (n =', n, ')'
    a_orig = a
    a_d_seed = a_d
    x_orig = x
    x_d_seed = x_d
    alpha_orig = alpha
    alpha_d_seed = alpha_d
    call dtbmv_d(uplo, trans, diag, nsize, ksize, a, a_d, lda_val, x, x_d, incx_val)
    write(*,*) 'Function calls completed successfully'
    call check_derivatives_numerically_band(n, lda_val, ksize, uplo, trans, diag, nsize, incx_val, a_orig, a_d_seed, x_orig, x_d_seed, x_d, passed)
    deallocate(a, a_d, a_orig, a_d_seed, x, x_d, x_orig, x_d_seed)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically_band(n, lda_val, ksize, uplo, trans, diag, nsize, incx_val, a_orig, a_d_seed, x_orig, x_d_seed, x_d_out, passed)
    implicit none
    integer, intent(in) :: n, lda_val, ksize, nsize, incx_val
    character, intent(in) :: uplo, trans, diag
    real(8), intent(in) :: a_orig(lda_val, n), a_d_seed(lda_val, n), x_orig(n), x_d_seed(n), x_d_out(n)
    logical, intent(out) :: passed
    real(8), parameter :: h = 1.0e-7
    real(8) :: abs_error, abs_ref, err_bound, max_error, relative_error
    real(8), dimension(n) :: x_fwd, x_bwd, x_t
    real(8), dimension(lda_val, n) :: a_t
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
    call dtbmv(uplo, trans, diag, nsize, ksize, a_t, lda_val, x_t, incx_val)
    x_fwd = x_t
    a_t = a_orig
    do j = 1, n
      do band_row = max(1, ksize+2-j), ksize+1
        a_t(band_row, j) = a_orig(band_row, j) - h * a_d_seed(band_row, j)
      end do
    end do
    x_t = x_orig - h * x_d_seed
    call dtbmv(uplo, trans, diag, nsize, ksize, a_t, lda_val, x_t, incx_val)
    x_bwd = x_t
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    do ii = 1, n
      abs_error = abs((x_fwd(ii) - x_bwd(ii)) / (2.0e0 * h) - x_d_out(ii))
      abs_ref = abs(x_d_out(ii))
      err_bound = 1.0e-5 + 1.0e-5 * abs_ref
      if (abs_error > err_bound) has_err = .true.
      relative_error = abs_error / max(abs_ref, 1.0e-10)
      if (relative_error > max_error) max_error = relative_error
    end do
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_err
    if (has_err) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
  end subroutine check_derivatives_numerically_band
end program test_dtbmv