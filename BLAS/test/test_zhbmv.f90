! Test program for ZHBMV differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size outlined run_test_for_size(n) - BLAS2 band (declarations in subroutines)

program test_zhbmv
  implicit none
  external :: zhbmv
  external :: zhbmv_d
  integer :: n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing ZHBMV (multi-size: n = 4)'
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
    complex(8) :: alpha, alpha_d, alpha_orig, alpha_d_seed
    complex(8) :: beta, beta_d, beta_orig, beta_d_seed
    complex(8), dimension(:,:), allocatable :: a, a_d, a_orig, a_d_seed
    complex(8), dimension(:), allocatable :: x, x_d, x_orig, x_d_seed
    complex(8), dimension(:), allocatable :: y, y_d, y_orig, y_d_seed
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
    allocate(y(n), y_d(n), y_orig(n), y_d_seed(n))
    ! Initialize a as Hermitian band matrix (upper band storage, real diagonal)
    do j = 1, n
    do band_row = max(1, ksize+2-j), ksize+1
    if (band_row .eq. ksize+1) then
    call random_number(temp_real)
    a(band_row, j) = cmplx(temp_real * 2.0 - 1.0, 0.0)  ! Real diagonal
    else
    call random_number(temp_real)
    call random_number(temp_imag)
    a(band_row, j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end if
    end do
    end do
    ! Keep direction consistent with Hermitian band: real diagonal, band entries only
    do j = 1, n
    do band_row = max(1, ksize+2-j), ksize+1
    if (band_row .eq. ksize+1) then
    call random_number(temp_real)
    a_d(band_row, j) = cmplx(temp_real * 2.0 - 1.0, 0.0d0)
    else
    call random_number(temp_real)
    call random_number(temp_imag)
    a_d(band_row, j) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
    end if
    end do
    end do
    call random_number(temp_real)
    call random_number(temp_imag)
    alpha = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(alpha))
    call random_number(temp_real)
    call random_number(temp_imag)
    alpha_d = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(alpha_d))
    call random_number(temp_real)
    call random_number(temp_imag)
    beta = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(beta))
    call random_number(temp_real)
    call random_number(temp_imag)
    beta_d = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(beta_d))
    do j = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      x(j) = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(x))
      x_d(j) = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(x_d))
    end do
    do j = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      y(j) = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(y))
      y_d(j) = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(y_d))
    end do
    write(*,*) 'Testing ZHBMV (n =', n, ')'
    a_orig = a
    a_d_seed = a_d
    x_orig = x
    x_d_seed = x_d
    alpha_orig = alpha
    alpha_d_seed = alpha_d
    y_orig = y
    y_d_seed = y_d
    beta_orig = beta
    beta_d_seed = beta_d
    call zhbmv_d(uplo, nsize, ksize, alpha, alpha_d, a, a_d, lda_val, x, x_d, incx_val, beta, beta_d, y, y_d, incy_val)
    call check_derivatives_numerically_band(n, lda_val, ksize, uplo, nsize, incx_val, incy_val, alpha_orig, alpha_d_seed, beta_orig, beta_d_seed, a_orig, a_d_seed, x_orig, x_d_seed, y_orig, y_d_seed, y_d, passed)
    deallocate(a, a_d, a_orig, a_d_seed, x, x_d, x_orig, x_d_seed)
    deallocate(y, y_d, y_orig, y_d_seed)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically_band(n, lda_val, ksize, uplo, nsize, incx_val, incy_val, alpha_orig, alpha_d_seed, beta_orig, beta_d_seed, a_orig, a_d_seed, x_orig, x_d_seed, y_orig, y_d_seed, y_d_out, passed)
    implicit none
    integer, intent(in) :: n, lda_val, ksize, nsize, incx_val, incy_val
    character, intent(in) :: uplo
    complex(8), intent(in) :: alpha_orig, alpha_d_seed, beta_orig, beta_d_seed
    complex(8), intent(in) :: a_orig(lda_val, n), a_d_seed(lda_val, n), x_orig(n), x_d_seed(n), y_orig(n), y_d_seed(n), y_d_out(n)
    logical, intent(out) :: passed
    real(8), parameter :: h = 1.0e-7
    real(8) :: abs_error, abs_ref, err_bound
    complex(8), dimension(n) :: y_fwd, y_bwd, y_t
    complex(8) :: alpha_t, beta_t
    complex(8), dimension(n) :: x_t
    complex(8), dimension(lda_val, n) :: a_t
    integer :: ii
    logical :: has_err
    has_err = .false.
    alpha_t = alpha_orig + h * alpha_d_seed
    beta_t = beta_orig + h * beta_d_seed
    a_t = a_orig + h * a_d_seed
    x_t = x_orig + h * x_d_seed
    y_t = y_orig + h * y_d_seed
    call zhbmv(uplo, nsize, ksize, alpha_t, a_t, lda_val, x_t, incx_val, beta_t, y_t, incy_val)
    y_fwd = y_t
    alpha_t = alpha_orig - h * alpha_d_seed
    beta_t = beta_orig - h * beta_d_seed
    a_t = a_orig - h * a_d_seed
    x_t = x_orig - h * x_d_seed
    y_t = y_orig - h * y_d_seed
    call zhbmv(uplo, nsize, ksize, alpha_t, a_t, lda_val, x_t, incx_val, beta_t, y_t, incy_val)
    y_bwd = y_t
    do ii = 1, min(3, n)
      abs_error = abs((y_fwd(ii) - y_bwd(ii)) / (2.0e0 * h) - y_d_out(ii))
      abs_ref = abs(y_d_out(ii))
      err_bound = 1.0e-5 + 1.0e-5 * abs_ref
      if (abs_error > err_bound) has_err = .true.
    end do
    passed = .not. has_err
    if (has_err) write(*,*) 'FAIL: Band scalar derivatives'
    if (.not. has_err) write(*,*) 'PASS: Band scalar derivatives'
  end subroutine check_derivatives_numerically_band
end program test_zhbmv