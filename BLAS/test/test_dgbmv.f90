! Test program for DGBMV differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size outlined run_test_for_size(n) - BLAS2 band (declarations in subroutines)

program test_dgbmv
  implicit none
  external :: dgbmv
  external :: dgbmv_d
  integer :: n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing DGBMV (multi-size: n = 4)'
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
    integer :: msize, kl, ku
    real(8) :: alpha, alpha_d, alpha_orig, alpha_d_seed
    real(8) :: beta, beta_d, beta_orig, beta_d_seed
    real(8), dimension(:,:), allocatable :: a, a_d, a_orig, a_d_seed
    real(8), dimension(:), allocatable :: x, x_d, x_orig, x_d_seed
    real(8), dimension(:), allocatable :: y, y_d, y_orig, y_d_seed
    integer :: band_row, j
    real(4) :: temp_real
    ksize = max(0, n - 1)
    msize = n
    nsize = n
    kl = 1
    ku = 1
    lda_val = kl + ku + 1
    incx_val = 1
    incy_val = 1
    uplo = 'U'
    trans = 'N'
    diag = 'N'
    allocate(a(lda_val, n), a_d(lda_val, n), a_orig(lda_val, n), a_d_seed(lda_val, n))
    allocate(x(n), x_d(n), x_orig(n), x_d_seed(n))
    allocate(y(n), y_d(n), y_orig(n), y_d_seed(n))
    ! Initialize a as general band matrix (kl, ku band storage)
    do j = 1, n
    do band_row = max(1, ku+2-j), min(kl+ku+1, ku+msize-j+1)
    call random_number(temp_real)
    a(band_row, j) = temp_real * 2.0 - 1.0
    end do
    end do
    ! Keep direction consistent with general band (kl, ku): only band entries used
    do j = 1, n
    do band_row = max(1, ku+2-j), min(kl+ku+1, ku+msize-j+1)
    call random_number(temp_real)
    a_d(band_row, j) = temp_real * 2.0 - 1.0
    end do
    end do
    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(alpha_d)
    alpha_d = alpha_d * 2.0d0 - 1.0d0
    call random_number(beta)
    beta = beta * 2.0d0 - 1.0d0
    call random_number(beta_d)
    beta_d = beta_d * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    call random_number(x_d)
    x_d = x_d * 2.0d0 - 1.0d0
    call random_number(y)
    y = y * 2.0d0 - 1.0d0
    call random_number(y_d)
    y_d = y_d * 2.0d0 - 1.0d0
    write(*,*) 'Testing DGBMV (n =', n, ')'
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
    call dgbmv_d(trans, msize, nsize, kl, ku, alpha, alpha_d, a, a_d, lda_val, x, x_d, incx_val, beta, beta_d, y, y_d, incy_val)
    call check_derivatives_numerically_band_gbmv(n, lda_val, msize, nsize, kl, ku, trans, incx_val, incy_val, alpha_orig, alpha_d_seed, beta_orig, beta_d_seed, a_orig, a_d_seed, x_orig, x_d_seed, y_orig, y_d_seed, y_d, passed)
    deallocate(a, a_d, a_orig, a_d_seed, x, x_d, x_orig, x_d_seed)
    deallocate(y, y_d, y_orig, y_d_seed)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically_band_gbmv(n, lda_val, msize, nsize, kl, ku, trans, incx_val, incy_val, alpha_orig, alpha_d_seed, beta_orig, beta_d_seed, a_orig, a_d_seed, x_orig, x_d_seed, y_orig, y_d_seed, y_d_out, passed)
    implicit none
    integer, intent(in) :: n, lda_val, msize, nsize, kl, ku, incx_val, incy_val
    character, intent(in) :: trans
    real(8), intent(in) :: alpha_orig, alpha_d_seed, beta_orig, beta_d_seed
    real(8), intent(in) :: a_orig(lda_val, n), a_d_seed(lda_val, n), x_orig(n), x_d_seed(n), y_orig(n), y_d_seed(n), y_d_out(n)
    logical, intent(out) :: passed
    real(8), parameter :: h = 1.0e-7
    real(8) :: abs_error, abs_ref, err_bound
    real(8), dimension(n) :: y_fwd, y_bwd, y_t
    real(8) :: alpha_t, beta_t
    real(8), dimension(n) :: x_t
    real(8), dimension(lda_val, n) :: a_t
    integer :: ii
    logical :: has_err
    has_err = .false.
    alpha_t = alpha_orig + h * alpha_d_seed
    beta_t = beta_orig + h * beta_d_seed
    a_t = a_orig + h * a_d_seed
    x_t = x_orig + h * x_d_seed
    y_t = y_orig + h * y_d_seed
    call dgbmv(trans, msize, nsize, kl, ku, alpha_t, a_t, lda_val, x_t, incx_val, beta_t, y_t, incy_val)
    y_fwd = y_t
    alpha_t = alpha_orig - h * alpha_d_seed
    beta_t = beta_orig - h * beta_d_seed
    a_t = a_orig - h * a_d_seed
    x_t = x_orig - h * x_d_seed
    y_t = y_orig - h * y_d_seed
    call dgbmv(trans, msize, nsize, kl, ku, alpha_t, a_t, lda_val, x_t, incx_val, beta_t, y_t, incy_val)
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
  end subroutine check_derivatives_numerically_band_gbmv
end program test_dgbmv