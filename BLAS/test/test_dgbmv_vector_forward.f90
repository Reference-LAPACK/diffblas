! Test program for DGBMV vector forward - BLAS2 band
! Generated automatically by run_tapenade_blas.py
! Multi-size outlined run_test_for_size(n, passed, nbdirs) - band

program test_dgbmv_vector_forward
  implicit none
  external :: dgbmv
  external :: dgbmv_dv
  integer :: nbdirs, n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing DGBMV (Vector Forward band, multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    nbdirs = test_sizes(i)
    call run_test_for_size(n_test, passed, nbdirs)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) write(*,*) 'PASS: Vector forward band - all sizes OK'
  if (.not. all_passed) write(*,*) 'FAIL: Vector forward band - errors'
contains
  subroutine run_test_for_size(n, passed, nbdirs)
    implicit none
    integer, intent(in) :: n, nbdirs
    logical, intent(out) :: passed
    character :: uplo, trans, diag
    integer :: nsize, ksize, lda_val, incx_val, incy_val
    integer :: msize, kl, ku
    real(8) :: alpha, beta
    real(8), dimension(:,:), allocatable :: a, a_orig
    real(8), dimension(:,:,:), allocatable :: a_dv, a_dv_seed
    real(8), dimension(:), allocatable :: x, y, x_orig, y_orig
    real(8), dimension(:,:), allocatable :: x_dv, y_dv, x_dv_seed, y_dv_seed
    real(8), dimension(nbdirs) :: alpha_dv, beta_dv, alpha_dv_seed, beta_dv_seed
    integer :: band_row, j, idir
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
    allocate(a(lda_val, n), a_orig(lda_val, n), a_dv(nbdirs, lda_val, n), a_dv_seed(nbdirs, lda_val, n), x(n), x_orig(n), x_dv(nbdirs, n), x_dv_seed(nbdirs, n), y(n), y_orig(n), y_dv(nbdirs, n), y_dv_seed(nbdirs, n))
    ! Initialize a as general band matrix (kl, ku band storage)
    do j = 1, n
    do band_row = max(1, ku+2-j), min(kl+ku+1, ku+msize-j+1)
    call random_number(temp_real)
    a(band_row, j) = temp_real * 2.0 - 1.0
    end do
    end do
    do idir = 1, nbdirs
      do j = 1, n
        do band_row = max(1, ku+2-j), min(kl+ku+1, ku+msize-j+1)
          call random_number(temp_real)
          a_dv(idir, band_row, j) = temp_real * 2.0d0 - 1.0d0
        end do
      end do
    end do
    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(beta)
    beta = beta * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    call random_number(y)
    y = y * 2.0d0 - 1.0d0
    call random_number(x_dv)
    x_dv = x_dv * 2.0d0 - 1.0d0
    call random_number(y_dv)
    y_dv = y_dv * 2.0d0 - 1.0d0
    do idir = 1, nbdirs
      call random_number(alpha_dv(idir))
      alpha_dv(idir) = alpha_dv(idir) * 2.0d0 - 1.0d0
      call random_number(beta_dv(idir))
      beta_dv(idir) = beta_dv(idir) * 2.0d0 - 1.0d0
    end do
    write(*,*) 'Testing DGBMV (Vector Forward band, n =', n, ')'
    a_orig = a
    x_orig = x
    a_dv_seed = a_dv
    x_dv_seed = x_dv
    y_orig = y
    y_dv_seed = y_dv
    alpha_dv_seed = alpha_dv
    beta_dv_seed = beta_dv
    call dgbmv_dv(trans, msize, nsize, kl, ku, alpha, alpha_dv, a, a_dv, lda_val, x, x_dv, incx_val, beta, beta_dv, y, y_dv, incy_val, nbdirs)
    call check_derivatives_numerically_band_gbmv(n, nbdirs, lda_val, msize, nsize, kl, ku, trans, incx_val, incy_val, alpha, alpha_dv_seed, beta, beta_dv_seed, a_orig, a_dv_seed, x_orig, x_dv_seed, y_orig, y_dv_seed, y_dv, passed)
    deallocate(a, a_orig, a_dv, a_dv_seed, x, x_orig, x_dv, x_dv_seed, y, y_orig, y_dv, y_dv_seed)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically_band_gbmv(n, nbdirs, lda_val, msize, nsize, kl, ku, trans, incx_val, incy_val, alpha, alpha_dv_seed, beta, beta_dv_seed, a_orig, a_dv_seed_mat, x_orig, x_dv_seed, y_orig, y_dv_seed, y_dv_out, passed)
    implicit none
    integer, intent(in) :: n, nbdirs, lda_val, msize, nsize, kl, ku, incx_val, incy_val
    character, intent(in) :: trans
    real(8), intent(in) :: alpha, beta
    real(8), intent(in) :: alpha_dv_seed(nbdirs), beta_dv_seed(nbdirs)
    real(8), intent(in) :: a_orig(lda_val, n), a_dv_seed_mat(nbdirs, lda_val, n), x_orig(n), x_dv_seed(nbdirs, n), y_orig(n), y_dv_seed(nbdirs, n), y_dv_out(nbdirs, n)
    logical, intent(out) :: passed
    real(8), parameter :: h = 1.0e-7
    real(8) :: abs_error, abs_ref, err_bound
    real(8) :: central_diff, ad_result
    logical :: has_err
    real(8), dimension(n) :: y_fwd, y_bwd, y_t
    real(8) :: alpha_t, beta_t
    real(8), dimension(n) :: x_t
    real(8), dimension(lda_val, n) :: a_t
    integer :: i, idir
    has_err = .false.
    do idir = 1, nbdirs
      alpha_t = alpha + h * alpha_dv_seed(idir)
      beta_t = beta + h * beta_dv_seed(idir)
      a_t = a_orig + h * a_dv_seed_mat(idir,:,:)
      x_t = x_orig + h * x_dv_seed(idir,:)
      y_t = y_orig + h * y_dv_seed(idir,:)
      call dgbmv(trans, msize, nsize, kl, ku, alpha_t, a_t, lda_val, x_t, incx_val, beta_t, y_t, incy_val)
      y_fwd = y_t
      alpha_t = alpha - h * alpha_dv_seed(idir)
      beta_t = beta - h * beta_dv_seed(idir)
      a_t = a_orig - h * a_dv_seed_mat(idir,:,:)
      x_t = x_orig - h * x_dv_seed(idir,:)
      y_t = y_orig - h * y_dv_seed(idir,:)
      call dgbmv(trans, msize, nsize, kl, ku, alpha_t, a_t, lda_val, x_t, incx_val, beta_t, y_t, incy_val)
      y_bwd = y_t
      do i = 1, min(3, n)
        central_diff = (y_fwd(i) - y_bwd(i)) / (2.0e0 * h)
        ad_result = y_dv_out(idir, i)
        abs_error = abs(central_diff - ad_result)
        abs_ref = abs(ad_result)
        err_bound = 1.0e-5 + 1.0e-5 * abs_ref
        if (abs_error > err_bound) has_err = .true.
      end do
    end do
    passed = .not. has_err
    if (has_err) write(*,*) 'FAIL: Band vector forward derivatives'
    if (.not. has_err) write(*,*) 'PASS: Band vector forward derivatives'
  end subroutine check_derivatives_numerically_band_gbmv
end program test_dgbmv_vector_forward