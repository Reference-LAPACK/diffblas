! Test program for ZTBMV vector forward - BLAS2 band
! Generated automatically by run_tapenade_blas.py
! Multi-size outlined run_test_for_size(n, passed, nbdirs) - band

program test_ztbmv_vector_forward
  implicit none
  external :: ztbmv
  external :: ztbmv_dv
  integer :: nbdirs, n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing ZTBMV (Vector Forward band, multi-size: n = 4)'
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
    complex(8) :: alpha, beta
    complex(8), dimension(:,:), allocatable :: a, a_orig
    complex(8), dimension(:,:,:), allocatable :: a_dv, a_dv_seed
    complex(8), dimension(:), allocatable :: x, y, x_orig, y_orig
    complex(8), dimension(:,:), allocatable :: x_dv, y_dv, x_dv_seed, y_dv_seed
    integer :: band_row, j, idir
    real(4) :: temp_real, temp_imag
    ksize = max(0, n - 1)
    nsize = n
    lda_val = ksize + 1
    incx_val = 1
    incy_val = 1
    uplo = 'U'
    trans = 'N'
    diag = 'N'
    allocate(a(lda_val, n), a_orig(lda_val, n), a_dv(nbdirs, lda_val, n), a_dv_seed(nbdirs, lda_val, n), x(n), x_orig(n), x_dv(nbdirs, n), x_dv_seed(nbdirs, n))
    ! Initialize a as triangular band matrix (upper band storage)
    do j = 1, n
    do band_row = max(1, ksize+2-j), ksize+1
    call random_number(temp_real)
    call random_number(temp_imag)
    a(band_row, j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
    end do
    do idir = 1, nbdirs
      do j = 1, n
        do band_row = max(1, ksize+2-j), ksize+1
          call random_number(temp_real)
          call random_number(temp_imag)
          a_dv(idir, band_row, j) = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(a_dv))
        end do
      end do
    end do
    call random_number(temp_real)
    call random_number(temp_imag)
    alpha = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(alpha))
    do j = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      x(j) = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(x))
    end do
    do idir = 1, nbdirs
      do j = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        x_dv(idir,j) = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(x_dv))
      end do
    end do
    write(*,*) 'Testing ZTBMV (Vector Forward band, n =', n, ')'
    a_orig = a
    x_orig = x
    a_dv_seed = a_dv
    x_dv_seed = x_dv
    call ztbmv_dv(uplo, trans, diag, nsize, ksize, a, a_dv, lda_val, x, x_dv, incx_val, nbdirs)
    call check_derivatives_numerically_band_tri(n, nbdirs, lda_val, ksize, uplo, trans, diag, nsize, incx_val, a_orig, a_dv_seed, x_orig, x_dv_seed, x_dv, passed)
    deallocate(a, a_orig, a_dv, a_dv_seed, x, x_orig, x_dv, x_dv_seed)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically_band_tri(n, nbdirs, lda_val, ksize, uplo, trans, diag, nsize, incx_val, a_orig, a_dv_seed, x_orig, x_dv_seed, x_dv_out, passed)
    implicit none
    integer, intent(in) :: n, nbdirs, lda_val, ksize, nsize, incx_val
    character, intent(in) :: uplo, trans, diag
    complex(8), intent(in) :: a_orig(lda_val, n), a_dv_seed(nbdirs, lda_val, n), x_orig(n), x_dv_seed(nbdirs, n), x_dv_out(nbdirs, n)
    logical, intent(out) :: passed
    real(8), parameter :: h = 1.0e-7
    real(8) :: abs_error, abs_ref, err_bound
    complex(8) :: central_diff, ad_result
    logical :: has_err
    complex(8), dimension(n) :: x_fwd, x_bwd, x_t
    complex(8), dimension(lda_val, n) :: a_t
    integer :: i, idir
    has_err = .false.
    do idir = 1, nbdirs
      a_t = a_orig + h * a_dv_seed(idir,:,:)
      x_t = x_orig + h * x_dv_seed(idir,:)
      call ztbmv(uplo, trans, diag, nsize, ksize, a_t, lda_val, x_t, incx_val)
      x_fwd = x_t
      a_t = a_orig - h * a_dv_seed(idir,:,:)
      x_t = x_orig - h * x_dv_seed(idir,:)
      call ztbmv(uplo, trans, diag, nsize, ksize, a_t, lda_val, x_t, incx_val)
      x_bwd = x_t
      do i = 1, min(3, n)
        central_diff = (x_fwd(i) - x_bwd(i)) / (2.0e0 * h)
        ad_result = x_dv_out(idir, i)
        abs_error = abs(central_diff - ad_result)
        abs_ref = abs(ad_result)
        err_bound = 1.0e-5 + 1.0e-5 * abs_ref
        if (abs_error > err_bound) has_err = .true.
      end do
    end do
    passed = .not. has_err
    if (has_err) write(*,*) 'FAIL: Band vector forward derivatives'
    if (.not. has_err) write(*,*) 'PASS: Band vector forward derivatives'
  end subroutine check_derivatives_numerically_band_tri
end program test_ztbmv_vector_forward