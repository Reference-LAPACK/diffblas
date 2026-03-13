! Test program for CTBMV vector forward - BLAS2 band
! Generated automatically by run_tapenade_blas.py
! Multi-size outlined run_test_for_size(n, passed, nbdirs) - band

program test_ctbmv_vector_forward
  implicit none
  external :: ctbmv
  external :: ctbmv_dv
  integer :: nbdirs, n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing CTBMV (Vector Forward band, multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    nbdirs = test_sizes(i)
    call run_test_for_size(n_test, passed, nbdirs)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) write(*,*) 'PASS: All sizes completed successfully'
  if (.not. all_passed) write(*,*) 'FAIL: One or more sizes had derivative errors'
contains
  subroutine run_test_for_size(n, passed, nbdirs)
    implicit none
    integer, intent(in) :: n, nbdirs
    logical, intent(out) :: passed
    character :: uplo, trans, diag
    integer :: nsize, ksize, lda_val, incx_val, incy_val
    complex(4) :: alpha, beta
    complex(4), dimension(:,:), allocatable :: a, a_orig
    complex(4), dimension(:,:,:), allocatable :: a_dv, a_dv_seed
    complex(4), dimension(:), allocatable :: x, y, x_orig, y_orig
    complex(4), dimension(:,:), allocatable :: x_dv, y_dv, x_dv_seed, y_dv_seed
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
    write(*,*) 'Testing CTBMV (Vector Forward band, n =', n, ')'
    a_orig = a
    x_orig = x
    a_dv_seed = a_dv
    x_dv_seed = x_dv
    call ctbmv_dv(uplo, trans, diag, nsize, ksize, a, a_dv, lda_val, x, x_dv, incx_val, nbdirs)
    write(*,*) 'Function calls completed successfully'
    call check_derivatives_numerically_band_tri(n, nbdirs, lda_val, ksize, uplo, trans, diag, nsize, incx_val, a_orig, a_dv_seed, x_orig, x_dv_seed, x_dv, passed)
    deallocate(a, a_orig, a_dv, a_dv_seed, x, x_orig, x_dv, x_dv_seed)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically_band_tri(n, nbdirs, lda_val, ksize, uplo, trans, diag, nsize, incx_val, a_orig, a_dv_seed, x_orig, x_dv_seed, x_dv_out, passed)
    implicit none
    integer, intent(in) :: n, nbdirs, lda_val, ksize, nsize, incx_val
    character, intent(in) :: uplo, trans, diag
    complex(4), intent(in) :: a_orig(lda_val, n), a_dv_seed(nbdirs, lda_val, n), x_orig(n), x_dv_seed(nbdirs, n), x_dv_out(nbdirs, n)
    logical, intent(out) :: passed
    real(4), parameter :: h = 1.0e-3
    real(4) :: abs_error, abs_ref, err_bound, max_error, relative_error
    complex(4) :: central_diff, ad_result
    logical :: has_err
    complex(4), dimension(n) :: x_fwd, x_bwd, x_t
    complex(4), dimension(lda_val, n) :: a_t
    integer :: i, idir, j, band_row
    has_err = .false.
    max_error = 0.0e0
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    do idir = 1, nbdirs
      a_t = a_orig
      do j = 1, n
        do band_row = max(1, ksize+2-j), ksize+1
          a_t(band_row, j) = a_orig(band_row, j) + h * a_dv_seed(idir, band_row, j)
        end do
      end do
      x_t = x_orig + h * x_dv_seed(idir,:)
      call ctbmv(uplo, trans, diag, nsize, ksize, a_t, lda_val, x_t, incx_val)
      x_fwd = x_t
      a_t = a_orig
      do j = 1, n
        do band_row = max(1, ksize+2-j), ksize+1
          a_t(band_row, j) = a_orig(band_row, j) - h * a_dv_seed(idir, band_row, j)
        end do
      end do
      x_t = x_orig - h * x_dv_seed(idir,:)
      call ctbmv(uplo, trans, diag, nsize, ksize, a_t, lda_val, x_t, incx_val)
      x_bwd = x_t
      do i = 1, min(3, n)
        central_diff = (x_fwd(i) - x_bwd(i)) / (2.0e0 * h)
        ad_result = x_dv_out(idir, i)
        abs_error = abs(central_diff - ad_result)
        abs_ref = abs(ad_result)
        err_bound = 1.0e-3 + 1.0e-3 * abs_ref
        if (abs_error > err_bound) has_err = .true.
        relative_error = abs_error / max(abs_ref, 1.0e-10)
        if (relative_error > max_error) max_error = relative_error
      end do
    end do
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    passed = .not. has_err
    if (has_err) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
  end subroutine check_derivatives_numerically_band_tri
end program test_ctbmv_vector_forward