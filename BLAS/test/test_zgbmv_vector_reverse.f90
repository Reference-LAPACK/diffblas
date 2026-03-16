! Test program for ZGBMV vector reverse - BLAS2 band
! Generated automatically by run_tapenade_blas.py
! Multi-size outlined run_test_for_size(n, passed, nbdirs)

program test_zgbmv_vector_reverse
  implicit none
  external :: zgbmv
  external :: zgbmv_bv
  integer :: nbdirs, n_test, seed_array(33), test_sizes(3), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing ZGBMV (Vector Reverse band, multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 3
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
    integer :: msize, kl, ku
    complex(8) :: alpha, beta
    complex(8), dimension(:), allocatable :: alphab, betab
    complex(8), dimension(:,:), allocatable :: a
    complex(8), dimension(:,:,:), allocatable :: ab
    complex(8), dimension(:), allocatable :: x, y
    complex(8), dimension(:,:), allocatable :: xb, yb, xb_seed, yb_seed
    integer :: band_row, j
    real(4) :: temp_real, temp_imag
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
    allocate(a(lda_val, n), ab(nbdirs, lda_val, n), x(n), xb(nbdirs, n), y(n), yb(nbdirs, n), yb_seed(nbdirs, n), alphab(nbdirs), betab(nbdirs))
    ! Initialize a as general band matrix (kl, ku band storage)
    do j = 1, n
    do band_row = max(1, ku+2-j), min(kl+ku+1, ku+msize-j+1)
    call random_number(temp_real)
    call random_number(temp_imag)
    a(band_row, j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
    end do
    call random_number(temp_real)
    call random_number(temp_imag)
    alpha = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(alpha))
    call random_number(temp_real)
    call random_number(temp_imag)
    beta = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(beta))
    do j = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      x(j) = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(x))
      call random_number(temp_real)
      call random_number(temp_imag)
      y(j) = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(y))
    end do
    ab = 0.0d0
    alphab = 0.0d0
    betab = 0.0d0
    xb = 0.0d0
    ! Seed for vector reverse: output adjoint yb is the seed per direction
    do j = 1, n
      do band_row = 1, nbdirs
        call random_number(temp_real)
        call random_number(temp_imag)
        yb(band_row, j) = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(yb))
      end do
    end do
    yb_seed = yb
    write(*,*) 'Testing ZGBMV (Vector Reverse band, n =', n, ')'
    call set_ISIZE1OFX(n)
    call set_ISIZE2OFA(n)
    call zgbmv_bv(trans, msize, nsize, kl, ku, alpha, alphab, a, ab, lda_val, x, xb, incx_val, beta, betab, y, yb, incy_val, nbdirs)
    call set_ISIZE1OFX(-1)
    call set_ISIZE2OFA(-1)
    write(*,*) 'Function calls completed successfully'
    call check_vjp_numerically_band_gbmv_vec(n, nbdirs, lda_val, msize, nsize, kl, ku, trans, incx_val, incy_val, alpha, alphab, beta, betab, a, ab, x, xb, y, yb_seed, yb, passed)
    if (allocated(a)) deallocate(a)
    if (allocated(ab)) deallocate(ab)
    if (allocated(x)) deallocate(x)
    if (allocated(xb)) deallocate(xb)
    if (allocated(y)) deallocate(y)
    if (allocated(yb)) deallocate(yb)
    if (allocated(yb_seed)) deallocate(yb_seed)
    if (allocated(alphab)) deallocate(alphab)
    if (allocated(betab)) deallocate(betab)
  end subroutine run_test_for_size

  subroutine check_vjp_numerically_band_gbmv_vec(n, nbdirs, lda_val, msize, nsize, kl, ku, trans, incx_val, incy_val, alpha, alphab, beta, betab, a, ab, x, xb, y, yb_seed, yb, passed)
    implicit none
    integer, intent(in) :: n, nbdirs, lda_val, msize, nsize, kl, ku, incx_val, incy_val
    character, intent(in) :: trans
    complex(8), intent(in) :: alpha, beta
    complex(8), intent(in) :: alphab(nbdirs), betab(nbdirs)
    complex(8), intent(in) :: a(lda_val, n), ab(nbdirs, lda_val, n), x(n), xb(nbdirs, n), y(n), yb_seed(nbdirs, n), yb(nbdirs, n)
    logical, intent(out) :: passed
    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_fd, vjp_ad, abs_error, abs_ref, err_bound, relative_error, max_re
    complex(8), dimension(n) :: y_plus, y_minus, y_t, y_central_diff
    complex(8) :: alpha_t, beta_t, alpha_dir, beta_dir
    complex(8), dimension(n) :: x_t, x_dir, y_dir
    complex(8), dimension(lda_val, n) :: a_t, a_dir
    real(8), dimension(:), allocatable :: temp_products
    real(kind(0.0d0)) :: tr, ti
    integer :: i, j, band_row, n_products, k
    logical :: has_err
    has_err = .false.
    max_re = 0.0d0
    allocate(temp_products(n + (kl+ku+1)*n + 2))
    do k = 1, nbdirs
    ! Random direction for this k (match BLAS1 scalar: direction^T @ adjoint)
    call random_number(tr)
    call random_number(ti)
    alpha_dir = cmplx(tr*2.0d0-1.0d0, ti*2.0d0-1.0d0, kind=kind(alpha_dir))
    call random_number(tr)
    call random_number(ti)
    beta_dir = cmplx(tr*2.0d0-1.0d0, ti*2.0d0-1.0d0, kind=kind(beta_dir))
    do j = 1, n
      do band_row = max(1, ku+2-j), min(kl+ku+1, ku+msize-j+1)
        call random_number(tr)
        call random_number(ti)
        a_dir(band_row, j) = cmplx(tr*2.0d0-1.0d0, ti*2.0d0-1.0d0, kind=kind(a_dir))
      end do
    end do
    do i = 1, n
      call random_number(tr)
      call random_number(ti)
      x_dir(i) = cmplx(tr*2.0d0-1.0d0, ti*2.0d0-1.0d0, kind=kind(x_dir))
      call random_number(tr)
      call random_number(ti)
      y_dir(i) = cmplx(tr*2.0d0-1.0d0, ti*2.0d0-1.0d0, kind=kind(y_dir))
    end do
    alpha_t = alpha + h * alpha_dir
    beta_t = beta + h * beta_dir
    a_t = a
    do j = 1, n
      do band_row = max(1, ku+2-j), min(kl+ku+1, ku+msize-j+1)
        a_t(band_row, j) = a(band_row, j) + h * a_dir(band_row, j)
      end do
    end do
    x_t = x + h * x_dir
    y_t = y + h * y_dir
    call zgbmv(trans, msize, nsize, kl, ku, alpha_t, a_t, lda_val, x_t, incx_val, beta_t, y_t, incy_val)
    y_plus = y_t
    alpha_t = alpha - h * alpha_dir
    beta_t = beta - h * beta_dir
    a_t = a
    do j = 1, n
      do band_row = max(1, ku+2-j), min(kl+ku+1, ku+msize-j+1)
        a_t(band_row, j) = a(band_row, j) - h * a_dir(band_row, j)
      end do
    end do
    x_t = x - h * x_dir
    y_t = y - h * y_dir
    call zgbmv(trans, msize, nsize, kl, ku, alpha_t, a_t, lda_val, x_t, incx_val, beta_t, y_t, incy_val)
    y_minus = y_t
    y_central_diff = (y_plus - y_minus) / (2.0d0 * h)
    vjp_fd = 0.0d0
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(yb_seed(k,i)) * y_central_diff(i))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do
    vjp_ad = 0.0d0
    vjp_ad = vjp_ad + real(conjg(alpha_dir) * alphab(k))
    vjp_ad = vjp_ad + real(conjg(beta_dir) * betab(k))
    n_products = 0
    do j = 1, n
      do band_row = max(1, ku+2-j), min(kl+ku+1, ku+msize-j+1)
        n_products = n_products + 1
        temp_products(n_products) = real(conjg(a_dir(band_row,j)) * ab(k,band_row,j))
      end do
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    do i = 1, n
      temp_products(i) = real(conjg(x_dir(i)) * xb(k,i))
    end do
    call sort_array(temp_products, n)
    do i = 1, n
      vjp_ad = vjp_ad + temp_products(i)
    end do
    do i = 1, n
      temp_products(i) = real(conjg(y_dir(i)) * yb(k,i))
    end do
    call sort_array(temp_products, n)
    do i = 1, n
      vjp_ad = vjp_ad + temp_products(i)
    end do
    abs_error = abs(vjp_fd - vjp_ad)
    abs_ref = abs(vjp_ad)
    err_bound = 1.0e-5 + 1.0e-5 * abs_ref
    if (abs_error > err_bound) has_err = .true.
    relative_error = abs_error / max(abs_ref, 1.0d-10)
    if (relative_error > max_re) max_re = relative_error
    end do
    deallocate(temp_products)
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Maximum relative error:', max_re
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_err
    if (.not. passed) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
  end subroutine check_vjp_numerically_band_gbmv_vec
  subroutine sort_array(arr, n)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(inout) :: arr
    integer :: i, j, min_idx
    real(8) :: temp
    do i = 1, n-1
      min_idx = i
      do j = i+1, n
        if (abs(arr(j)) < abs(arr(min_idx))) min_idx = j
      end do
      if (min_idx /= i) then
        temp = arr(i)
        arr(i) = arr(min_idx)
        arr(min_idx) = temp
      end if
    end do
  end subroutine sort_array
end program test_zgbmv_vector_reverse