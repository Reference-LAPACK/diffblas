! Test program for DGBMV reverse mode (adjoint) - BLAS2 band
! Generated automatically by run_tapenade_blas.py
! Multi-size outlined run_test_for_size(n) - band (declarations in subroutines)

program test_dgbmv_reverse
  implicit none
  external :: dgbmv
  external :: dgbmv_b
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
    real(8) :: alpha, alphab
    real(8) :: beta, betab
    real(8), dimension(:,:), allocatable :: a, ab
    real(8), dimension(:), allocatable :: x, xb
    real(8), dimension(:), allocatable :: y, yb, yb_seed
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
    allocate(a(lda_val, n), ab(lda_val, n), x(n), xb(n))
    allocate(y(n), yb(n), yb_seed(n))
    ! Initialize a as general band matrix (kl, ku band storage)
    do j = 1, n
    do band_row = max(1, ku+2-j), min(kl+ku+1, ku+msize-j+1)
    call random_number(temp_real)
    a(band_row, j) = temp_real * 2.0 - 1.0
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
    alphab = 0.0d0
    ab = 0.0d0
    xb = 0.0d0
    ! Seed for reverse mode: output adjoint yb is the seed (d(scalar)/d(y))
    call random_number(yb)
    yb = yb * 2.0d0 - 1.0d0
    yb_seed = yb
    write(*,*) 'Testing DGBMV (n =', n, ')'
    call set_ISIZE1OFX(n)
    call set_ISIZE2OFA(lda_val)
    call dgbmv_b(trans, msize, nsize, kl, ku, alpha, alphab, a, ab, lda_val, x, xb, incx_val, beta, betab, y, yb, incy_val)
    call set_ISIZE1OFX(-1)
    call set_ISIZE2OFA(-1)
    write(*,*) 'Function calls completed successfully'
    call check_vjp_numerically_band_gbmv(n, lda_val, msize, nsize, kl, ku, trans, incx_val, incy_val, alpha, alphab, beta, betab, a, ab, x, xb, y, yb_seed, yb, passed)
    deallocate(a, ab, x, xb)
    deallocate(y, yb, yb_seed)
  end subroutine run_test_for_size

  subroutine check_vjp_numerically_band_gbmv(n, lda_val, msize, nsize, kl, ku, trans, incx_val, incy_val, alpha, alphab, beta, betab, a, ab, x, xb, y, yb_seed, yb, passed)
    implicit none
    integer, intent(in) :: n, lda_val, msize, nsize, kl, ku, incx_val, incy_val
    character, intent(in) :: trans
    real(8), intent(in) :: alpha, alphab, beta, betab
    real(8), intent(in) :: a(lda_val, n), ab(lda_val, n), x(n), xb(n), y(n), yb_seed(n), yb(n)
    logical, intent(out) :: passed
    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_fd, vjp_ad, abs_error, abs_ref, err_bound, relative_error
    real(8), dimension(n) :: y_plus, y_minus, y_t, y_central_diff
    real(8) :: alpha_t, beta_t, alpha_dir, beta_dir
    real(8), dimension(n) :: x_t, x_dir, y_dir
    real(8), dimension(lda_val, n) :: a_t, a_dir
    real(8), dimension(:), allocatable :: temp_products
    integer :: i, j, band_row, n_products
    allocate(temp_products(n + (kl+ku+1)*n + 2))
    ! Random direction for FD (match BLAS1 reference: direction^T @ adjoint)
    call random_number(alpha_dir)
    alpha_dir = alpha_dir * 2.0d0 - 1.0d0
    call random_number(beta_dir)
    beta_dir = beta_dir * 2.0d0 - 1.0d0
    do j = 1, n
      do band_row = max(1, ku+2-j), min(kl+ku+1, ku+msize-j+1)
        call random_number(a_dir(band_row, j))
        a_dir(band_row, j) = a_dir(band_row, j) * 2.0d0 - 1.0d0
      end do
    end do
    call random_number(x_dir)
    x_dir = x_dir * 2.0d0 - 1.0d0
    call random_number(y_dir)
    y_dir = y_dir * 2.0d0 - 1.0d0
    ! Forward perturbation: f(x + h*direction)
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
    call dgbmv(trans, msize, nsize, kl, ku, alpha_t, a_t, lda_val, x_t, incx_val, beta_t, y_t, incy_val)
    y_plus = y_t
    ! Backward perturbation: f(x - h*direction)
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
    call dgbmv(trans, msize, nsize, kl, ku, alpha_t, a_t, lda_val, x_t, incx_val, beta_t, y_t, incy_val)
    y_minus = y_t
    y_central_diff = (y_plus - y_minus) / (2.0d0 * h)
    vjp_fd = 0.0d0
    n_products = n
    do i = 1, n
      temp_products(i) = yb_seed(i) * y_central_diff(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do
    ! VJP(AD) = direction^T @ adjoint (BLAS1 reference)
    vjp_ad = 0.0d0
    vjp_ad = vjp_ad + alpha_dir * alphab
    vjp_ad = vjp_ad + beta_dir * betab
    n_products = 0
    do j = 1, n
      do band_row = max(1, ku+2-j), min(kl+ku+1, ku+msize-j+1)
        n_products = n_products + 1
        temp_products(n_products) = a_dir(band_row,j) * ab(band_row,j)
      end do
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    do i = 1, n
      temp_products(i) = x_dir(i) * xb(i)
    end do
    call sort_array(temp_products, n)
    do i = 1, n
      vjp_ad = vjp_ad + temp_products(i)
    end do
    do i = 1, n
      temp_products(i) = y_dir(i) * yb(i)
    end do
    call sort_array(temp_products, n)
    do i = 1, n
      vjp_ad = vjp_ad + temp_products(i)
    end do
    abs_error = abs(vjp_fd - vjp_ad)
    abs_ref = abs(vjp_ad)
    err_bound = 1.0e-5 + 1.0e-5 * abs_ref
    relative_error = 0.0d0
    if (abs_ref > 1.0d-10) relative_error = abs_error / abs_ref
    deallocate(temp_products)
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Maximum relative error:', relative_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = abs_error <= err_bound
    if (.not. passed) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
  end subroutine check_vjp_numerically_band_gbmv
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
end program test_dgbmv_reverse