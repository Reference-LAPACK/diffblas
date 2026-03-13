! Test program for SSBMV vector reverse - BLAS2 band
! Generated automatically by run_tapenade_blas.py
! Multi-size outlined run_test_for_size(n, passed, nbdirs)

program test_ssbmv_vector_reverse
  implicit none
  external :: ssbmv
  external :: ssbmv_bv
  integer :: nbdirs, n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing SSBMV (Vector Reverse band, multi-size: n = 4)'
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
    real(4) :: alpha, beta
    real(4), dimension(:), allocatable :: alphab, betab
    real(4), dimension(:,:), allocatable :: a
    real(4), dimension(:,:,:), allocatable :: ab
    real(4), dimension(:), allocatable :: x, y
    real(4), dimension(:,:), allocatable :: xb, yb, xb_seed, yb_seed
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
    allocate(a(lda_val, n), ab(nbdirs, lda_val, n), x(n), xb(nbdirs, n), y(n), yb(nbdirs, n), yb_seed(nbdirs, n), alphab(nbdirs), betab(nbdirs))
    ! Initialize a as symmetric band matrix (upper band storage)
    ! A(band_row, j) = full(i,j) with band_row = ksize+1+i-j, i = max(1,j-ksize)..j
    do j = 1, n
    do band_row = max(1, ksize+2-j), ksize+1
    call random_number(temp_real)
    a(band_row, j) = temp_real * 2.0 - 1.0  ! Scale to [-1,1]
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
    ab = 0.0d0
    alphab = 0.0d0
    betab = 0.0d0
    xb = 0.0d0
    ! Seed for vector reverse: output adjoint yb is the seed per direction
    call random_number(yb)
    yb = yb * 2.0d0 - 1.0d0
    yb_seed = yb
    write(*,*) 'Testing SSBMV (Vector Reverse band, n =', n, ')'
    call set_ISIZE1OFX(n)
    call set_ISIZE2OFA(n)
    call ssbmv_bv(uplo, nsize, ksize, alpha, alphab, a, ab, lda_val, x, xb, incx_val, beta, betab, y, yb, incy_val, nbdirs)
    call set_ISIZE1OFX(-1)
    call set_ISIZE2OFA(-1)
    write(*,*) 'Function calls completed successfully'
    call check_vjp_numerically_band_vec(n, nbdirs, lda_val, ksize, uplo, nsize, incx_val, incy_val, alpha, alphab, beta, betab, a, ab, x, xb, y, yb_seed, yb, passed)
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

  subroutine check_vjp_numerically_band_vec(n, nbdirs, lda_val, ksize, uplo, nsize, incx_val, incy_val, alpha, alphab, beta, betab, a, ab, x, xb, y, yb_seed, yb, passed)
    implicit none
    integer, intent(in) :: n, nbdirs, lda_val, ksize, nsize, incx_val, incy_val
    character, intent(in) :: uplo
    real(4), intent(in) :: alpha, beta
    real(4), intent(in) :: alphab(nbdirs), betab(nbdirs)
    real(4), intent(in) :: a(lda_val, n), ab(nbdirs, lda_val, n), x(n), xb(nbdirs, n), y(n), yb_seed(nbdirs, n), yb(nbdirs, n)
    logical, intent(out) :: passed
    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_fd, vjp_ad, abs_error, abs_ref, err_bound, relative_error, max_re
    real(4), dimension(n) :: y_plus, y_minus, y_t, y_central_diff
    real(4) :: alpha_t, beta_t, alpha_dir, beta_dir
    real(4), dimension(n) :: x_t, x_dir, y_dir
    real(4), dimension(lda_val, n) :: a_t, a_dir
    real(4), dimension(:), allocatable :: temp_products
    integer :: i, j, band_row, n_products, k
    logical :: has_err
    has_err = .false.
    max_re = 0.0d0
    allocate(temp_products(n + n + n + (ksize+1)*n + 2))
    do k = 1, nbdirs
    ! Random direction for this k
    call random_number(alpha_dir)
    alpha_dir = alpha_dir * 2.0d0 - 1.0d0
    call random_number(beta_dir)
    beta_dir = beta_dir * 2.0d0 - 1.0d0
    do j = 1, n
      do band_row = max(1, ksize+2-j), ksize+1
        call random_number(a_dir(band_row, j))
        a_dir(band_row, j) = a_dir(band_row, j) * 2.0d0 - 1.0d0
      end do
    end do
    call random_number(x_dir)
    x_dir = x_dir * 2.0d0 - 1.0d0
    call random_number(y_dir)
    y_dir = y_dir * 2.0d0 - 1.0d0
    ! Forward perturbation: f(inputs + h*direction)
    alpha_t = alpha + h * alpha_dir
    beta_t = beta + h * beta_dir
    a_t = a
    do j = 1, n
      do band_row = max(1, ksize+2-j), ksize+1
        a_t(band_row, j) = a(band_row, j) + h * a_dir(band_row, j)
      end do
    end do
    x_t = x + h * x_dir
    y_t = y + h * y_dir
    call ssbmv(uplo, nsize, ksize, alpha_t, a_t, lda_val, x_t, incx_val, beta_t, y_t, incy_val)
    y_plus = y_t
    ! Backward perturbation: f(inputs - h*direction)
    alpha_t = alpha - h * alpha_dir
    beta_t = beta - h * beta_dir
    a_t = a
    do j = 1, n
      do band_row = max(1, ksize+2-j), ksize+1
        a_t(band_row, j) = a(band_row, j) - h * a_dir(band_row, j)
      end do
    end do
    x_t = x - h * x_dir
    y_t = y - h * y_dir
    call ssbmv(uplo, nsize, ksize, alpha_t, a_t, lda_val, x_t, incx_val, beta_t, y_t, incy_val)
    y_minus = y_t
    y_central_diff = (y_plus - y_minus) / (2.0d0 * h)
    vjp_fd = 0.0d0
    n_products = n
    do i = 1, n
      temp_products(i) = yb_seed(k,i) * y_central_diff(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do
    ! VJP(AD) = direction^T @ adjoint
    vjp_ad = 0.0d0
    vjp_ad = vjp_ad + alpha_dir * alphab(k)
    vjp_ad = vjp_ad + beta_dir * betab(k)
    n_products = 0
    do j = 1, n
      do band_row = max(1, ksize+2-j), ksize+1
        n_products = n_products + 1
        temp_products(n_products) = a_dir(band_row,j) * ab(k,band_row,j)
      end do
    end do
    do i = 1, n
      n_products = n_products + 1
      temp_products(n_products) = x_dir(i) * xb(k,i)
    end do
    do i = 1, n
      n_products = n_products + 1
      temp_products(n_products) = y_dir(i) * yb(k,i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    abs_error = abs(vjp_fd - vjp_ad)
    abs_ref = abs(vjp_ad)
    err_bound = 1.0e-3 + 1.0e-3 * abs_ref
    if (abs_error > err_bound) has_err = .true.
    relative_error = abs_error / max(abs_ref, 1.0d-10)
    if (relative_error > max_re) max_re = relative_error
    end do
    deallocate(temp_products)
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Maximum relative error:', max_re
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    passed = .not. has_err
    if (.not. passed) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
  end subroutine check_vjp_numerically_band_vec
  subroutine sort_array(arr, n)
    implicit none
    integer, intent(in) :: n
    real(4), dimension(n), intent(inout) :: arr
    integer :: i, j, min_idx
    real(4) :: temp
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
end program test_ssbmv_vector_reverse