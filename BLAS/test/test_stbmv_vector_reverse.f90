! Test program for STBMV vector reverse - BLAS2 band
! Generated automatically by run_tapenade_blas.py
! Multi-size outlined run_test_for_size(n, passed, nbdirs)

program test_stbmv_vector_reverse
  implicit none
  external :: stbmv
  external :: stbmv_bv
  integer :: nbdirs, n_test, seed_array(33), test_sizes(3), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing STBMV (Vector Reverse band, multi-size: n = 4)'
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
    real(4) :: alpha, alphab, beta, betab
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
    allocate(a(lda_val, n), ab(nbdirs, lda_val, n), x(n), xb(nbdirs, n), xb_seed(nbdirs, n))
    ! Initialize a as triangular band matrix (upper band storage)
    ! A(band_row, j) = full(i,j) with band_row = ksize+1+i-j, i = max(1,j-ksize)..j
    do j = 1, n
    do band_row = max(1, ksize+2-j), ksize+1
    call random_number(temp_real)
    a(band_row, j) = temp_real * 2.0 - 1.0  ! Scale to [-1,1]
    end do
    end do
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    ab = 0.0d0
    ! Seed for vector reverse: output adjoint xb is the seed per direction
    call random_number(xb)
    xb = xb * 2.0d0 - 1.0d0
    xb_seed = xb
    write(*,*) 'Testing STBMV (Vector Reverse band, n =', n, ')'
    call set_ISIZE2OFA(n)
    call stbmv_bv(uplo, trans, diag, nsize, ksize, a, ab, lda_val, x, xb, incx_val, nbdirs)
    call set_ISIZE2OFA(-1)
    write(*,*) 'Function calls completed successfully'
    call check_vjp_numerically_band_vec(n, nbdirs, lda_val, ksize, uplo, trans, diag, nsize, incx_val, a, ab, x, xb_seed, xb, passed)
    if (allocated(a)) deallocate(a)
    if (allocated(ab)) deallocate(ab)
    if (allocated(x)) deallocate(x)
    if (allocated(xb)) deallocate(xb)
    if (allocated(xb_seed)) deallocate(xb_seed)
  end subroutine run_test_for_size

  subroutine check_vjp_numerically_band_vec(n, nbdirs, lda_val, ksize, uplo, trans, diag, nsize, incx_val, a, ab, x, xb_seed, xb, passed)
    implicit none
    integer, intent(in) :: n, nbdirs, lda_val, ksize, nsize, incx_val
    character, intent(in) :: uplo, trans, diag
    real(4), intent(in) :: a(lda_val, n), ab(nbdirs, lda_val, n), x(n), xb_seed(nbdirs, n), xb(nbdirs, n)
    logical, intent(out) :: passed
    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_fd, vjp_ad, abs_error, abs_ref, err_bound, relative_error, max_re
    real(4), dimension(n) :: x_plus, x_minus, x_t, x_dir
    real(4), dimension(lda_val, n) :: a_t, a_dir
    real(4), dimension(:), allocatable :: temp_products
    integer :: i, j, band_row, n_products, k
    logical :: has_err
    has_err = .false.
    max_re = 0.0d0
    allocate(temp_products(n + n + (ksize+1)*n))
    do k = 1, nbdirs
    vjp_fd = 0.0d0
    ! Random direction for this k
    do j = 1, n
      do band_row = max(1, ksize+2-j), ksize+1
        call random_number(a_dir(band_row, j))
        a_dir(band_row, j) = a_dir(band_row, j) * 2.0d0 - 1.0d0
      end do
    end do
    call random_number(x_dir)
    x_dir = x_dir * 2.0d0 - 1.0d0
    ! Forward perturbation: f(a + h*a_dir, x + h*x_dir)
    a_t = a
    do j = 1, n
      do band_row = max(1, ksize+2-j), ksize+1
        a_t(band_row, j) = a(band_row, j) + h * a_dir(band_row, j)
      end do
    end do
    x_t = x + h * x_dir
    call stbmv(uplo, trans, diag, nsize, ksize, a_t, lda_val, x_t, incx_val)
    x_plus = x_t
    ! Backward perturbation: f(a - h*a_dir, x - h*x_dir)
    a_t = a
    do j = 1, n
      do band_row = max(1, ksize+2-j), ksize+1
        a_t(band_row, j) = a(band_row, j) - h * a_dir(band_row, j)
      end do
    end do
    x_t = x - h * x_dir
    call stbmv(uplo, trans, diag, nsize, ksize, a_t, lda_val, x_t, incx_val)
    x_minus = x_t
    n_products = n
    do i = 1, n
      temp_products(i) = xb_seed(k,i) * ((x_plus(i) - x_minus(i)) / (2.0d0 * h))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do
    ! VJP(AD) = direction^T @ adjoint
    vjp_ad = 0.0d0
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
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    abs_error = abs(vjp_fd - vjp_ad)
    abs_ref = abs(vjp_ad)
    err_bound = 2.0e-3 + 2.0e-3 * abs_ref
    if (abs_error > err_bound) has_err = .true.
    relative_error = abs_error / max(abs_ref, 1.0d-10)
    if (relative_error > max_re) max_re = relative_error
    end do
    deallocate(temp_products)
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Maximum relative error:', max_re
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
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
end program test_stbmv_vector_reverse