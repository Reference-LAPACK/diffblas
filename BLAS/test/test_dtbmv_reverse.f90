! Test program for DTBMV reverse mode (adjoint) - BLAS2 band
! Generated automatically by run_tapenade_blas.py
! Multi-size outlined run_test_for_size(n) - band (declarations in subroutines)

program test_dtbmv_reverse
  implicit none
  external :: dtbmv
  external :: dtbmv_b
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
    real(8) :: alpha, alphab
    real(8), dimension(:,:), allocatable :: a, ab
    real(8), dimension(:), allocatable :: x, xb
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
    allocate(a(lda_val, n), ab(lda_val, n), x(n), xb(n))
    ! Initialize a as triangular band matrix (upper band storage)
    ! A(band_row, j) = full(i,j) with band_row = ksize+1+i-j, i = max(1,j-ksize)..j
    do j = 1, n
    do band_row = max(1, ksize+2-j), ksize+1
    call random_number(temp_real)
    a(band_row, j) = temp_real * 2.0 - 1.0  ! Scale to [-1,1]
    end do
    end do
    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    alphab = 0.0d0
    xb = 0.0d0
    ab = 0.0d0
    write(*,*) 'Testing DTBMV (n =', n, ')'
    call set_ISIZE2OFA(lda_val)
    call dtbmv_b(uplo, trans, diag, nsize, ksize, a, ab, lda_val, x, xb, incx_val)
    call set_ISIZE2OFA(-1)
    call check_vjp_numerically_band(n, lda_val, ksize, uplo, trans, diag, nsize, incx_val, a, ab, x, xb, passed)
    deallocate(a, ab, x, xb)
  end subroutine run_test_for_size

  subroutine check_vjp_numerically_band(n, lda_val, ksize, uplo, trans, diag, nsize, incx_val, a, ab, x, xb, passed)
    implicit none
    integer, intent(in) :: n, lda_val, ksize, nsize, incx_val
    character, intent(in) :: uplo, trans, diag
    real(8), intent(in) :: a(lda_val, n), ab(lda_val, n), x(n), xb(n)
    logical, intent(out) :: passed
    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_fd, vjp_ad, abs_error, abs_ref, err_bound
    real(8), dimension(n) :: x_plus, x_minus, x_t
    real(8), dimension(lda_val, n) :: a_t
    real(8), dimension(:), allocatable :: temp_products
    integer :: i, j, band_row, n_products
    allocate(temp_products(n + (ksize+1)*n))
    vjp_fd = 0.0d0
    a_t = a + h * ab
    x_t = x + h * xb
    call dtbmv(uplo, trans, diag, nsize, ksize, a_t, lda_val, x_t, incx_val)
    x_plus = x_t
    a_t = a - h * ab
    x_t = x - h * xb
    call dtbmv(uplo, trans, diag, nsize, ksize, a_t, lda_val, x_t, incx_val)
    x_minus = x_t
    n_products = n
    do i = 1, n
      temp_products(i) = xb(i) * ((x_plus(i) - x_minus(i)) / (2.0d0 * h))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do
    vjp_ad = 0.0d0
    do i = 1, n
      vjp_ad = vjp_ad + xb(i) * xb(i)
    end do
    n_products = 0
    do j = 1, n
      do band_row = max(1, ksize+2-j), ksize+1
        n_products = n_products + 1
        temp_products(n_products) = ab(band_row,j) * ab(band_row,j)
      end do
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    deallocate(temp_products)
    abs_error = abs(vjp_fd - vjp_ad)
    abs_ref = abs(vjp_ad)
    err_bound = 1.0e-5 + 1.0e-5 * abs_ref
    passed = abs_error <= err_bound
    if (.not. passed) write(*,*) 'FAIL: Band VJP error'
    if (passed) write(*,*) 'PASS: Band VJP within tolerance'
  end subroutine check_vjp_numerically_band
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
end program test_dtbmv_reverse