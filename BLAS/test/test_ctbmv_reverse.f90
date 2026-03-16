! Test program for CTBMV reverse mode (adjoint) - BLAS2 band
! Generated automatically by run_tapenade_blas.py
! Multi-size outlined run_test_for_size(n) - band (declarations in subroutines)

program test_ctbmv_reverse
  implicit none
  external :: ctbmv
  external :: ctbmv_b
  integer :: n_test, seed_array(33), test_sizes(3), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing CTBMV (multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 3
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
    complex(4) :: alpha, alphab
    complex(4), dimension(:,:), allocatable :: a, ab
    complex(4), dimension(:), allocatable :: x, xb
    complex(4), dimension(:), allocatable :: xb_seed
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
    allocate(a(lda_val, n), ab(lda_val, n), x(n), xb(n))
    allocate(xb_seed(n))
    ! Initialize a as triangular band matrix (upper band storage)
    do j = 1, n
    do band_row = max(1, ksize+2-j), ksize+1
    call random_number(temp_real)
    call random_number(temp_imag)
    a(band_row, j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
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
    alphab = 0.0d0
    ab = 0.0d0
    ! Seed for reverse mode: output adjoint xb is the seed (d(scalar)/d(x))
    do j = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      xb(j) = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(xb))
    end do
    xb_seed = xb
    write(*,*) 'Testing CTBMV (n =', n, ')'
    call set_ISIZE2OFA(lda_val)
    call ctbmv_b(uplo, trans, diag, nsize, ksize, a, ab, lda_val, x, xb, incx_val)
    call set_ISIZE2OFA(-1)
    write(*,*) 'Function calls completed successfully'
    call check_vjp_numerically_band(n, lda_val, ksize, uplo, trans, diag, nsize, incx_val, a, ab, x, xb_seed, xb, passed)
    deallocate(a, ab, x, xb)
    deallocate(xb_seed)
  end subroutine run_test_for_size

  subroutine check_vjp_numerically_band(n, lda_val, ksize, uplo, trans, diag, nsize, incx_val, a, ab, x, xb_seed, xb, passed)
    implicit none
    integer, intent(in) :: n, lda_val, ksize, nsize, incx_val
    character, intent(in) :: uplo, trans, diag
    complex(4), intent(in) :: a(lda_val, n), ab(lda_val, n), x(n), xb_seed(n), xb(n)
    logical, intent(out) :: passed
    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_fd, vjp_ad, abs_error, abs_ref, err_bound, relative_error
    complex(4), dimension(n) :: x_plus, x_minus, x_t, x_dir
    complex(4), dimension(lda_val, n) :: a_t, a_dir
    real(4), dimension(:), allocatable :: temp_products
    real(kind(0.0d0)) :: tr, ti
    integer :: i, j, band_row, n_products
    allocate(temp_products(n + n + (ksize+1)*n))
    ! Random direction for FD (direction^T @ adjoint)
    do j = 1, n
      do band_row = max(1, ksize+2-j), ksize+1
        call random_number(tr)
        call random_number(ti)
        a_dir(band_row, j) = cmplx(tr*2.0d0-1.0d0, ti*2.0d0-1.0d0, kind=kind(a_dir))
      end do
    end do
    do i = 1, n
      call random_number(tr)
      call random_number(ti)
      x_dir(i) = cmplx(tr*2.0d0-1.0d0, ti*2.0d0-1.0d0, kind=kind(x_dir))
    end do
    ! Forward perturbation: f(a + h*a_dir, x + h*x_dir)
    a_t = a
    do j = 1, n
      do band_row = max(1, ksize+2-j), ksize+1
        a_t(band_row, j) = a(band_row, j) + h * a_dir(band_row, j)
      end do
    end do
    x_t = x + h * x_dir
    call ctbmv(uplo, trans, diag, nsize, ksize, a_t, lda_val, x_t, incx_val)
    x_plus = x_t
    ! Backward perturbation: f(a - h*a_dir, x - h*x_dir)
    a_t = a
    do j = 1, n
      do band_row = max(1, ksize+2-j), ksize+1
        a_t(band_row, j) = a(band_row, j) - h * a_dir(band_row, j)
      end do
    end do
    x_t = x - h * x_dir
    call ctbmv(uplo, trans, diag, nsize, ksize, a_t, lda_val, x_t, incx_val)
    x_minus = x_t
    ! VJP(FD) = xb_seed^T @ (x_plus-x_minus)/(2h)
    vjp_fd = 0.0d0
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(xb_seed(i)) * ((x_plus(i) - x_minus(i)) / (2.0d0 * h)))
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
        temp_products(n_products) = real(conjg(a_dir(band_row,j)) * ab(band_row,j))
      end do
    end do
    do i = 1, n
      n_products = n_products + 1
      temp_products(n_products) = real(conjg(x_dir(i)) * xb(i))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    deallocate(temp_products)
    abs_error = abs(vjp_fd - vjp_ad)
    abs_ref = abs(vjp_ad)
    err_bound = 1.0e-3 + 1.0e-3 * abs_ref
    relative_error = 0.0d0
    if (abs_ref > 1.0d-10) relative_error = abs_error / abs_ref
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Maximum relative error:', relative_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    passed = abs_error <= err_bound
    if (.not. passed) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
  end subroutine check_vjp_numerically_band
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
end program test_ctbmv_reverse