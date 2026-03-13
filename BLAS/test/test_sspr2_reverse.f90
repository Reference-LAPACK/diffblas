! Test program for SSPR2 reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size outlined run_test_for_size(n) - SPR/SPR2 packed (declarations in subroutines)

program test_sspr2_reverse
  implicit none
  external :: sspr2
  external :: sspr2_b
  integer :: n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing SSPR2 (multi-size: n = 4)'
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
    character :: uplo
    integer :: nsize, incx_val, incy_val, npack
    real(4) :: alpha, alphab
    real(4), dimension(n) :: x, xb
    real(4), allocatable :: ap(:), apb(:)
    real(4) :: alpha_orig
    real(4), dimension(n) :: x_orig
    real(4), allocatable :: ap_orig(:), ap_plus(:), ap_minus(:), apb_orig(:)
    real(4), dimension(n) :: y, yb, y_orig
    integer :: ii
    uplo = 'U'
    nsize = n
    incx_val = 1
    incy_val = 1
    npack = (n * (n + 1)) / 2
    allocate(ap(npack), apb(npack), ap_orig(npack), ap_plus(npack), ap_minus(npack), apb_orig(npack))
    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    call random_number(y)
    y = y * 2.0d0 - 1.0d0
    call random_number(ap)
    ap = ap * 2.0d0 - 1.0d0
    alpha_orig = alpha
    x_orig = x
    ap_orig = ap
    y_orig = y
    call random_number(apb)
    apb = apb * 2.0d0 - 1.0d0
    apb_orig = apb
    alphab = 0.0d0
    xb = 0.0d0
    write(*,*) 'Testing SSPR2 (n =', n, ')'
    call set_ISIZE1OFX(n)
    call set_ISIZE1OFY(n)
    call sspr2_b(uplo, nsize, alpha, alphab, x, xb, incx_val, y, yb, incy_val, ap, apb)
    call set_ISIZE1OFX(-1)
    call set_ISIZE1OFY(-1)
    write(*,*) 'Function calls completed successfully'
    call check_vjp_numerically(n, npack, uplo, nsize, incx_val, incy_val, alpha_orig, x_orig, ap_orig, apb_orig, alphab, xb, apb, passed, y_orig, yb)
    deallocate(ap, apb, ap_orig, ap_plus, ap_minus, apb_orig)
  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, npack, uplo, nsize, incx_val, incy_val, alpha_orig, x_orig, ap_orig, apb_orig, alphab, xb, apb, passed, y_orig, yb)
    implicit none
    integer, intent(in) :: n, npack
    character, intent(in) :: uplo
    integer, intent(in) :: nsize, incx_val, incy_val
    real(4), intent(in) :: alpha_orig
    real(4), intent(in) :: x_orig(n), ap_orig(npack), apb_orig(npack)
    real(4), intent(in) :: alphab, xb(n), apb(npack)
    logical, intent(out) :: passed
    real(4), intent(in), optional :: y_orig(n), yb(n)
    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_fd, vjp_ad, abs_error, abs_reference, error_bound, relative_error
    real(4) :: alpha_dir
    real(4), dimension(n) :: x_dir, x_t
    real(4), dimension(npack) :: ap_dir, ap_t, ap_plus, ap_minus, ap_central_diff
    real(4), dimension(npack) :: temp_products
    real(4), dimension(n) :: y_dir, y_t
    real(4) :: alpha_t
    integer :: i, n_products
    call random_number(alpha_dir)
    alpha_dir = alpha_dir * 2.0d0 - 1.0d0
    call random_number(x_dir)
    x_dir = x_dir * 2.0d0 - 1.0d0
    if (present(y_orig)) call random_number(y_dir)
    if (present(y_orig)) y_dir = y_dir * 2.0d0 - 1.0d0
    call random_number(ap_dir)
    ap_dir = ap_dir * 2.0d0 - 1.0d0
    alpha_t = alpha_orig + h * alpha_dir
    x_t = x_orig + h * x_dir
    ap_t = ap_orig + h * ap_dir
    if (present(y_orig)) y_t = y_orig + h * y_dir
    if (present(y_orig)) then
      call sspr2(uplo, nsize, alpha_t, x_t, incx_val, y_t, incy_val, ap_t)
    else
      call sspr2(uplo, nsize, alpha_t, x_t, incx_val, ap_t)
    end if
    ap_plus = ap_t
    alpha_t = alpha_orig - h * alpha_dir
    x_t = x_orig - h * x_dir
    ap_t = ap_orig - h * ap_dir
    if (present(y_orig)) y_t = y_orig - h * y_dir
    if (present(y_orig)) then
      call sspr2(uplo, nsize, alpha_t, x_t, incx_val, y_t, incy_val, ap_t)
    else
      call sspr2(uplo, nsize, alpha_t, x_t, incx_val, ap_t)
    end if
    ap_minus = ap_t
    ap_central_diff = (ap_plus - ap_minus) / (2.0d0 * h)
    vjp_fd = 0.0d0
    n_products = npack
    do i = 1, n_products
      temp_products(i) = apb_orig(i) * ap_central_diff(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do
    vjp_ad = alpha_dir * alphab
    n_products = n
    do i = 1, n
      temp_products(i) = x_dir(i) * xb(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    n_products = npack
    do i = 1, n_products
      temp_products(i) = ap_dir(i) * apb(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    if (present(y_orig)) then
      n_products = n
      do i = 1, n
        temp_products(i) = y_dir(i) * yb(i)
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_ad = vjp_ad + temp_products(i)
      end do
    end if
    abs_error = abs(vjp_fd - vjp_ad)
    abs_reference = abs(vjp_ad)
    relative_error = 0.0d0
    if (abs_reference > 1.0d-10) relative_error = abs_error / abs_reference
    error_bound = 1.0e-3 + 1.0e-3 * abs_reference
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Maximum relative error:', relative_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    passed = abs_error <= error_bound
    if (.not. passed) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
  end subroutine check_vjp_numerically

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
end program test_sspr2_reverse