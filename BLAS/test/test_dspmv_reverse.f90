! Test program for DSPMV reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size outlined - SPMV (symmetric packed matrix-vector)

program test_dspmv_reverse
  implicit none
  external :: dspmv
  external :: dspmv_b
  integer :: n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing DSPMV (multi-size: n = 4)'
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
    real(8) :: alpha, alphab, beta, betab, alpha_orig, beta_orig
    real(8), dimension(n) :: x, xb, y, yb, y_orig, yb_orig
    real(8), dimension(:), allocatable :: ap, apb, ap_orig, x_orig
    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_fd, vjp_ad, re, err_bnd, max_error
    integer :: ii
    write(*,*) 'Testing DSPMV (n =', n, ')'
    uplo = 'U'
    nsize = n
    incx_val = 1
    incy_val = 1
    npack = (n * (n + 1)) / 2
    allocate(ap(npack), apb(npack), ap_orig(npack), x_orig(n))
    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(beta)
    beta = beta * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    call random_number(y)
    y = y * 2.0d0 - 1.0d0
    call random_number(ap)
    ap = ap * 2.0d0 - 1.0d0
    call random_number(yb)
    yb = yb * 2.0d0 - 1.0d0
    alpha_orig = alpha
    beta_orig = beta
    ap_orig = ap
    x_orig = x
    y_orig = y
    yb_orig = yb
    alphab = 0.0d0
    betab = 0.0d0
    xb = 0.0d0
    apb = 0.0d0
    call set_ISIZE1OFAp(npack)
    call set_ISIZE1OFX(n)
    call dspmv_b(uplo, nsize, alpha, alphab, ap, apb, x, xb, incx_val, beta, betab, y, yb, incy_val)
    call set_ISIZE1OFAp(-1)
    call set_ISIZE1OFX(-1)
    call check_vjp_spmv(n, npack, uplo, nsize, incx_val, incy_val, alpha_orig, alphab, ap_orig, apb, x_orig, xb, beta_orig, betab, y_orig, yb_orig, yb, passed)
    deallocate(ap, apb, ap_orig, x_orig)
  end subroutine run_test_for_size

  subroutine check_vjp_spmv(n, npack, uplo, nsize, incx_val, incy_val, alpha_orig, alphab, ap_orig, apb, x_orig, xb, beta_orig, betab, y_orig, yb_seed, yb, passed)
    implicit none
    integer, intent(in) :: n, npack, nsize, incx_val, incy_val
    character, intent(in) :: uplo
    real(8), intent(in) :: alpha_orig, beta_orig
    real(8), intent(in) :: ap_orig(npack), x_orig(n), y_orig(n)
    real(8), intent(in) :: alphab, betab, apb(npack), xb(n), yb_seed(n), yb(n)
    logical, intent(out) :: passed
    real(8) :: alpha_t, beta_t, ap_t(npack), x_t(n), y_t(n)
    real(8) :: vjp_fd, vjp_ad, re, err_bnd, relative_error
    real(8), parameter :: h = 1.0e-7
    integer :: i
    vjp_fd = 0.0d0
    vjp_ad = 0.0d0
    alpha_t = alpha_orig + h * alphab
    beta_t = beta_orig + h * betab
    ap_t = ap_orig + h * apb
    x_t = x_orig + h * xb
    y_t = y_orig + h * yb_seed
    call dspmv(uplo, nsize, alpha_t, ap_t, x_t, incx_val, beta_t, y_t, incy_val)
    vjp_fd = vjp_fd + sum(yb_seed * y_t)
    alpha_t = alpha_orig - h * alphab
    beta_t = beta_orig - h * betab
    ap_t = ap_orig - h * apb
    x_t = x_orig - h * xb
    y_t = y_orig - h * yb_seed
    call dspmv(uplo, nsize, alpha_t, ap_t, x_t, incx_val, beta_t, y_t, incy_val)
    vjp_fd = (vjp_fd - sum(yb_seed * y_t)) / (2.0d0 * h)
    vjp_ad = alphab*alphab + betab*betab + sum(apb*apb) + sum(xb*xb) + sum(yb_seed*yb)
    re = abs(vjp_fd - vjp_ad)
    err_bnd = 1.0e-5 + 1.0e-5 * abs(vjp_ad)
    relative_error = 0.0d0
    if (abs(vjp_ad) > 1.0d-10) relative_error = re / abs(vjp_ad)
    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Maximum relative error:', relative_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = (re <= err_bnd)
    if (.not. passed) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
  end subroutine check_vjp_spmv
end program test_dspmv_reverse