! Test program for DSPMV vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Multi-size outlined - SPMV vector reverse

program test_dspmv_vector_reverse
  implicit none
  external :: dspmv
  external :: dspmv_bv
  integer :: nbdirs, n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing DSPMV (Vector Reverse, multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    nbdirs = n_test
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
    character :: uplo
    integer :: nsize, incx_val, incy_val, npack, k
    real(8) :: alpha, alphab(nbdirs), beta, betab(nbdirs)
    real(8), dimension(n) :: x, y, y_orig
    real(8), dimension(nbdirs,n) :: xb, yb, yb_seed
    real(8), dimension(:), allocatable :: ap
    real(8), dimension(:,:), allocatable :: apb
    real(8), dimension(:), allocatable :: ap_orig, ap_t, x_orig
    real(8), dimension(n) :: y_plus, y_minus
    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_fd, vjp_ad, re, err_bnd
    integer :: ii
    write(*,*) 'Testing DSPMV (Vector Reverse, n =', n, ')'
    uplo = 'U'
    nsize = n
    incx_val = 1
    incy_val = 1
    npack = (n * (n + 1)) / 2
    allocate(ap(npack), apb(nbdirs, npack), ap_orig(npack), ap_t(npack), x_orig(n))
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
    ap_orig = ap
    x_orig = x
    y_orig = y
    yb_seed = yb
    alphab = 0.0d0
    betab = 0.0d0
    xb = 0.0d0
    apb = 0.0d0
    call set_ISIZE1OFAp(npack)
    call set_ISIZE1OFX(n)
    call dspmv_bv(uplo, nsize, alpha, alphab, ap, apb, x, xb, incx_val, beta, betab, y, yb, incy_val, nbdirs)
    call set_ISIZE1OFAp(-1)
    call set_ISIZE1OFX(-1)
    re = 0.0d0
    do k = 1, nbdirs
      y_plus = y_orig + h * yb_seed(k,:)
      ap_t = ap_orig + h * apb(k,:)
      call dspmv(uplo, nsize, alpha + h*alphab(k), ap_t, x_orig + h*xb(k,:), incx_val, beta + h*betab(k), y_plus, incy_val)
      y_minus = y_orig - h * yb_seed(k,:)
      ap_t = ap_orig - h * apb(k,:)
      call dspmv(uplo, nsize, alpha - h*alphab(k), ap_t, x_orig - h*xb(k,:), incx_val, beta - h*betab(k), y_minus, incy_val)
      vjp_fd = sum(yb_seed(k,:) * (y_plus - y_minus)) / (2.0d0 * h)
      vjp_ad = alphab(k)*alphab(k) + betab(k)*betab(k) + sum(apb(k,:)*apb(k,:)) + sum(xb(k,:)*xb(k,:)) + sum(yb_seed(k,:)*yb(k,:))
      re = max(re, abs(vjp_fd - vjp_ad))
    end do
    err_bnd = 1.0e-5 + 1.0e-5 * 1.0d0
    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Maximum relative error:', re
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = (re <= err_bnd)
    if (.not. passed) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
    deallocate(ap, apb, ap_orig, ap_t, x_orig)
  end subroutine run_test_for_size
end program test_dspmv_vector_reverse