! Test program for STPMV reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size outlined run_test_for_size(n) - TPMV/TPSV packed triangular

program test_stpmv_reverse
  implicit none
  external :: stpmv
  external :: stpmv_b
  integer :: n_test, seed_array(33), test_sizes(3), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing STPMV (multi-size: n = 4)'
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
    integer :: nsize, incx_val, npack
    real(4), allocatable :: ap(:), apb(:), x(:), xb(:)
    real(4), allocatable :: ap_orig(:), ap_plus(:), ap_minus(:), x_orig(:), x_plus(:), x_minus(:), xb_dir(:), apb_dir(:)
    integer :: ii
    write(*,*) 'Testing STPMV (n =', n, ')'
    uplo = 'U'
    trans = 'N'
    diag = 'N'
    nsize = n
    incx_val = 1
    npack = (n * (n + 1)) / 2
    allocate(ap(npack), apb(npack), x(n), xb(n))
    allocate(ap_orig(npack), ap_plus(npack), ap_minus(npack), x_orig(n), x_plus(n), x_minus(n), xb_dir(n), apb_dir(npack))
    call random_number(ap)
    ap = ap * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    ap_orig = ap
    x_orig = x
    call random_number(xb)
    xb = xb * 2.0d0 - 1.0d0
    call random_number(apb)
    apb = apb * 2.0d0 - 1.0d0
    xb_dir = xb
    apb_dir = apb
    call set_ISIZE1OFAp(npack)
    call stpmv_b(uplo, trans, diag, nsize, ap, apb, x, xb, incx_val)
    call set_ISIZE1OFAp(-1)
    call check_vjp_numerically(n, npack, uplo, trans, diag, nsize, incx_val, ap_orig, x_orig, xb_dir, apb_dir, xb, apb, passed)
    deallocate(ap, apb, x, xb, ap_orig, ap_plus, ap_minus, x_orig, x_plus, x_minus, xb_dir, apb_dir)
  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, npack, uplo, trans, diag, nsize, incx_val, ap_orig, x_orig, xb_dir, apb_dir, xb_adj, apb_adj, passed)
    implicit none
    integer, intent(in) :: n, npack, nsize, incx_val
    character, intent(in) :: uplo, trans, diag
    real(4), intent(in) :: ap_orig(npack), x_orig(n), xb_dir(n), apb_dir(npack), xb_adj(n), apb_adj(npack)
    logical, intent(out) :: passed
    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_fd, vjp_ad, abs_error, abs_reference, error_bound, relative_error
    real(4) :: ap_t(npack), x_t(n), x_plus(n), x_minus(n)
    integer :: i, j
    vjp_fd = 0.0d0
    do i = 1, n
      x_t = x_orig
      ap_t = ap_orig
      x_t(i) = x_orig(i) + h * xb_dir(i)
      call stpmv(uplo, trans, diag, nsize, ap_t, x_t, incx_val)
      x_plus = x_t
      x_t = x_orig
      ap_t = ap_orig
      x_t(i) = x_orig(i) - h * xb_dir(i)
      call stpmv(uplo, trans, diag, nsize, ap_t, x_t, incx_val)
      x_minus = x_t
      do j = 1, n
        vjp_fd = vjp_fd + xb_dir(j) * (x_plus(j) - x_minus(j)) / (2.0d0 * h)
      end do
    end do
    do i = 1, npack
      x_t = x_orig
      ap_t = ap_orig
      ap_t(i) = ap_orig(i) + h * apb_dir(i)
      call stpmv(uplo, trans, diag, nsize, ap_t, x_t, incx_val)
      x_plus = x_t
      x_t = x_orig
      ap_t = ap_orig
      ap_t(i) = ap_orig(i) - h * apb_dir(i)
      call stpmv(uplo, trans, diag, nsize, ap_t, x_t, incx_val)
      x_minus = x_t
      do j = 1, n
        vjp_fd = vjp_fd + xb_dir(j) * (x_plus(j) - x_minus(j)) / (2.0d0 * h)
      end do
    end do
    vjp_ad = 0.0d0
    do i = 1, n
      vjp_ad = vjp_ad + xb_dir(i) * xb_adj(i)
    end do
    do i = 1, npack
      vjp_ad = vjp_ad + apb_dir(i) * apb_adj(i)
    end do
    abs_error = abs(vjp_fd - vjp_ad)
    abs_reference = abs(vjp_ad)
    error_bound = 2.0e-3 + 2.0e-3 * abs_reference
    relative_error = 0.0d0
    if (abs_reference > 1.0d-10) then
      relative_error = abs_error / abs_reference
    end if
    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Maximum relative error:', relative_error
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
    passed = abs_error <= error_bound
    if (.not. passed) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
  end subroutine check_vjp_numerically
end program test_stpmv_reverse