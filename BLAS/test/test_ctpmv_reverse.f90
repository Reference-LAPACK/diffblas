! Test program for CTPMV reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size outlined run_test_for_size(n) - TPMV/TPSV packed triangular

program test_ctpmv_reverse
  implicit none
  external :: ctpmv
  external :: ctpmv_b
  integer :: n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing CTPMV (multi-size: n = 4)'
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
    integer :: nsize, incx_val, npack
    complex(4), allocatable :: ap(:), apb(:), x(:), xb(:)
    complex(4), allocatable :: ap_orig(:), ap_plus(:), ap_minus(:), x_orig(:), x_plus(:), x_minus(:), xb_dir(:), apb_dir(:)
    integer :: ii
    real(4) :: tr, ti
    write(*,*) 'Testing CTPMV (n =', n, ')'
    uplo = 'U'
    trans = 'N'
    diag = 'N'
    nsize = n
    incx_val = 1
    npack = (n * (n + 1)) / 2
    allocate(ap(npack), apb(npack), x(n), xb(n))
    allocate(ap_orig(npack), ap_plus(npack), ap_minus(npack), x_orig(n), x_plus(n), x_minus(n), xb_dir(n), apb_dir(npack))
    do ii = 1, npack
      call random_number(tr)
      call random_number(ti)
      ap(ii) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(ap))
    end do
    do ii = 1, n
      call random_number(tr)
      call random_number(ti)
      x(ii) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(x))
    end do
    ap_orig = ap
    x_orig = x
    do ii = 1, n
      call random_number(tr)
      call random_number(ti)
      xb(ii) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(xb))
    end do
    do ii = 1, npack
      call random_number(tr)
      call random_number(ti)
      apb(ii) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(apb))
    end do
    xb_dir = xb
    apb_dir = apb
    call set_ISIZE1OFAp(npack)
    call ctpmv_b(uplo, trans, diag, nsize, ap, apb, x, xb, incx_val)
    call set_ISIZE1OFAp(-1)
    call check_vjp_numerically(n, npack, uplo, trans, diag, nsize, incx_val, ap_orig, x_orig, xb_dir, apb_dir, xb, apb, passed)
    deallocate(ap, apb, x, xb, ap_orig, ap_plus, ap_minus, x_orig, x_plus, x_minus, xb_dir, apb_dir)
  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, npack, uplo, trans, diag, nsize, incx_val, ap_orig, x_orig, xb_dir, apb_dir, xb_adj, apb_adj, passed)
    implicit none
    integer, intent(in) :: n, npack, nsize, incx_val
    character, intent(in) :: uplo, trans, diag
    complex(4), intent(in) :: ap_orig(npack), x_orig(n), xb_dir(n), apb_dir(npack), xb_adj(n), apb_adj(npack)
    logical, intent(out) :: passed
    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_fd, vjp_ad, abs_error, abs_reference, error_bound, relative_error
    complex(4) :: ap_t(npack), x_t(n), x_plus(n), x_minus(n)
    integer :: i, j
    vjp_fd = 0.0d0
    do i = 1, n
      x_t = x_orig
      ap_t = ap_orig
      x_t(i) = x_orig(i) + h * xb_dir(i)
      call ctpmv(uplo, trans, diag, nsize, ap_t, x_t, incx_val)
      x_plus = x_t
      x_t = x_orig
      ap_t = ap_orig
      x_t(i) = x_orig(i) - h * xb_dir(i)
      call ctpmv(uplo, trans, diag, nsize, ap_t, x_t, incx_val)
      x_minus = x_t
      do j = 1, n
        vjp_fd = vjp_fd + real(conjg(xb_dir(j)) * (x_plus(j) - x_minus(j)) / (2.0d0 * h))
      end do
    end do
    do i = 1, npack
      x_t = x_orig
      ap_t = ap_orig
      ap_t(i) = ap_orig(i) + h * apb_dir(i)
      call ctpmv(uplo, trans, diag, nsize, ap_t, x_t, incx_val)
      x_plus = x_t
      x_t = x_orig
      ap_t = ap_orig
      ap_t(i) = ap_orig(i) - h * apb_dir(i)
      call ctpmv(uplo, trans, diag, nsize, ap_t, x_t, incx_val)
      x_minus = x_t
      do j = 1, n
        vjp_fd = vjp_fd + real(conjg(xb_dir(j)) * (x_plus(j) - x_minus(j)) / (2.0d0 * h))
      end do
    end do
    vjp_ad = 0.0d0
    do i = 1, n
      vjp_ad = vjp_ad + real(conjg(xb_dir(i)) * xb_adj(i))
    end do
    do i = 1, npack
      vjp_ad = vjp_ad + real(conjg(apb_dir(i)) * apb_adj(i))
    end do
    abs_error = abs(vjp_fd - vjp_ad)
    abs_reference = abs(vjp_ad)
    error_bound = 1.0e-3 + 1.0e-3 * abs_reference
    relative_error = 0.0d0
    if (abs_reference > 1.0d-10) then
      relative_error = abs_error / abs_reference
    end if
    write(*,*) 'Function calls completed successfully'
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
end program test_ctpmv_reverse