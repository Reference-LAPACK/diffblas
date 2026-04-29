! Test program for DSPR vector reverse mode differentiation
! Multi-size outlined run_test_for_size(n) - SPR/SPR2 packed

program test_dspr_vector_reverse
  implicit none
  external :: dspr
  external :: dspr_bv
  integer :: nbdirs, n_test, seed_array(33), test_sizes(3), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing DSPR (Vector Reverse, multi-size: n =', test_sizes(1), ')'
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
    integer, intent(in) :: n, nbdirs
    logical, intent(out) :: passed
    character :: uplo
    integer :: nsize, incx_val, incy_val, npack
    real(8) :: alpha
    real(8), dimension(n) :: x
    real(8), allocatable :: ap(:)
    real(8), dimension(nbdirs) :: alphab
    real(8), dimension(nbdirs,n) :: xb
    real(8), allocatable :: apb(:,:)
    real(8), allocatable :: apb_orig(:,:)
    integer :: k, ii
    real(4) :: tr, ti
    uplo = 'L'
    nsize = n
    incx_val = 1
    incy_val = 1
    npack = (n * (n + 1)) / 2
    allocate(ap(npack), apb(nbdirs, npack), apb_orig(nbdirs, npack))
    call random_number(tr)
    alpha = tr * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    call random_number(ap)
    ap = ap * 2.0d0 - 1.0d0
    do k = 1, nbdirs
      call random_number(apb(k,:))
      apb(k,:) = apb(k,:) * 2.0d0 - 1.0d0
    end do
    apb_orig = apb
    alphab = 0.0d0
    xb = 0.0d0
    write(*,*) 'Testing DSPR (Vector Reverse, n =', n, ')'
    call set_ISIZE1OFX(n)
    call dspr_bv(uplo, nsize, alpha, alphab, x, xb, incx_val, ap, apb, nbdirs)
    call set_ISIZE1OFX(-1)
    write(*,*) 'Function calls completed successfully'
    call check_vjp_spr_spr2(n, npack, nbdirs, uplo, nsize, incx_val, incy_val, alpha, x, ap, apb_orig, alphab, xb, apb, passed)
    deallocate(ap, apb, apb_orig)
  end subroutine run_test_for_size
  subroutine check_vjp_spr_spr2(n, npack, nbdirs, uplo, nsize, incx_val, incy_val, alpha, x, ap, apb_orig, alphab, xb, apb, passed, y, yb)
    integer, intent(in) :: n, npack, nbdirs
    character, intent(in) :: uplo
    integer, intent(in) :: nsize, incx_val, incy_val
    real(8), intent(in) :: alpha, x(n)
    real(8), intent(in) :: ap(npack)
    real(8), intent(in) :: apb_orig(nbdirs,npack)
    real(8), intent(in) :: alphab(nbdirs), xb(nbdirs,n)
    real(8), intent(in) :: apb(nbdirs,npack)
    logical, intent(out) :: passed
    real(8), intent(in), optional :: y(n), yb(nbdirs,n)
    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_fd, vjp_ad, re, err_bnd, max_re
    real(4) :: tr, ti
    real(8) :: alpha_dir
    real(8), dimension(n) :: x_dir, x_t
    real(8), dimension(npack) :: ap_dir, ap_t, ap_plus, ap_minus, ap_cdiff
    real(8), dimension(n) :: y_dir, y_t
    integer :: k, ii
    logical :: has_err
    has_err = .false.
    max_re = 0.0d0
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    do k = 1, nbdirs
      call random_number(tr)
      call random_number(ti)
      alpha_dir = tr * 2.0d0 - 1.0d0
      call random_number(x_dir)
      x_dir = x_dir * 2.0d0 - 1.0d0
      call random_number(ap_dir)
      ap_dir = ap_dir * 2.0d0 - 1.0d0
      ap_t = ap + h * ap_dir
      x_t = x + h * x_dir
      call dspr(uplo, nsize, alpha + h*alpha_dir, x_t, incx_val, ap_t)
      ap_plus = ap_t
      ap_t = ap - h * ap_dir
      x_t = x - h * x_dir
      call dspr(uplo, nsize, alpha - h*alpha_dir, x_t, incx_val, ap_t)
      ap_minus = ap_t
      ap_cdiff = (ap_plus - ap_minus) / (2.0e0 * h)
      vjp_fd = sum(apb_orig(k,:) * ap_cdiff)
      vjp_ad = alpha_dir * alphab(k)
      vjp_ad = vjp_ad + sum(x_dir*xb(k,:))
      vjp_ad = vjp_ad + sum(ap_dir*apb(k,:))
      re = abs(vjp_fd - vjp_ad)
      if (re > max_re) max_re = re
      err_bnd = 1.0e-5 + 1.0e-5 * abs(vjp_ad)
      if (re > err_bnd) has_err = .true.
    end do
    write(*,*) 'Maximum relative error:', max_re
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_err
    if (has_err) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
  end subroutine check_vjp_spr_spr2
end program test_dspr_vector_reverse