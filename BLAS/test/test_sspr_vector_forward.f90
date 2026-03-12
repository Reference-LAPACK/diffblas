! Test program for SSPR vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Multi-size outlined run_test_for_size(n) - SPR/SPR2 packed

program test_sspr_vector_forward
  implicit none
  external :: sspr
  external :: sspr_dv
  integer :: nbdirs, n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing SSPR (Vector Forward, multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    nbdirs = test_sizes(i)
    call run_test_for_size(n_test, passed, nbdirs)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) write(*,*) 'PASS: Vector forward - all sizes OK'
  if (.not. all_passed) write(*,*) 'FAIL: Vector forward - derivative errors'
contains
  subroutine run_test_for_size(n, passed, nbdirs)
    implicit none
    integer, intent(in) :: n, nbdirs
    logical, intent(out) :: passed
    character :: uplo
    integer :: nsize, incx_val, incy_val, npack
    real(4) :: alpha
    real(4), dimension(n) :: x
    real(4), allocatable :: ap(:), ap_orig(:)
    real(4), dimension(nbdirs) :: alpha_dv
    real(4), dimension(nbdirs,n) :: x_dv
    real(4), allocatable :: ap_dv(:,:), ap_dv_seed(:,:)
    integer :: idir, ii
    real(4) :: tr, ti
    uplo = 'U'
    nsize = n
    incx_val = 1
    incy_val = 1
    npack = (n * (n + 1)) / 2
    allocate(ap(npack), ap_orig(npack), ap_dv(nbdirs, npack), ap_dv_seed(nbdirs, npack))
    call random_number(tr)
    alpha = tr * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    do idir = 1, nbdirs
      call random_number(tr)
      alpha_dv(idir) = tr * 2.0d0 - 1.0d0
    end do
    do idir = 1, nbdirs
      call random_number(x_dv(idir,:))
      x_dv(idir,:) = x_dv(idir,:) * 2.0d0 - 1.0d0
    end do
    call random_number(ap)
    ap = ap * 2.0d0 - 1.0d0
    do idir = 1, nbdirs
      call random_number(ap_dv(idir,:))
      ap_dv(idir,:) = ap_dv(idir,:) * 2.0d0 - 1.0d0
    end do

    write(*,*) 'Testing SSPR (Vector Forward, n =', n, ')'
    ap_orig = ap
    ap_dv_seed = ap_dv
    call sspr_dv(uplo, nsize, alpha, alpha_dv, x, x_dv, incx_val, ap, ap_dv, nbdirs)
    call check_derivatives_numerically(n, npack, nbdirs, uplo, nsize, incx_val, alpha, alpha_dv, x, x_dv, ap_orig, ap_dv, ap_dv_seed, passed)
    deallocate(ap, ap_orig, ap_dv, ap_dv_seed)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, npack, nbdirs, uplo, nsize, incx_val, alpha, alpha_dv, x, x_dv, ap_orig, ap_dv, ap_dv_seed, passed)
    integer, intent(in) :: n, npack, nbdirs
    character, intent(in) :: uplo
    integer, intent(in) :: nsize, incx_val
    real(4), intent(in) :: alpha
    real(4), intent(in) :: alpha_dv(nbdirs), x(n), x_dv(nbdirs,n)
    real(4), intent(in) :: ap_orig(npack), ap_dv(nbdirs,npack), ap_dv_seed(nbdirs,npack)
    logical, intent(out) :: passed
    real(4), parameter :: h = 1.0e-3
    real(4) :: abs_error, abs_ref, err_bound
    real(4), dimension(npack) :: ap_fwd, ap_bwd, ap_t
    real(4) :: alpha_t
    real(4), dimension(n) :: x_t
    integer :: idir, ii
    logical :: has_err
    has_err = .false.
    do idir = 1, nbdirs
      alpha_t = alpha + h * alpha_dv(idir)
      x_t = x + h * x_dv(idir,:)
      ap_t = ap_orig + h * ap_dv_seed(idir,:)
      call sspr(uplo, nsize, alpha_t, x_t, incx_val, ap_t)
      ap_fwd = ap_t
      alpha_t = alpha - h * alpha_dv(idir)
      x_t = x - h * x_dv(idir,:)
      ap_t = ap_orig - h * ap_dv_seed(idir,:)
      call sspr(uplo, nsize, alpha_t, x_t, incx_val, ap_t)
      ap_bwd = ap_t
      do ii = 1, min(3, npack)
        abs_error = abs((ap_fwd(ii) - ap_bwd(ii)) / (2.0e0 * h) - ap_dv(idir,ii))
        abs_ref = abs(ap_dv(idir,ii))
        err_bound = 1.0e-3 + 1.0e-3 * abs_ref
        if (abs_error > err_bound) has_err = .true.
      end do
    end do
    passed = .not. has_err
    if (has_err) write(*,*) 'FAIL: SPR/SPR2 vector derivatives'
    if (.not. has_err) write(*,*) 'PASS: SPR/SPR2 vector derivatives'
  end subroutine check_derivatives_numerically

end program test_sspr_vector_forward