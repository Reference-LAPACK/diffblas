! Test program for DSPR vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Multi-size outlined run_test_for_size(n) - SPR/SPR2 packed

program test_dspr_vector_forward
  implicit none
  external :: dspr
  external :: dspr_dv
  integer :: nbdirs, n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing DSPR (Vector Forward, multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 1
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
    character :: uplo
    integer :: nsize, incx_val, incy_val, npack
    real(8) :: alpha
    real(8), dimension(n) :: x
    real(8), allocatable :: ap(:), ap_orig(:)
    real(8), dimension(nbdirs) :: alpha_dv
    real(8), dimension(nbdirs,n) :: x_dv
    real(8), allocatable :: ap_dv(:,:), ap_dv_seed(:,:)
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

    write(*,*) 'Testing DSPR (Vector Forward, n =', n, ')'
    ap_orig = ap
    ap_dv_seed = ap_dv
    call dspr_dv(uplo, nsize, alpha, alpha_dv, x, x_dv, incx_val, ap, ap_dv, nbdirs)
    call check_derivatives_numerically(n, npack, nbdirs, uplo, nsize, incx_val, alpha, alpha_dv, x, x_dv, ap_orig, ap_dv, ap_dv_seed, passed)
    deallocate(ap, ap_orig, ap_dv, ap_dv_seed)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, npack, nbdirs, uplo, nsize, incx_val, alpha, alpha_dv, x, x_dv, ap_orig, ap_dv, ap_dv_seed, passed)
    integer, intent(in) :: n, npack, nbdirs
    character, intent(in) :: uplo
    integer, intent(in) :: nsize, incx_val
    real(8), intent(in) :: alpha
    real(8), intent(in) :: alpha_dv(nbdirs), x(n), x_dv(nbdirs,n)
    real(8), intent(in) :: ap_orig(npack), ap_dv(nbdirs,npack), ap_dv_seed(nbdirs,npack)
    logical, intent(out) :: passed
    real(8), parameter :: h = 1.0e-7
    real(8) :: abs_error, abs_ref, err_bound, max_error, relative_error
    real(8), dimension(npack) :: ap_fwd, ap_bwd, ap_t
    real(8) :: alpha_t
    real(8), dimension(n) :: x_t
    integer :: idir, ii
    logical :: has_err
    has_err = .false.
    max_error = 0.0e0
    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    do idir = 1, nbdirs
      alpha_t = alpha + h * alpha_dv(idir)
      x_t = x + h * x_dv(idir,:)
      ap_t = ap_orig + h * ap_dv_seed(idir,:)
      call dspr(uplo, nsize, alpha_t, x_t, incx_val, ap_t)
      ap_fwd = ap_t
      alpha_t = alpha - h * alpha_dv(idir)
      x_t = x - h * x_dv(idir,:)
      ap_t = ap_orig - h * ap_dv_seed(idir,:)
      call dspr(uplo, nsize, alpha_t, x_t, incx_val, ap_t)
      ap_bwd = ap_t
      do ii = 1, min(3, npack)
        abs_error = abs((ap_fwd(ii) - ap_bwd(ii)) / (2.0e0 * h) - ap_dv(idir,ii))
        abs_ref = abs(ap_dv(idir,ii))
        err_bound = 1.0e-5 + 1.0e-5 * abs_ref
        if (abs_error > err_bound) has_err = .true.
        relative_error = abs_error / max(abs_ref, 1.0e-10)
        if (relative_error > max_error) max_error = relative_error
      end do
    end do
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_err
    if (has_err) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
  end subroutine check_derivatives_numerically

end program test_dspr_vector_forward