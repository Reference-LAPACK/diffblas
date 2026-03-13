! Test program for DTPMV vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Multi-size outlined run_test_for_size(n) - TPMV/TPSV packed triangular

program test_dtpmv_vector_forward
  implicit none
  external :: dtpmv
  external :: dtpmv_dv
  integer :: nbdirs, n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing DTPMV (Vector Forward, multi-size: n = 4)'
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
    character :: uplo, trans, diag
    integer :: nsize, incx_val, npack
    real(8), allocatable :: ap(:), x(:)
    real(8), allocatable :: ap_dv(:,:), x_dv(:,:)
    real(8), allocatable :: ap_orig(:), x_orig(:)
    real(8), allocatable :: ap_dv_seed(:,:), x_dv_seed(:,:)
    integer :: idir, ii
    real(4) :: tr, ti
    uplo = 'U'
    trans = 'N'
    diag = 'N'
    nsize = n
    incx_val = 1
    npack = (n * (n + 1)) / 2
    allocate(ap(npack), x(n), ap_dv(nbdirs, npack), x_dv(nbdirs, n))
    allocate(ap_orig(npack), x_orig(n), ap_dv_seed(nbdirs, npack), x_dv_seed(nbdirs, n))
    call random_number(ap)
    ap = ap * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    do idir = 1, nbdirs
      call random_number(ap_dv(idir,:))
      ap_dv(idir,:) = ap_dv(idir,:) * 2.0d0 - 1.0d0
      call random_number(x_dv(idir,:))
      x_dv(idir,:) = x_dv(idir,:) * 2.0d0 - 1.0d0
    end do

    write(*,*) 'Testing DTPMV (Vector Forward, n =', n, ')'
    ap_orig = ap
    x_orig = x
    ap_dv_seed = ap_dv
    x_dv_seed = x_dv
    call dtpmv_dv(uplo, trans, diag, nsize, ap, ap_dv, x, x_dv, incx_val, nbdirs)
    write(*,*) 'Function calls completed successfully'
    call check_derivatives_numerically(n, npack, nbdirs, uplo, trans, diag, nsize, incx_val, ap_orig, ap_dv_seed, x_orig, x_dv_seed, x_dv, passed)
    deallocate(ap, x, ap_dv, x_dv, ap_orig, x_orig, ap_dv_seed, x_dv_seed)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, npack, nbdirs, uplo, trans, diag, nsize, incx_val, ap_orig, ap_dv_seed, x_orig, x_dv_seed, x_dv, passed)
    integer, intent(in) :: n, npack, nbdirs
    character, intent(in) :: uplo, trans, diag
    integer, intent(in) :: nsize, incx_val
    real(8), intent(in) :: ap_orig(npack), ap_dv_seed(nbdirs,npack), x_orig(n), x_dv_seed(nbdirs,n), x_dv(nbdirs,n)
    logical, intent(out) :: passed
    real(8), parameter :: h = 1.0e-7
    real(8) :: abs_error, abs_ref, err_bound, relative_error, max_error
    real(8), dimension(npack) :: ap_t
    real(8), dimension(n) :: x_t, x_plus, x_minus
    integer :: idir, ii
    logical :: has_err
    has_err = .false.
    max_error = 0.0d0
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    do idir = 1, nbdirs
      ap_t = ap_orig + h * ap_dv_seed(idir,:)
      x_t = x_orig + h * x_dv_seed(idir,:)
      call dtpmv(uplo, trans, diag, nsize, ap_t, x_t, incx_val)
      x_plus = x_t
      ap_t = ap_orig - h * ap_dv_seed(idir,:)
      x_t = x_orig - h * x_dv_seed(idir,:)
      call dtpmv(uplo, trans, diag, nsize, ap_t, x_t, incx_val)
      x_minus = x_t
      do ii = 1, min(2, n)
        abs_error = abs((x_plus(ii) - x_minus(ii)) / (2.0d0 * h) - x_dv(idir,ii))
        abs_ref = abs(x_dv(idir,ii))
        err_bound = 1.0e-5 + 1.0e-5 * abs_ref
        if (abs_error > err_bound) then
          has_err = .true.
          relative_error = abs_error / max(abs_ref, 1.0e-10)
          write(*,*) 'Large error direction', idir, ' X(', ii, '): abs_err=', abs_error, ' rel_err=', relative_error
        end if
        relative_error = abs_error / max(abs_ref, 1.0e-10)
        max_error = max(max_error, relative_error)
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
end program test_dtpmv_vector_forward