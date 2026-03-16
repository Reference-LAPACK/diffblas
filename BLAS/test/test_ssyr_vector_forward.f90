! Test program for SSYR vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=n
! Multi-size outlined run_test_for_size(n) - SYR/SYR2

program test_ssyr_vector_forward
  implicit none

  external :: ssyr
  external :: ssyr_dv

  integer :: nbdirs, n_test, seed_array(33), test_sizes(3), i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing SSYR (Vector Forward, multi-size: n = 4)'
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
    implicit none
    integer, intent(in) :: n, nbdirs
    logical, intent(out) :: passed
    character :: uplo
    integer :: nsize, lda_val, incx_val, incy_val
    real(4) :: alpha
    real(4), dimension(n) :: x
    real(4), dimension(n,n) :: a
    real(4), dimension(nbdirs) :: alpha_dv
    real(4), dimension(nbdirs,n) :: x_dv
    real(4), dimension(nbdirs,n,n) :: a_dv
    real(4) :: alpha_orig
    real(4), dimension(nbdirs) :: alpha_dv_seed
    real(4), dimension(n) :: x_orig
    real(4), dimension(nbdirs,n) :: x_dv_seed
    real(4), dimension(n,n) :: a_orig
    real(4), dimension(nbdirs,n,n) :: a_dv_seed
    integer :: ii, jj, idir
    real(4) :: temp_real, temp_imag

    uplo = 'U'
    nsize = n
    lda_val = n
    incx_val = 1
    incy_val = 1

    call random_number(temp_real)
    alpha = temp_real * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    call random_number(a)
    a = a * 2.0d0 - 1.0d0
    do jj = 1, n
      do ii = jj + 1, n
        a(ii,jj) = a(jj,ii)
      end do
    end do
    do idir = 1, nbdirs
      call random_number(temp_real)
      alpha_dv(idir) = temp_real * 2.0d0 - 1.0d0
    end do
    do idir = 1, nbdirs
      call random_number(x_dv(idir,:))
      x_dv(idir,:) = x_dv(idir,:) * 2.0d0 - 1.0d0
    end do
    do idir = 1, nbdirs
      call random_number(a_dv(idir,:,:))
      a_dv(idir,:,:) = a_dv(idir,:,:) * 2.0d0 - 1.0d0
      do jj = 1, n
        do ii = jj + 1, n
          a_dv(idir,ii,jj) = a_dv(idir,jj,ii)
        end do
      end do
    end do

    write(*,*) 'Testing SSYR (Vector Forward, n =', n, ')'
    alpha_orig = alpha
    alpha_dv_seed = alpha_dv
    x_orig = x
    x_dv_seed = x_dv
    a_orig = a
    a_dv_seed = a_dv

    call ssyr_dv(uplo, nsize, alpha, alpha_dv, x, x_dv, incx_val, a, a_dv, lda_val, nbdirs)

    call check_derivatives_numerically(n, nbdirs, uplo, nsize, lda_val, incx_val, alpha_orig, alpha_dv_seed, x_orig, x_dv_seed, a_orig, a_dv_seed, a_dv, passed)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, nbdirs, uplo, nsize, lda_val, incx_val, alpha_orig, alpha_dv_seed, x_orig, x_dv_seed, a_orig, a_dv_seed, a_dv, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    character, intent(in) :: uplo
    integer, intent(in) :: nsize, lda_val, incx_val
    real(4), intent(in) :: alpha_orig
    real(4), intent(in) :: alpha_dv_seed(nbdirs), x_orig(n), x_dv_seed(nbdirs,n)
    real(4), intent(in) :: a_orig(n,n), a_dv_seed(nbdirs,n,n), a_dv(nbdirs,n,n)
    logical, intent(out) :: passed
    real(4), parameter :: h = 1.0e-3
    real(4), dimension(n,n) :: a_fwd, a_bwd
    real(4) :: alpha_t
    real(4), dimension(n) :: x_t
    real(4), dimension(n,n) :: a_t
    integer :: idir, i, j
    logical :: has_err
    real(4) :: abs_error, abs_ref, err_bound, max_error, relative_error
    has_err = .false.
    max_error = 0.0d0
    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    do idir = 1, nbdirs
      alpha_t = alpha_orig + h * alpha_dv_seed(idir)
      x_t = x_orig + h * x_dv_seed(idir,:)
      a_t = a_orig + h * a_dv_seed(idir,:,:)
      call ssyr(uplo, nsize, alpha_t, x_t, incx_val, a_t, lda_val)
      a_fwd = a_t
      alpha_t = alpha_orig - h * alpha_dv_seed(idir)
      x_t = x_orig - h * x_dv_seed(idir,:)
      a_t = a_orig - h * a_dv_seed(idir,:,:)
      call ssyr(uplo, nsize, alpha_t, x_t, incx_val, a_t, lda_val)
      a_bwd = a_t
      do j = 1, min(2, n)
        do i = 1, min(2, n)
          abs_error = abs((a_fwd(i,j) - a_bwd(i,j)) / (2.0e0 * h) - a_dv(idir,i,j))
          abs_ref = abs(a_dv(idir,i,j))
          err_bound = 2.0e-3 + 2.0e-3 * abs_ref
          if (abs_error > err_bound) has_err = .true.
          relative_error = 0.0d0
          if (abs_ref > 1.0d-10) relative_error = abs_error / abs_ref
          if (relative_error > max_error) max_error = relative_error
        end do
      end do
    end do
    passed = .not. has_err
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
    if (has_err) write(*,*) 'FAIL: Derivatives are outside tolerance'
    if (.not. has_err) write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
  end subroutine check_derivatives_numerically

end program test_ssyr_vector_forward