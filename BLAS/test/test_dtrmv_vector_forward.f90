! Test program for DTRMV vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirs=n
! Multi-size test with outlined run_test_for_size(n) - TRMV/TRSV

program test_dtrmv_vector_forward
  implicit none

  external :: dtrmv
  external :: dtrmv_dv

  integer :: nbdirs
  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing DTRMV (Vector Forward, multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    nbdirs = test_sizes(i)
    call run_test_for_size(n_test, passed, nbdirs)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) then
    write(*,*) 'PASS: All sizes completed successfully'
  else
    write(*,*) 'FAIL: One or more sizes had derivative errors'
  end if

contains

  subroutine run_test_for_size(n, passed, nbdirs)
    implicit none
    integer, intent(in) :: n
    logical, intent(out) :: passed
    integer, intent(in) :: nbdirs

    character :: uplo, trans, diag
    integer :: nsize, lda_val, incx_val
    real(8), dimension(n,n) :: a
    real(8), dimension(n) :: x
    real(8), dimension(nbdirs,n,n) :: a_dv
    real(8), dimension(nbdirs,n) :: x_dv
    real(8), dimension(n,n) :: a_orig
    real(8), dimension(nbdirs,n,n) :: a_dv_orig
    real(8), dimension(n) :: x_orig
    real(8), dimension(nbdirs,n) :: x_dv_orig
    integer :: idir, ii, jj
    real(4) :: temp_real, temp_imag

    uplo = 'L'
    trans = 'N'
    diag = 'N'
    nsize = n
    lda_val = n
    incx_val = 1

    ! Lower triangular A (non-unit)
    do jj = 1, n
      do ii = jj, n
        call random_number(temp_real)
        a(ii,jj) = temp_real * 2.0d0 - 1.0d0
      end do
    end do
    do jj = 1, n
      do ii = 1, jj - 1
        a(ii,jj) = 0.0d0
      end do
    end do
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    do idir = 1, nbdirs
      do jj = 1, n
        do ii = jj, n
          call random_number(temp_real)
          a_dv(idir,ii,jj) = temp_real * 2.0d0 - 1.0d0
        end do
      end do
      do jj = 1, n
        do ii = 1, jj - 1
          a_dv(idir,ii,jj) = 0.0d0
        end do
      end do
      call random_number(x_dv(idir,:))
      x_dv(idir,:) = x_dv(idir,:) * 2.0d0 - 1.0d0
    end do

    a_orig = a
    a_dv_orig = a_dv
    x_orig = x
    x_dv_orig = x_dv

    write(*,*) 'Testing DTRMV (Vector Forward, n =', n, ')'

    call dtrmv_dv(uplo, trans, diag, nsize, a, a_dv, lda_val, x, x_dv, incx_val, nbdirs)

    call check_derivatives_numerically(n, nbdirs, uplo, trans, diag, nsize, lda_val, incx_val, a_orig, a_dv_orig, x_orig, x_dv_orig, x_dv, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, nbdirs, uplo, trans, diag, nsize, lda_val, incx_val, a_orig, a_dv_orig, x_orig, x_dv_orig, x_dv, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    character, intent(in) :: uplo, trans, diag
    integer, intent(in) :: nsize, lda_val, incx_val
    real(8), intent(in) :: a_orig(n,n), a_dv_orig(nbdirs,n,n)
    real(8), intent(in) :: x_orig(n), x_dv_orig(nbdirs,n)
    real(8), intent(in) :: x_dv(nbdirs,n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: relative_error, max_error, abs_error, abs_reference, error_bound
    real(8) :: central_diff, ad_result
    real(8), dimension(n) :: x_forward, x_backward
    real(8), dimension(n,n) :: a
    real(8), dimension(n) :: x
    integer :: i, idir
    logical :: has_large_errors

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    do idir = 1, nbdirs
      a = a_orig + h * a_dv_orig(idir,:,:)
      x = x_orig + h * x_dv_orig(idir,:)
      call dtrmv(uplo, trans, diag, nsize, a, lda_val, x, incx_val)
      x_forward = x
      a = a_orig - h * a_dv_orig(idir,:,:)
      x = x_orig - h * x_dv_orig(idir,:)
      call dtrmv(uplo, trans, diag, nsize, a, lda_val, x, incx_val)
      x_backward = x
      do i = 1, min(4, n)
        central_diff = (x_forward(i) - x_backward(i)) / (2.0e0 * h)
        ad_result = x_dv(idir,i)
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 1.0e-5 + 1.0e-5 * abs_reference
        if (abs_error > error_bound) has_large_errors = .true.
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        max_error = max(max_error, relative_error)
      end do
    end do

    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_dtrmv_vector_forward