! Test program for SSPMV vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Multi-size outlined - SPMV vector forward

program test_sspmv_vector_forward
  implicit none
  external :: sspmv
  external :: sspmv_dv
  integer :: nbdirs, n_test, seed_array(33), test_sizes(3), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing SSPMV (Vector Forward, multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 3
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
    real(4) :: alpha, beta
    real(4), dimension(n) :: x, y, y_orig, y_plus, y_minus
    real(4), dimension(nbdirs) :: alpha_dv, beta_dv
    real(4), dimension(nbdirs,n) :: x_dv, y_dv, y_dv_seed
    real(4), dimension(:), allocatable :: ap
    real(4), dimension(:,:), allocatable :: ap_dv
    real(4), dimension(:), allocatable :: ap_orig, ap_t
    real(4), parameter :: h = 1.0e-3
    real(4) :: max_err, abs_ref
    integer :: ii
    write(*,*) 'Testing SSPMV (Vector Forward, n =', n, ')'
    uplo = 'U'
    nsize = n
    incx_val = 1
    incy_val = 1
    npack = (n * (n + 1)) / 2
    allocate(ap(npack), ap_dv(nbdirs, npack), ap_orig(npack), ap_t(npack))
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
    do k = 1, nbdirs
      call random_number(alpha_dv(k))
      alpha_dv(k) = alpha_dv(k) * 2.0d0 - 1.0d0
      call random_number(beta_dv(k))
      beta_dv(k) = beta_dv(k) * 2.0d0 - 1.0d0
      call random_number(x_dv(k,:))
      x_dv(k,:) = x_dv(k,:) * 2.0d0 - 1.0d0
      call random_number(y_dv(k,:))
      y_dv(k,:) = y_dv(k,:) * 2.0d0 - 1.0d0
      call random_number(ap_dv(k,:))
      ap_dv(k,:) = ap_dv(k,:) * 2.0d0 - 1.0d0
    end do
    ap_orig = ap
    y_orig = y
    y_dv_seed = y_dv
    call sspmv_dv(uplo, nsize, alpha, alpha_dv, ap, ap_dv, x, x_dv, incx_val, beta, beta_dv, y, y_dv, incy_val, nbdirs)
    max_err = 0.0d0
    do k = 1, nbdirs
      y_plus = y_orig + h * y_dv_seed(k,:)
      y_minus = y_orig - h * y_dv_seed(k,:)
      ap_t = ap_orig + h * ap_dv(k,:)
      call sspmv(uplo, nsize, alpha + h*alpha_dv(k), ap_t, x + h*x_dv(k,:), incx_val, beta + h*beta_dv(k), y_plus, incy_val)
      ap_t = ap_orig - h * ap_dv(k,:)
      call sspmv(uplo, nsize, alpha - h*alpha_dv(k), ap_t, x - h*x_dv(k,:), incx_val, beta - h*beta_dv(k), y_minus, incy_val)
      do ii = 1, n
        max_err = max(max_err, abs((y_plus(ii) - y_minus(ii)) / (2.0d0 * h) - y_dv(k,ii)))
      end do
    end do
    abs_ref = maxval(abs(y_dv)) + 1.0d0
    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Maximum relative error:', max_err / abs_ref
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
    passed = (max_err <= 2.0e-3 * abs_ref)
    if (.not. passed) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
    deallocate(ap, ap_dv, ap_orig, ap_t)
  end subroutine run_test_for_size
end program test_sspmv_vector_forward