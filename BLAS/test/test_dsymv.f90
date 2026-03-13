! Test program for DSYMV differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_dsymv
  implicit none

  external :: dsymv
  external :: dsymv_d

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing DSYMV (multi-size: n = 4)'
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

    character :: uplo
    integer :: nsize
    real(8) :: alpha
    real(8), dimension(n,n) :: a
    integer :: lda_val
    real(8), dimension(n) :: x
    integer :: incx
    real(8) :: beta
    real(8), dimension(n) :: y
    integer :: incy

    ! Derivative variables
    real(8) :: beta_d
    real(8) :: alpha_d
    real(8), dimension(n,n) :: a_d
    real(8), dimension(n) :: x_d
    real(8), dimension(n) :: y_d

    ! Array restoration and derivative storage
    real(8) :: beta_orig, beta_d_orig
    real(8) :: alpha_orig, alpha_d_orig
    real(8), dimension(n,n) :: a_orig, a_d_orig
    real(8), dimension(n) :: x_orig, x_d_orig
    real(8), dimension(n) :: y_orig, y_d_orig
    integer :: i, j

    uplo = 'U'
    nsize = n
    lda_val = n
    incx = 1
    incy = 1

    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0  ! Scale to [-1,1]
    call random_number(a)
    a = a * 2.0d0 - 1.0d0  ! Scale to [-1,1]
    call random_number(x)
    x = x * 2.0d0 - 1.0d0  ! Scale to [-1,1]
    call random_number(beta)
    beta = beta * 2.0d0 - 1.0d0  ! Scale to [-1,1]
    call random_number(y)
    y = y * 2.0d0 - 1.0d0  ! Scale to [-1,1]

    ! Initialize input derivatives
    call random_number(beta_d)
    beta_d = beta_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
    call random_number(alpha_d)
    alpha_d = alpha_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
    call random_number(a_d)
    a_d = a_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
    call random_number(x_d)
    x_d = x_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
    call random_number(y_d)
    y_d = y_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]

    ! Store _orig and _d_orig
    beta_d_orig = beta_d
    alpha_d_orig = alpha_d
    a_d_orig = a_d
    x_d_orig = x_d
    y_d_orig = y_d
    beta_orig = beta
    alpha_orig = alpha
    a_orig = a
    x_orig = x
    y_orig = y

    write(*,*) 'Testing DSYMV (n =', n, ')'
    y_orig = y

    ! Call the differentiated function
    call dsymv_d(uplo, nsize, alpha, alpha_d, a, a_d, lda_val, x, x_d, 1, beta, beta_d, y, y_d, 1)

    write(*,*) 'Function calls completed successfully'

    ! Numerical differentiation check
    call check_derivatives_numerically(n, uplo, nsize, lda_val, beta_orig, alpha_orig, a_orig, x_orig, y_orig, beta_d_orig, alpha_d_orig, a_d_orig, x_d_orig, y_d_orig, y_d, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, uplo, nsize, lda_val, beta_orig, alpha_orig, a_orig, x_orig, y_orig, beta_d_orig, alpha_d_orig, a_d_orig, x_d_orig, y_d_orig, y_d, passed)
    implicit none
    integer, intent(in) :: n
    character, intent(in) :: uplo
    integer, intent(in) :: nsize
    integer, intent(in) :: lda_val
    real(8), intent(in) :: beta_orig, beta_d_orig
    real(8), intent(in) :: alpha_orig, alpha_d_orig
    real(8), intent(in) :: a_orig(n,n), a_d_orig(n,n)
    real(8), intent(in) :: x_orig(n), x_d_orig(n)
    real(8), intent(in) :: y_orig(n), y_d_orig(n)
    real(8), intent(in) :: y_d(n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-6  ! Step size for finite differences
    real(8) :: relative_error, max_error
    real(8) :: abs_error, abs_reference, error_bound
    real(8) :: central_diff, ad_result
    logical :: has_large_errors
    real(8), dimension(n) :: y_forward, y_backward
    integer :: i, j
    real(8) :: beta
    real(8) :: alpha
    real(8), dimension(n,n) :: a
    real(8), dimension(n) :: x
    real(8), dimension(n) :: y

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    ! Forward perturbation: f(x + h)
    beta = beta_orig + h * beta_d_orig
    alpha = alpha_orig + h * alpha_d_orig
    a = a_orig + h * a_d_orig
    x = x_orig + h * x_d_orig
    y = y_orig + h * y_d_orig
    call dsymv(uplo, nsize, alpha, a, lda_val, x, 1, beta, y, 1)
    y_forward = y

    ! Backward perturbation: f(x - h)
    beta = beta_orig - h * beta_d_orig
    alpha = alpha_orig - h * alpha_d_orig
    a = a_orig - h * a_d_orig
    x = x_orig - h * x_d_orig
    y = y_orig - h * y_d_orig
    call dsymv(uplo, nsize, alpha, a, lda_val, x, 1, beta, y, 1)
    y_backward = y

    ! Compute central differences and compare with AD results
    do i = 1, n
        central_diff = (y_forward(i) - y_backward(i)) / (2.0e0 * h)
        ad_result = y_d(i)
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 1.0e-5 + 1.0e-5 * abs_reference
        if (abs_error > error_bound) then
          has_large_errors = .true.
          relative_error = abs_error / max(abs_reference, 1.0e-10)
          write(*,*) 'Large error in output Y(', i, '):'
          write(*,*) '  Central diff: ', central_diff
          write(*,*) '  AD result:   ', ad_result
          write(*,*) '  Absolute error:', abs_error
          write(*,*) '  Error bound:', error_bound
          write(*,*) '  Relative error:', relative_error
        end if
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        max_error = max(max_error, relative_error)
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

end program test_dsymv