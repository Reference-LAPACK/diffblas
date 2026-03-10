! Test program for DTRSM differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_dtrsm
  implicit none

  external :: dtrsm
  external :: dtrsm_d

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing DTRSM (multi-size: n = 4)'
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

    character :: side
    character :: uplo
    character :: transa
    character :: diag
    integer :: msize
    integer :: nsize
    real(8) :: alpha
    real(8), dimension(n,n) :: a
    integer :: lda_val
    real(8), dimension(n,n) :: b
    integer :: ldb_val

    ! Derivative variables
    real(8), dimension(n,n) :: a_d
    real(8), dimension(n,n) :: b_d
    real(8) :: alpha_d

    ! Array restoration and derivative storage
    real(8), dimension(n,n) :: a_orig, a_d_orig
    real(8), dimension(n,n) :: b_orig, b_d_orig
    real(8) :: alpha_orig, alpha_d_orig
    integer :: i, j

    side = 'L'
    uplo = 'U'
    transa = 'N'
    diag = 'N'
    msize = n
    nsize = n
    lda_val = n
    ldb_val = n

    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0  ! Scale to [-1,1]
    call random_number(a)
    a = a * 2.0d0 - 1.0d0  ! Scale to [-1,1]
    call random_number(b)
    b = b * 2.0d0 - 1.0d0  ! Scale to [-1,1]

    ! Initialize input derivatives
    call random_number(a_d)
    a_d = a_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
    call random_number(b_d)
    b_d = b_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
    call random_number(alpha_d)
    alpha_d = alpha_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]

    ! Store _orig and _d_orig
    a_d_orig = a_d
    b_d_orig = b_d
    alpha_d_orig = alpha_d
    a_orig = a
    b_orig = b
    alpha_orig = alpha

    write(*,*) 'Testing DTRSM (n =', n, ')'
    b_orig = b

    ! Call the differentiated function
    call dtrsm_d(side, uplo, transa, diag, msize, nsize, alpha, alpha_d, a, a_d, lda_val, b, b_d, ldb_val)

    write(*,*) 'Function calls completed successfully'

    ! Numerical differentiation check
    call check_derivatives_numerically(n, transa, uplo, side, diag, msize, nsize, lda_val, ldb_val, a_orig, alpha_orig, b_orig, a_d_orig, alpha_d_orig, b_d_orig, b_d, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, transa, uplo, side, diag, msize, nsize, lda_val, ldb_val, a_orig, alpha_orig, b_orig, a_d_orig, alpha_d_orig, b_d_orig, b_d, passed)
    implicit none
    integer, intent(in) :: n
    character, intent(in) :: transa
    character, intent(in) :: uplo
    character, intent(in) :: side
    character, intent(in) :: diag
    integer, intent(in) :: msize
    integer, intent(in) :: nsize
    integer, intent(in) :: lda_val
    integer, intent(in) :: ldb_val
    real(8), intent(in) :: a_orig(n,n), a_d_orig(n,n)
    real(8), intent(in) :: alpha_orig, alpha_d_orig
    real(8), intent(in) :: b_orig(n,n), b_d_orig(n,n)
    real(8), intent(in) :: b_d(n,n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-6  ! Step size for finite differences
    real(8) :: relative_error, max_error
    real(8) :: abs_error, abs_reference, error_bound
    real(8) :: central_diff, ad_result
    logical :: has_large_errors
    real(8), dimension(n,n) :: b_forward, b_backward
    integer :: i, j
    real(8), dimension(n,n) :: a
    real(8) :: alpha
    real(8), dimension(n,n) :: b

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    ! Forward perturbation: f(x + h)
    a = a_orig + h * a_d_orig
    alpha = alpha_orig + h * alpha_d_orig
    b = b_orig + h * b_d_orig
    call dtrsm(side, uplo, transa, diag, msize, nsize, alpha, a, lda_val, b, ldb_val)
    b_forward = b

    ! Backward perturbation: f(x - h)
    a = a_orig - h * a_d_orig
    alpha = alpha_orig - h * alpha_d_orig
    b = b_orig - h * b_d_orig
    call dtrsm(side, uplo, transa, diag, msize, nsize, alpha, a, lda_val, b, ldb_val)
    b_backward = b

    ! Compute central differences and compare with AD results
    do j = 1, min(2, n)
      do i = 1, min(2, n)
        central_diff = (b_forward(i,j) - b_backward(i,j)) / (2.0e0 * h)
        ad_result = b_d(i,j)
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 1.0e-5 + 1.0e-5 * abs_reference
        if (abs_error > error_bound) then
          has_large_errors = .true.
          relative_error = abs_error / max(abs_reference, 1.0e-10)
          write(*,*) 'Large error in output B(', i, ',', j, '):'
          write(*,*) '  Central diff: ', central_diff
          write(*,*) '  AD result:   ', ad_result
          write(*,*) '  Absolute error:', abs_error
          write(*,*) '  Error bound:', error_bound
          write(*,*) '  Relative error:', relative_error
        end if
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        max_error = max(max_error, relative_error)
      end do
    end do

    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_dtrsm