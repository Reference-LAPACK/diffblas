! Test program for ZGEMM differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_zgemm
  implicit none

  external :: zgemm
  external :: zgemm_d

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(3)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing ZGEMM (multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 3
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

    character :: transa
    character :: transb
    integer :: msize
    integer :: nsize
    integer :: ksize
    complex(8) :: alpha
    complex(8), dimension(n,n) :: a
    integer :: lda_val
    complex(8), dimension(n,n) :: b
    integer :: ldb_val
    complex(8) :: beta
    complex(8), dimension(n,n) :: c
    integer :: ldc_val

    ! Derivative variables
    complex(8) :: alpha_d
    complex(8), dimension(n,n) :: c_d
    complex(8), dimension(n,n) :: a_d
    complex(8), dimension(n,n) :: b_d
    complex(8) :: beta_d

    ! Array restoration and derivative storage
    complex(8) :: alpha_orig, alpha_d_orig
    complex(8), dimension(n,n) :: c_orig, c_d_orig
    complex(8), dimension(n,n) :: a_orig, a_d_orig
    complex(8), dimension(n,n) :: b_orig, b_d_orig
    complex(8) :: beta_orig, beta_d_orig
    real(8) :: temp_re, temp_im  ! For complex random init
    integer :: i, j

    transa = 'N'
    transb = 'N'
    msize = n
    nsize = n
    ksize = n
    lda_val = n
    ldb_val = n
    ldc_val = n

    call random_number(temp_re)
    call random_number(temp_im)
    alpha = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    call random_number(temp_re)
    call random_number(temp_im)
    a = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    call random_number(temp_re)
    call random_number(temp_im)
    b = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    call random_number(temp_re)
    call random_number(temp_im)
    beta = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    call random_number(temp_re)
    call random_number(temp_im)
    c = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)

    ! Initialize input derivatives
    call random_number(temp_re)
    call random_number(temp_im)
    alpha_d = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    call random_number(temp_re)
    call random_number(temp_im)
    c_d = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    call random_number(temp_re)
    call random_number(temp_im)
    a_d = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    call random_number(temp_re)
    call random_number(temp_im)
    b_d = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    call random_number(temp_re)
    call random_number(temp_im)
    beta_d = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)

    ! Store _orig and _d_orig
    alpha_d_orig = alpha_d
    c_d_orig = c_d
    a_d_orig = a_d
    b_d_orig = b_d
    beta_d_orig = beta_d
    alpha_orig = alpha
    c_orig = c
    a_orig = a
    b_orig = b
    beta_orig = beta

    write(*,*) 'Testing ZGEMM (n =', n, ')'
    c_orig = c

    ! Call the differentiated function
    call zgemm_d(transa, transb, msize, nsize, ksize, alpha, alpha_d, a, a_d, lda_val, b, b_d, ldb_val, beta, beta_d, c, c_d, ldc_val)
    alpha_d = alpha_d_orig
    a_d = a_d_orig
    b_d = b_d_orig
    beta_d = beta_d_orig

    write(*,*) 'Function calls completed successfully'

    ! Numerical differentiation check
    call check_derivatives_numerically(n, transa, transb, msize, nsize, ksize, lda_val, ldb_val, ldc_val, alpha_orig, c_orig, beta_orig, b_orig, a_orig, alpha_d_orig, c_d_orig, beta_d_orig, b_d_orig, a_d_orig, c_d, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, transa, transb, msize, nsize, ksize, lda_val, ldb_val, ldc_val, alpha_orig, c_orig, beta_orig, b_orig, a_orig, alpha_d_orig, c_d_orig, beta_d_orig, b_d_orig, a_d_orig, c_d, passed)
    implicit none
    integer, intent(in) :: n
    character, intent(in) :: transa
    character, intent(in) :: transb
    integer, intent(in) :: msize
    integer, intent(in) :: nsize
    integer, intent(in) :: ksize
    integer, intent(in) :: lda_val
    integer, intent(in) :: ldb_val
    integer, intent(in) :: ldc_val
    complex(8), intent(in) :: alpha_orig, alpha_d_orig
    complex(8), intent(in) :: c_orig(n,n), c_d_orig(n,n)
    complex(8), intent(in) :: beta_orig, beta_d_orig
    complex(8), intent(in) :: b_orig(n,n), b_d_orig(n,n)
    complex(8), intent(in) :: a_orig(n,n), a_d_orig(n,n)
    complex(8), intent(in) :: c_d(n,n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-6  ! Step size for finite differences
    real(8) :: relative_error, max_error
    real(8) :: abs_error, abs_reference, error_bound
    real(8) :: central_diff, ad_result
    logical :: has_large_errors
    complex(8), dimension(n,n) :: c_forward, c_backward
    integer :: i, j
    complex(8) :: alpha
    complex(8), dimension(n,n) :: c
    complex(8) :: beta
    complex(8), dimension(n,n) :: b
    complex(8), dimension(n,n) :: a

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    ! Forward perturbation: f(x + h)
    alpha = alpha_orig + h * alpha_d_orig
    c = c_orig + h * c_d_orig
    beta = beta_orig + h * beta_d_orig
    b = b_orig + h * b_d_orig
    a = a_orig + h * a_d_orig
    call zgemm(transa, transb, msize, nsize, ksize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
    c_forward = c

    ! Backward perturbation: f(x - h)
    alpha = alpha_orig - h * alpha_d_orig
    c = c_orig - h * c_d_orig
    beta = beta_orig - h * beta_d_orig
    b = b_orig - h * b_d_orig
    a = a_orig - h * a_d_orig
    call zgemm(transa, transb, msize, nsize, ksize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
    c_backward = c

    ! Compute central differences and compare with AD results
    do j = 1, min(2, n)
      do i = 1, min(2, n)
        central_diff = (c_forward(i,j) - c_backward(i,j)) / (2.0e0 * h)
        ad_result = c_d(i,j)
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 1.0e-5 + 1.0e-5 * abs_reference
        if (abs_error > error_bound) then
          has_large_errors = .true.
          relative_error = abs_error / max(abs_reference, 1.0e-10)
          write(*,*) 'Large error in output C(', i, ',', j, '):'
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
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_zgemm