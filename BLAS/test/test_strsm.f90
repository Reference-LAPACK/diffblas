! Test program for STRSM differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision

program test_strsm
  implicit none

  external :: strsm
  external :: strsm_d

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension (rows/cols of matrices)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions

  character :: side
  character :: uplo
  character :: transa
  character :: diag
  integer :: msize
  integer :: nsize
  real(4) :: alpha
  real(4), dimension(max_size,max_size) :: a
  integer :: lda_val
  real(4), dimension(max_size,max_size) :: b
  integer :: ldb_val

  ! Derivative variables
  real(4) :: alpha_d
  real(4), dimension(max_size,max_size) :: a_d
  real(4), dimension(max_size,max_size) :: b_d

  ! Storage variables for inout parameters
  real(4), dimension(max_size,max_size) :: b_output

  ! Array restoration variables for numerical differentiation
  real(4) :: alpha_orig
  real(4), dimension(max_size,max_size) :: a_orig
  real(4), dimension(max_size,max_size) :: b_orig

  ! Variables for central difference computation
  real(4), dimension(max_size,max_size) :: b_forward, b_backward
  ! Scalar variables for central difference computation
  real(4) :: central_diff, ad_result
  logical :: has_large_errors

  ! Variables for storing original derivative values
  real(4) :: alpha_d_orig
  real(4), dimension(max_size,max_size) :: a_d_orig
  real(4), dimension(max_size,max_size) :: b_d_orig

  ! Temporary variables for matrix initialization
  real(4) :: temp_real, temp_imag
  integer :: i, j

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  integer :: seed_array(33)
  seed_array = 42
  call random_seed(put=seed_array)

  side = 'L'
  uplo = 'U'
  transa = 'N'
  diag = 'N'
  msize = n
  nsize = n
  call random_number(alpha)
  alpha = alpha * 2.0 - 1.0  ! Scale to [-1,1]
  call random_number(a)
  a = a * 2.0d0 - 1.0d0  ! Scale to [-1,1]
  lda_val = lda  ! LDA must be at least max( 1
  call random_number(b)
  b = b * 2.0d0 - 1.0d0  ! Scale to [-1,1]
  ldb_val = ldb

  ! Initialize input derivatives to random values
  call random_number(alpha_d)
  alpha_d = alpha_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
  call random_number(a_d)
  a_d = a_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
  call random_number(b_d)
  b_d = b_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]

  ! Store initial derivative values after random initialization
  alpha_d_orig = alpha_d
  a_d_orig = a_d
  b_d_orig = b_d

  ! Store original values for central difference computation
  alpha_orig = alpha
  a_orig = a
  b_orig = b

  write(*,*) 'Testing STRSM'
  ! Store input values of inout parameters before first function call
  b_orig = b

  ! Re-initialize data for differentiated function
  ! Only reinitialize inout parameters - keep input-only parameters unchanged

  ! side already has correct value from original call
  ! uplo already has correct value from original call
  ! transa already has correct value from original call
  ! diag already has correct value from original call
  msize = n
  nsize = n
  ! alpha already has correct value from original call
  ! a already has correct value from original call
  lda_val = lda  ! LDA must be at least max( 1
  b = b_orig
  ldb_val = ldb

  ! Call the differentiated function
  call strsm_d(side, uplo, transa, diag, msize, nsize, alpha, alpha_d, a, a_d, lda_val, b, b_d, ldb_val)

  ! Print results and compare
  write(*,*) 'Function calls completed successfully'

  ! Numerical differentiation check
  call check_derivatives_numerically()

  write(*,*) 'Test completed successfully'

contains

  subroutine check_derivatives_numerically()
    implicit none
    real(4), parameter :: h = 1.0e-3  ! Step size for finite differences
    real(4) :: relative_error, max_error
    real(4) :: output_orig, output_pert
    real(4) :: numerical_result, analytical_result
    real(4) :: abs_error, abs_reference, error_bound
    integer :: i, j
    
    max_error = 0.0e0
    has_large_errors = .false.
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3
    
    ! Original values already stored in main program
    
    ! Central difference computation: f(x + h) - f(x - h) / (2h)
    ! Forward perturbation: f(x + h)
    alpha = alpha_orig + h * alpha_d_orig
    a = a_orig + h * a_d_orig
    b = b_orig + h * b_d_orig
    call strsm(side, uplo, transa, diag, msize, nsize, alpha, a, lda_val, b, ldb_val)
    ! Store forward perturbation results
    b_forward = b
    
    ! Backward perturbation: f(x - h)
    alpha = alpha_orig - h * alpha_d_orig
    a = a_orig - h * a_d_orig
    b = b_orig - h * b_d_orig
    call strsm(side, uplo, transa, diag, msize, nsize, alpha, a, lda_val, b, ldb_val)
    ! Store backward perturbation results
    b_backward = b
    
    ! Compute central differences and compare with AD results
    ! Check derivatives for output B
    do j = 1, min(2, n)  ! Check only first few elements
      do i = 1, min(2, n)
        ! Central difference: (f(x+h) - f(x-h)) / (2h)
        central_diff = (b_forward(i,j) - b_backward(i,j)) / (2.0e0 * h)
        ! AD result
        ad_result = b_d(i,j)
        ! Error check: |a - b| > atol + rtol * |b|
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 2.0e-3 + 2.0e-3 * abs_reference
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
        ! Track max error for reporting (normalized)
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        max_error = max(max_error, relative_error)
      end do
    end do
    
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_strsm