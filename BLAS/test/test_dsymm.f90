! Test program for DSYMM differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision

program test_dsymm
  implicit none

  external :: dsymm
  external :: dsymm_d

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension (rows/cols of matrices)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions

  character :: side
  character :: uplo
  integer :: msize
  integer :: nsize
  real(8) :: alpha
  real(8), dimension(max_size,max_size) :: a
  integer :: lda_val
  real(8), dimension(max_size,max_size) :: b
  integer :: ldb_val
  real(8) :: beta
  real(8), dimension(max_size,max_size) :: c
  integer :: ldc_val

  ! Derivative variables
  real(8) :: alpha_d
  real(8), dimension(max_size,max_size) :: a_d
  real(8), dimension(max_size,max_size) :: b_d
  real(8) :: beta_d
  real(8), dimension(max_size,max_size) :: c_d

  ! Storage variables for inout parameters
  real(8), dimension(max_size,max_size) :: c_output

  ! Array restoration variables for numerical differentiation
  real(8) :: alpha_orig
  real(8) :: beta_orig
  real(8), dimension(max_size,max_size) :: a_orig
  real(8), dimension(max_size,max_size) :: c_orig
  real(8), dimension(max_size,max_size) :: b_orig

  ! Variables for central difference computation
  real(8), dimension(max_size,max_size) :: c_forward, c_backward
  ! Scalar variables for central difference computation
  real(8) :: central_diff, ad_result
  logical :: has_large_errors

  ! Variables for storing original derivative values
  real(8) :: alpha_d_orig
  real(8) :: beta_d_orig
  real(8), dimension(max_size,max_size) :: a_d_orig
  real(8), dimension(max_size,max_size) :: c_d_orig
  real(8), dimension(max_size,max_size) :: b_d_orig

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
  msize = n
  nsize = n
  call random_number(alpha)
  alpha = alpha * 2.0d0 - 1.0d0  ! Scale to [-1,1]
  ! Initialize a as symmetric matrix
  ! Fill upper triangle with random numbers
  do i = 1, lda
    do j = i, lda
      call random_number(temp_real)
      a(i,j) = temp_real * 2.0 - 1.0  ! Scale to [-1,1]
    end do
  end do
  
  ! Fill lower triangle with symmetric values
  do i = 2, lda
    do j = 1, i-1
      a(i,j) = a(j,i)  ! A(i,j) = A(j,i)
    end do
  end do
  lda_val = lda  ! LDA must be at least max( 1, m )
  call random_number(b)
  b = b * 2.0d0 - 1.0d0  ! Scale to [-1,1]
  ldb_val = ldb
  call random_number(beta)
  beta = beta * 2.0d0 - 1.0d0  ! Scale to [-1,1]
  call random_number(c)
  c = c * 2.0d0 - 1.0d0  ! Scale to [-1,1]
  ldc_val = ldc

  ! Initialize input derivatives to random values
  call random_number(alpha_d)
  alpha_d = alpha_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
  call random_number(beta_d)
  beta_d = beta_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
  ! Initialize a_d as symmetric matrix
  ! Fill upper triangle with random numbers
  do i = 1, lda
    do j = i, lda
      call random_number(temp_real)
      a_d(i,j) = temp_real * 2.0 - 1.0  ! Scale to [-1,1]
    end do
  end do
  
  ! Fill lower triangle with symmetric values
  do i = 2, lda
    do j = 1, i-1
      a_d(i,j) = a_d(j,i)  ! A(i,j) = A(j,i)
    end do
  end do
  call random_number(c_d)
  c_d = c_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
  call random_number(b_d)
  b_d = b_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]

  ! Store initial derivative values after random initialization
  alpha_d_orig = alpha_d
  beta_d_orig = beta_d
  a_d_orig = a_d
  c_d_orig = c_d
  b_d_orig = b_d

  ! Store original values for central difference computation
  alpha_orig = alpha
  beta_orig = beta
  a_orig = a
  c_orig = c
  b_orig = b

  write(*,*) 'Testing DSYMM'
  ! Store input values of inout parameters before first function call
  c_orig = c

  ! Re-initialize data for differentiated function
  ! Only reinitialize inout parameters - keep input-only parameters unchanged

  ! side already has correct value from original call
  ! uplo already has correct value from original call
  msize = n
  nsize = n
  ! alpha already has correct value from original call
  ! a already has correct value from original call
  lda_val = lda  ! LDA must be at least max( 1, m )
  ! b already has correct value from original call
  ldb_val = ldb
  ! beta already has correct value from original call
  c = c_orig
  ldc_val = ldc

  ! Call the differentiated function
  call dsymm_d(side, uplo, msize, nsize, alpha, alpha_d, a, a_d, lda_val, b, b_d, ldb_val, beta, beta_d, c, c_d, ldc_val)

  ! Print results and compare
  write(*,*) 'Function calls completed successfully'

  ! Numerical differentiation check
  call check_derivatives_numerically()

  write(*,*) 'Test completed successfully'

contains

  subroutine check_derivatives_numerically()
    implicit none
    real(8), parameter :: h = 1.0e-6  ! Step size for finite differences
    real(8) :: relative_error, max_error
    real(8) :: output_orig, output_pert
    real(8) :: numerical_result, analytical_result
    real(8) :: abs_error, abs_reference, error_bound
    integer :: i, j
    
    max_error = 0.0e0
    has_large_errors = .false.
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5
    
    ! Original values already stored in main program
    
    ! Central difference computation: f(x + h) - f(x - h) / (2h)
    ! Forward perturbation: f(x + h)
    alpha = alpha_orig + h * alpha_d_orig
    beta = beta_orig + h * beta_d_orig
    a = a_orig + h * a_d_orig
    c = c_orig + h * c_d_orig
    b = b_orig + h * b_d_orig
    call dsymm(side, uplo, msize, nsize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
    ! Store forward perturbation results
    c_forward = c
    
    ! Backward perturbation: f(x - h)
    alpha = alpha_orig - h * alpha_d_orig
    beta = beta_orig - h * beta_d_orig
    a = a_orig - h * a_d_orig
    c = c_orig - h * c_d_orig
    b = b_orig - h * b_d_orig
    call dsymm(side, uplo, msize, nsize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
    ! Store backward perturbation results
    c_backward = c
    
    ! Compute central differences and compare with AD results
    ! Check derivatives for output C
    do j = 1, min(2, n)  ! Check only first few elements
      do i = 1, min(2, n)
        ! Central difference: (f(x+h) - f(x-h)) / (2h)
        central_diff = (c_forward(i,j) - c_backward(i,j)) / (2.0e0 * h)
        ! AD result
        ad_result = c_d(i,j)
        ! Error check: |a - b| > atol + rtol * |b|
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
        ! Track max error for reporting (normalized)
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        max_error = max(max_error, relative_error)
      end do
    end do
    
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_dsymm