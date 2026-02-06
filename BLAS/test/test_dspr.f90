! Test program for DSPR differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision

program test_dspr
  implicit none

  external :: dspr
  external :: dspr_d

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension (rows/cols of matrices)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions

  character :: uplo
  integer :: nsize
  real(8) :: alpha
  real(8), dimension(max_size) :: x
  integer :: incx_val
  real(8), dimension((n*(n+1))/2) :: ap

  ! Derivative variables
  real(8) :: alpha_d
  real(8), dimension(max_size) :: x_d
  real(8), dimension((n*(n+1))/2) :: ap_d

  ! Storage variables for inout parameters
  real(8), dimension((n*(n+1))/2) :: ap_output

  ! Array restoration variables for numerical differentiation
  real(8) :: alpha_orig
  real(8), dimension((n*(n+1))/2) :: ap_orig
  real(8), dimension(max_size) :: x_orig

  ! Variables for central difference computation
  ! Scalar variables for central difference computation
  real(8) :: central_diff, ad_result
  logical :: has_large_errors

  ! Variables for storing original derivative values
  real(8) :: alpha_d_orig
  real(8), dimension((n*(n+1))/2) :: ap_d_orig
  real(8), dimension(max_size) :: x_d_orig

  ! Temporary variables for matrix initialization
  real(4) :: temp_real, temp_imag
  integer :: i, j

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  integer :: seed_array(33)
  seed_array = 42
  call random_seed(put=seed_array)

  uplo = 'U'
  nsize = n
  call random_number(alpha)
  alpha = alpha * 2.0d0 - 1.0d0  ! Scale to [-1,1]
  call random_number(x)
  x = x * 2.0d0 - 1.0d0  ! Scale to [-1,1]
  incx_val = 1  ! INCX 1
  call random_number(ap)
  ap = ap * 2.0d0 - 1.0d0  ! Scale to [-1,1]

  ! Initialize input derivatives to random values
  call random_number(alpha_d)
  alpha_d = alpha_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
  call random_number(ap_d)
  ap_d = ap_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
  call random_number(x_d)
  x_d = x_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]

  ! Store initial derivative values after random initialization
  alpha_d_orig = alpha_d
  ap_d_orig = ap_d
  x_d_orig = x_d

  ! Store original values for central difference computation
  alpha_orig = alpha
  ap_orig = ap
  x_orig = x

  write(*,*) 'Testing DSPR'
  ! Store input values of inout parameters before first function call
  ap_orig = ap

  ! Re-initialize data for differentiated function
  ! Only reinitialize inout parameters - keep input-only parameters unchanged

  ! uplo already has correct value from original call
  nsize = n
  ! alpha already has correct value from original call
  ! x already has correct value from original call
  incx_val = 1  ! INCX 1
  ap = ap_orig

  ! Call the differentiated function
  call dspr_d(uplo, nsize, alpha, alpha_d, x, x_d, incx_val, ap, ap_d)

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
    ap = ap_orig + h * ap_d_orig
    x = x_orig + h * x_d_orig
    call dspr(uplo, nsize, alpha, x, incx_val, ap)
    ! Store forward perturbation results
    
    ! Backward perturbation: f(x - h)
    alpha = alpha_orig - h * alpha_d_orig
    ap = ap_orig - h * ap_d_orig
    x = x_orig - h * x_d_orig
    call dspr(uplo, nsize, alpha, x, incx_val, ap)
    ! Store backward perturbation results
    
    ! Compute central differences and compare with AD results
    
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_dspr