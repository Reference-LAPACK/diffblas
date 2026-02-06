! Test program for DASUM differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision

program test_dasum
  implicit none

  real(8), external :: dasum
  real(8), external :: dasum_d

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension (rows/cols of matrices)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions

  integer :: nsize
  real(8), dimension(4) :: dx
  integer :: incx_val

  ! Derivative variables
  real(8), dimension(4) :: dx_d

  ! Storage variables for inout parameters

  ! Array restoration variables for numerical differentiation
  real(8), dimension(4) :: dx_orig
  real(8) :: dasum_orig

  ! Variables for central difference computation
  ! Scalar variables for central difference computation
  real(8) :: central_diff, ad_result
  logical :: has_large_errors
  real(8) :: dasum_result, dasum_d_result
  real(8) :: dasum_forward, dasum_backward

  ! Variables for storing original derivative values
  real(8), dimension(4) :: dx_d_orig

  ! Temporary variables for matrix initialization
  real(4) :: temp_real, temp_imag
  integer :: i, j

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  integer :: seed_array(33)
  seed_array = 42
  call random_seed(put=seed_array)

  nsize = n
  call random_number(dx)
  dx = dx * 2.0d0 - 1.0d0  ! Scale to [-1,1]
  incx_val = 1

  ! Initialize input derivatives to random values
  call random_number(dx_d)
  dx_d = dx_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]

  ! Store initial derivative values after random initialization
  dx_d_orig = dx_d

  ! Store original values for central difference computation
  dx_orig = dx

  write(*,*) 'Testing DASUM'
  ! Store input values of inout parameters before first function call

  ! Call the original function
  dasum_result = dasum(nsize, dx, incx_val)

  ! Store output values of inout parameters after first function call

  ! Re-initialize data for differentiated function
  ! Only reinitialize inout parameters - keep input-only parameters unchanged

  nsize = n
  ! dx already has correct value from original call
  incx_val = 1

  ! Call the differentiated function
  dasum_d_result = dasum_d(nsize, dx, dx_d, incx_val, dasum_result)

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
    dx = dx_orig + h * dx_d_orig
    dasum_forward = dasum(nsize, dx, incx_val)
    ! Store forward perturbation results
    ! dasum_forward already captured above
    
    ! Backward perturbation: f(x - h)
    dx = dx_orig - h * dx_d_orig
    dasum_backward = dasum(nsize, dx, incx_val)
    ! Store backward perturbation results
    ! dasum_backward already captured above
    
    ! Compute central differences and compare with AD results
    ! Check derivatives for function DASUM
    ! Central difference: (f(x+h) - f(x-h)) / (2h)
    central_diff = (dasum_forward - dasum_backward) / (2.0e0 * h)
    ! AD result
    ad_result = dasum_d_result
    ! Error check: |a - b| > atol + rtol * |b|
    abs_error = abs(central_diff - ad_result)
    abs_reference = abs(ad_result)
    error_bound = 1.0e-5 + 1.0e-5 * abs_reference
    if (abs_error > error_bound) then
      has_large_errors = .true.
      relative_error = abs_error / max(abs_reference, 1.0e-10)
      write(*,*) 'Large error in function DASUM:'
      write(*,*) '  Central diff: ', central_diff
      write(*,*) '  AD result:   ', ad_result
      write(*,*) '  Absolute error:', abs_error
      write(*,*) '  Error bound:', error_bound
      write(*,*) '  Relative error:', relative_error
    end if
    ! Track max error for reporting (normalized)
    relative_error = abs_error / max(abs_reference, 1.0e-10)
    max_error = max(max_error, relative_error)
    
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_dasum