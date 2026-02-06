! Test program for SNRM2 differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision

program test_snrm2
  implicit none

  real(4), external :: snrm2
  real(4), external :: snrm2_d

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension (rows/cols of matrices)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions

  integer :: nsize
  real(4), dimension(4) :: x
  integer :: incx_val

  ! Derivative variables
  real(4), dimension(4) :: x_d

  ! Storage variables for inout parameters

  ! Array restoration variables for numerical differentiation
  real(4), dimension(4) :: x_orig
  real(4) :: snrm2_orig

  ! Variables for central difference computation
  ! Scalar variables for central difference computation
  real(4) :: central_diff, ad_result
  logical :: has_large_errors
  real(4) :: snrm2_result, snrm2_d_result
  real(4) :: snrm2_forward, snrm2_backward

  ! Variables for storing original derivative values
  real(4), dimension(4) :: x_d_orig

  ! Temporary variables for matrix initialization
  real(4) :: temp_real, temp_imag
  integer :: i, j

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  integer :: seed_array(33)
  seed_array = 42
  call random_seed(put=seed_array)

  nsize = n
  call random_number(x)
  x = x * 2.0 - 1.0  ! Scale to [-1,1]
  incx_val = 1

  ! Initialize input derivatives to random values
  call random_number(x_d)
  x_d = x_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]

  ! Store initial derivative values after random initialization
  x_d_orig = x_d

  ! Store original values for central difference computation
  x_orig = x

  write(*,*) 'Testing SNRM2'
  ! Store input values of inout parameters before first function call

  ! Call the original function
  snrm2_result = snrm2(nsize, x, incx_val)

  ! Store output values of inout parameters after first function call

  ! Re-initialize data for differentiated function
  ! Only reinitialize inout parameters - keep input-only parameters unchanged

  nsize = n
  ! x already has correct value from original call
  incx_val = 1

  ! Call the differentiated function
  snrm2_d_result = snrm2_d(nsize, x, x_d, incx_val, snrm2_result)

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
    x = x_orig + h * x_d_orig
    snrm2_forward = snrm2(nsize, x, incx_val)
    ! Store forward perturbation results
    ! snrm2_forward already captured above
    
    ! Backward perturbation: f(x - h)
    x = x_orig - h * x_d_orig
    snrm2_backward = snrm2(nsize, x, incx_val)
    ! Store backward perturbation results
    ! snrm2_backward already captured above
    
    ! Compute central differences and compare with AD results
    ! Check derivatives for function SNRM2
    ! Central difference: (f(x+h) - f(x-h)) / (2h)
    central_diff = (snrm2_forward - snrm2_backward) / (2.0e0 * h)
    ! AD result
    ad_result = snrm2_d_result
    ! Error check: |a - b| > atol + rtol * |b|
    abs_error = abs(central_diff - ad_result)
    abs_reference = abs(ad_result)
    error_bound = 2.0e-3 + 2.0e-3 * abs_reference
    if (abs_error > error_bound) then
      has_large_errors = .true.
      relative_error = abs_error / max(abs_reference, 1.0e-10)
      write(*,*) 'Large error in function SNRM2:'
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
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_snrm2