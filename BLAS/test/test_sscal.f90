! Test program for SSCAL differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision

program test_sscal
  implicit none

  external :: sscal
  external :: sscal_d

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension (rows/cols of matrices)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions

  integer :: nsize
  real(4) :: sa
  real(4), dimension(max_size) :: sx
  integer :: incx_val

  ! Derivative variables
  real(4) :: sa_d
  real(4), dimension(max_size) :: sx_d

  ! Storage variables for inout parameters
  real(4), dimension(max_size) :: sx_output

  ! Array restoration variables for numerical differentiation
  real(4) :: sa_orig
  real(4), dimension(max_size) :: sx_orig

  ! Variables for central difference computation
  real(4), dimension(max_size) :: sx_forward, sx_backward
  ! Scalar variables for central difference computation
  real(4) :: central_diff, ad_result
  logical :: has_large_errors

  ! Variables for storing original derivative values
  real(4) :: sa_d_orig
  real(4), dimension(max_size) :: sx_d_orig

  ! Temporary variables for matrix initialization
  real(4) :: temp_real, temp_imag
  integer :: i, j

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  integer :: seed_array(33)
  seed_array = 42
  call random_seed(put=seed_array)

  nsize = n
  call random_number(sa)
  sa = sa * 2.0 - 1.0  ! Scale to [-1,1]
  call random_number(sx)
  sx = sx * 2.0 - 1.0  ! Scale to [-1,1]
  incx_val = 1

  ! Initialize input derivatives to random values
  call random_number(sa_d)
  sa_d = sa_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
  call random_number(sx_d)
  sx_d = sx_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]

  ! Store initial derivative values after random initialization
  sa_d_orig = sa_d
  sx_d_orig = sx_d

  ! Store original values for central difference computation
  sa_orig = sa
  sx_orig = sx

  write(*,*) 'Testing SSCAL'
  ! Store input values of inout parameters before first function call
  sx_orig = sx

  ! Re-initialize data for differentiated function
  ! Only reinitialize inout parameters - keep input-only parameters unchanged

  nsize = n
  ! sa already has correct value from original call
  sx = sx_orig
  incx_val = 1

  ! Call the differentiated function
  call sscal_d(nsize, sa, sa_d, sx, sx_d, incx_val)

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
    sa = sa_orig + h * sa_d_orig
    sx = sx_orig + h * sx_d_orig
    call sscal(nsize, sa, sx, incx_val)
    ! Store forward perturbation results
    sx_forward = sx
    
    ! Backward perturbation: f(x - h)
    sa = sa_orig - h * sa_d_orig
    sx = sx_orig - h * sx_d_orig
    call sscal(nsize, sa, sx, incx_val)
    ! Store backward perturbation results
    sx_backward = sx
    
    ! Compute central differences and compare with AD results
    ! Check derivatives for output SX
    do i = 1, min(2, n)  ! Check only first few elements
      ! Central difference: (f(x+h) - f(x-h)) / (2h)
      central_diff = (sx_forward(i) - sx_backward(i)) / (2.0e0 * h)
      ! AD result
      ad_result = sx_d(i)
      ! Error check: |a - b| > atol + rtol * |b|
      abs_error = abs(central_diff - ad_result)
      abs_reference = abs(ad_result)
      error_bound = 2.0e-3 + 2.0e-3 * abs_reference
      if (abs_error > error_bound) then
        has_large_errors = .true.
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        write(*,*) 'Large error in output SX(', i, '):'
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
    
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_sscal