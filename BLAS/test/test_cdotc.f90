! Test program for CDOTC differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision

program test_cdotc
  implicit none

  complex(4), external :: cdotc
  complex(4), external :: cdotc_d

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension (rows/cols of matrices)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions

  integer :: nsize
  complex(4), dimension(4) :: cx
  integer :: incx_val
  complex(4), dimension(4) :: cy
  integer :: incy_val

  ! Derivative variables
  complex(4), dimension(4) :: cx_d
  complex(4), dimension(4) :: cy_d

  ! Storage variables for inout parameters

  ! Array restoration variables for numerical differentiation
  complex(4), dimension(4) :: cx_orig
  complex(4), dimension(4) :: cy_orig
  complex(4) :: cdotc_orig

  ! Variables for central difference computation
  ! Scalar variables for central difference computation
  complex(4) :: central_diff, ad_result
  logical :: has_large_errors
  complex(4) :: cdotc_result, cdotc_d_result
  complex(4) :: cdotc_forward, cdotc_backward

  ! Variables for storing original derivative values
  complex(4), dimension(4) :: cx_d_orig
  complex(4), dimension(4) :: cy_d_orig

  ! Temporary variables for matrix initialization
  real(4) :: temp_real, temp_imag
  integer :: i, j

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  integer :: seed_array(33)
  seed_array = 42
  call random_seed(put=seed_array)

  nsize = n
  do i = 1, n
    call random_number(temp_real)
    call random_number(temp_imag)
    cx(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  end do
  incx_val = 1
  do i = 1, n
    call random_number(temp_real)
    call random_number(temp_imag)
    cy(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  end do
  incy_val = 1

  ! Initialize input derivatives to random values
  do i = 1, n
    call random_number(temp_real)
    call random_number(temp_imag)
    cx_d(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  end do
  do i = 1, n
    call random_number(temp_real)
    call random_number(temp_imag)
    cy_d(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  end do

  ! Store initial derivative values after random initialization
  cx_d_orig = cx_d
  cy_d_orig = cy_d

  ! Store original values for central difference computation
  cx_orig = cx
  cy_orig = cy

  write(*,*) 'Testing CDOTC'
  ! Store input values of inout parameters before first function call

  ! Call the original function
  cdotc_result = cdotc(nsize, cx, incx_val, cy, incy_val)

  ! Store output values of inout parameters after first function call

  ! Re-initialize data for differentiated function
  ! Only reinitialize inout parameters - keep input-only parameters unchanged

  nsize = n
  ! cx already has correct value from original call
  incx_val = 1
  ! cy already has correct value from original call
  incy_val = 1

  ! Call the differentiated function
  cdotc_d_result = cdotc_d(nsize, cx, cx_d, incx_val, cy, cy_d, incy_val, cdotc_result)

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
    
    ! Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3
    
    ! Original values already stored in main program
    
    ! Central difference computation: f(x + h) - f(x - h) / (2h)
    ! Forward perturbation: f(x + h)
    cx = cx_orig + cmplx(h, 0.0) * cx_d_orig
    cy = cy_orig + cmplx(h, 0.0) * cy_d_orig
    cdotc_forward = cdotc(nsize, cx, incx_val, cy, incy_val)
    ! Store forward perturbation results
    ! cdotc_forward already captured above
    
    ! Backward perturbation: f(x - h)
    cx = cx_orig - cmplx(h, 0.0) * cx_d_orig
    cy = cy_orig - cmplx(h, 0.0) * cy_d_orig
    cdotc_backward = cdotc(nsize, cx, incx_val, cy, incy_val)
    ! Store backward perturbation results
    ! cdotc_backward already captured above
    
    ! Compute central differences and compare with AD results
    ! Check derivatives for function CDOTC
    ! Central difference: (f(x+h) - f(x-h)) / (2h)
    central_diff = (cdotc_forward - cdotc_backward) / (2.0e0 * h)
    ! AD result
    ad_result = cdotc_d_result
    ! Error check: |a - b| > atol + rtol * |b|
    abs_error = abs(central_diff - ad_result)
    abs_reference = abs(ad_result)
    error_bound = 1.0e-3 + 1.0e-3 * abs_reference
    if (abs_error > error_bound) then
      has_large_errors = .true.
      relative_error = abs_error / max(abs_reference, 1.0e-10)
      write(*,*) 'Large error in function CDOTC:'
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
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_cdotc