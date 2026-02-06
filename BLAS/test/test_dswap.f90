! Test program for DSWAP differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision

program test_dswap
  implicit none

  external :: dswap
  external :: dswap_d

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension (rows/cols of matrices)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions

  integer :: nsize
  real(8), dimension(max_size) :: dx
  integer :: incx_val
  real(8), dimension(max_size) :: dy
  integer :: incy_val

  ! Derivative variables
  real(8), dimension(max_size) :: dx_d
  real(8), dimension(max_size) :: dy_d

  ! Storage variables for inout parameters
  real(8), dimension(max_size) :: dx_output
  real(8), dimension(max_size) :: dy_output

  ! Array restoration variables for numerical differentiation
  real(8), dimension(max_size) :: dy_orig
  real(8), dimension(max_size) :: dx_orig

  ! Variables for central difference computation
  real(8), dimension(max_size) :: dy_forward, dy_backward
  real(8), dimension(max_size) :: dx_forward, dx_backward
  ! Scalar variables for central difference computation
  real(8) :: central_diff, ad_result
  logical :: has_large_errors

  ! Variables for storing original derivative values
  real(8), dimension(max_size) :: dy_d_orig
  real(8), dimension(max_size) :: dx_d_orig

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
  call random_number(dy)
  dy = dy * 2.0d0 - 1.0d0  ! Scale to [-1,1]
  incy_val = 1

  ! Initialize input derivatives to random values
  call random_number(dy_d)
  dy_d = dy_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
  call random_number(dx_d)
  dx_d = dx_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]

  ! Store initial derivative values after random initialization
  dy_d_orig = dy_d
  dx_d_orig = dx_d

  ! Store original values for central difference computation
  dy_orig = dy
  dx_orig = dx

  write(*,*) 'Testing DSWAP'
  ! Store input values of inout parameters before first function call
  dx_orig = dx
  dy_orig = dy

  ! Re-initialize data for differentiated function
  ! Only reinitialize inout parameters - keep input-only parameters unchanged

  nsize = n
  dx = dx_orig
  incx_val = 1
  dy = dy_orig
  incy_val = 1

  ! Call the differentiated function
  call dswap_d(nsize, dx, dx_d, incx_val, dy, dy_d, incy_val)

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
    dy = dy_orig + h * dy_d_orig
    dx = dx_orig + h * dx_d_orig
    call dswap(nsize, dx, incx_val, dy, incy_val)
    ! Store forward perturbation results
    dy_forward = dy
    dx_forward = dx
    
    ! Backward perturbation: f(x - h)
    dy = dy_orig - h * dy_d_orig
    dx = dx_orig - h * dx_d_orig
    call dswap(nsize, dx, incx_val, dy, incy_val)
    ! Store backward perturbation results
    dy_backward = dy
    dx_backward = dx
    
    ! Compute central differences and compare with AD results
    ! Check derivatives for output DY
    do i = 1, min(2, n)  ! Check only first few elements
      ! Central difference: (f(x+h) - f(x-h)) / (2h)
      central_diff = (dy_forward(i) - dy_backward(i)) / (2.0e0 * h)
      ! AD result
      ad_result = dy_d(i)
      ! Error check: |a - b| > atol + rtol * |b|
      abs_error = abs(central_diff - ad_result)
      abs_reference = abs(ad_result)
      error_bound = 1.0e-5 + 1.0e-5 * abs_reference
      if (abs_error > error_bound) then
        has_large_errors = .true.
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        write(*,*) 'Large error in output DY(', i, '):'
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
    ! Check derivatives for output DX
    do i = 1, min(2, n)  ! Check only first few elements
      ! Central difference: (f(x+h) - f(x-h)) / (2h)
      central_diff = (dx_forward(i) - dx_backward(i)) / (2.0e0 * h)
      ! AD result
      ad_result = dx_d(i)
      ! Error check: |a - b| > atol + rtol * |b|
      abs_error = abs(central_diff - ad_result)
      abs_reference = abs(ad_result)
      error_bound = 1.0e-5 + 1.0e-5 * abs_reference
      if (abs_error > error_bound) then
        has_large_errors = .true.
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        write(*,*) 'Large error in output DX(', i, '):'
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
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_dswap