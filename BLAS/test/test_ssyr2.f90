! Test program for SSYR2 differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision

program test_ssyr2
  implicit none

  external :: ssyr2
  external :: ssyr2_d

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension (rows/cols of matrices)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions

  character :: uplo
  integer :: nsize
  real(4) :: alpha
  real(4), dimension(max_size) :: x
  integer :: incx_val
  real(4), dimension(max_size) :: y
  integer :: incy_val
  real(4), dimension(max_size,max_size) :: a
  integer :: lda_val

  ! Derivative variables
  real(4) :: alpha_d
  real(4), dimension(max_size) :: x_d
  real(4), dimension(max_size) :: y_d
  real(4), dimension(max_size,max_size) :: a_d

  ! Storage variables for inout parameters
  real(4), dimension(max_size,max_size) :: a_output

  ! Array restoration variables for numerical differentiation
  real(4) :: alpha_orig
  real(4), dimension(max_size,max_size) :: a_orig
  real(4), dimension(max_size) :: x_orig
  real(4), dimension(max_size) :: y_orig

  ! Variables for central difference computation
  real(4), dimension(max_size,max_size) :: a_forward, a_backward
  ! Scalar variables for central difference computation
  real(4) :: central_diff, ad_result
  logical :: has_large_errors

  ! Variables for storing original derivative values
  real(4) :: alpha_d_orig
  real(4), dimension(max_size,max_size) :: a_d_orig
  real(4), dimension(max_size) :: x_d_orig
  real(4), dimension(max_size) :: y_d_orig

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
  alpha = alpha * 2.0 - 1.0  ! Scale to [-1,1]
  call random_number(x)
  x = x * 2.0 - 1.0  ! Scale to [-1,1]
  incx_val = 1  ! INCX 1
  call random_number(y)
  y = y * 2.0 - 1.0  ! Scale to [-1,1]
  incy_val = 1  ! INCY 1
  call random_number(a)
  a = a * 2.0d0 - 1.0d0  ! Scale to [-1,1]
  lda_val = lda  ! LDA must be at least max( 1

  ! Initialize input derivatives to random values
  call random_number(alpha_d)
  alpha_d = alpha_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
  call random_number(a_d)
  a_d = a_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
  call random_number(x_d)
  x_d = x_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]
  call random_number(y_d)
  y_d = y_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]

  ! Store initial derivative values after random initialization
  alpha_d_orig = alpha_d
  a_d_orig = a_d
  x_d_orig = x_d
  y_d_orig = y_d

  ! Store original values for central difference computation
  alpha_orig = alpha
  a_orig = a
  x_orig = x
  y_orig = y

  write(*,*) 'Testing SSYR2'
  ! Store input values of inout parameters before first function call
  a_orig = a

  ! Re-initialize data for differentiated function
  ! Only reinitialize inout parameters - keep input-only parameters unchanged

  ! uplo already has correct value from original call
  nsize = n
  ! alpha already has correct value from original call
  ! x already has correct value from original call
  incx_val = 1  ! INCX 1
  ! y already has correct value from original call
  incy_val = 1  ! INCY 1
  a = a_orig
  lda_val = lda  ! LDA must be at least max( 1

  ! Call the differentiated function
  call ssyr2_d(uplo, nsize, alpha, alpha_d, x, x_d, incx_val, y, y_d, incy_val, a, a_d, lda_val)

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
    x = x_orig + h * x_d_orig
    y = y_orig + h * y_d_orig
    call ssyr2(uplo, nsize, alpha, x, incx_val, y, incy_val, a, lda_val)
    ! Store forward perturbation results
    a_forward = a
    
    ! Backward perturbation: f(x - h)
    alpha = alpha_orig - h * alpha_d_orig
    a = a_orig - h * a_d_orig
    x = x_orig - h * x_d_orig
    y = y_orig - h * y_d_orig
    call ssyr2(uplo, nsize, alpha, x, incx_val, y, incy_val, a, lda_val)
    ! Store backward perturbation results
    a_backward = a
    
    ! Compute central differences and compare with AD results
    ! Check derivatives for output A
    do j = 1, min(2, n)  ! Check only first few elements
      do i = 1, min(2, n)
        ! Central difference: (f(x+h) - f(x-h)) / (2h)
        central_diff = (a_forward(i,j) - a_backward(i,j)) / (2.0e0 * h)
        ! AD result
        ad_result = a_d(i,j)
        ! Error check: |a - b| > atol + rtol * |b|
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 2.0e-3 + 2.0e-3 * abs_reference
        if (abs_error > error_bound) then
          has_large_errors = .true.
          relative_error = abs_error / max(abs_reference, 1.0e-10)
          write(*,*) 'Large error in output A(', i, ',', j, '):'
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

end program test_ssyr2