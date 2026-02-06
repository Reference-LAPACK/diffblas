! Test program for STBMV differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision

program test_stbmv
  implicit none

  external :: stbmv
  external :: stbmv_d

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension (rows/cols of matrices)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions

  character :: uplo
  character :: trans
  character :: diag
  integer :: nsize
  integer :: ksize
  real(4), dimension(max_size,n) :: a  ! Band storage (k+1) x n
  integer :: lda_val
  real(4), dimension(max_size) :: x
  integer :: incx_val

  ! Derivative variables
  real(4), dimension(max_size,max_size) :: a_d
  real(4), dimension(max_size) :: x_d

  ! Storage variables for inout parameters
  real(4), dimension(max_size) :: x_output

  ! Array restoration variables for numerical differentiation
  real(4), dimension(max_size,n) :: a_orig  ! Band storage
  real(4), dimension(max_size) :: x_orig

  ! Variables for central difference computation
  real(4), dimension(max_size) :: x_forward, x_backward
  ! Scalar variables for central difference computation
  real(4) :: central_diff, ad_result
  logical :: has_large_errors

  ! Variables for storing original derivative values
  real(4), dimension(max_size,max_size) :: a_d_orig
  real(4), dimension(max_size) :: x_d_orig

  ! Temporary variables for matrix initialization
  real(4) :: temp_real, temp_imag
  integer :: i, j, band_row

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  integer :: seed_array(33)
  seed_array = 42
  call random_seed(put=seed_array)

  uplo = 'U'
  trans = 'N'
  diag = 'N'
  nsize = n
  ksize = max(0, n - 1)  ! Band width: 0 <= K <= N-1
  ! Initialize a as triangular band matrix (upper band storage)
  ! A(band_row, j) = full(i,j) with band_row = ksize+1+i-j, i = max(1,j-ksize)..j
  do j = 1, n
    do band_row = max(1, ksize+2-j), ksize+1
      call random_number(temp_real)
      a(band_row, j) = temp_real * 2.0 - 1.0  ! Scale to [-1,1]
    end do
  end do
  lda_val = lda  ! LDA must be at least ( k + 1 )
  call random_number(x)
  x = x * 2.0 - 1.0  ! Scale to [-1,1]
  incx_val = 1  ! INCX 1

  ! Initialize input derivatives to random values
  ! Initialize a_d as triangular band matrix (upper band storage)
  ! A(band_row, j) = full(i,j) with band_row = ksize+1+i-j, i = max(1,j-ksize)..j
  do j = 1, n
    do band_row = max(1, ksize+2-j), ksize+1
      call random_number(temp_real)
      a_d(band_row, j) = temp_real * 2.0 - 1.0  ! Scale to [-1,1]
    end do
  end do
  call random_number(x_d)
  x_d = x_d * 2.0e0 - 1.0e0  ! Scale to [-1,1]

  ! Store initial derivative values after random initialization
  a_d_orig = a_d
  x_d_orig = x_d

  ! Store original values for central difference computation
  a_orig = a
  x_orig = x

  write(*,*) 'Testing STBMV'
  ! Store input values of inout parameters before first function call
  x_orig = x

  ! Re-initialize data for differentiated function
  ! Only reinitialize inout parameters - keep input-only parameters unchanged

  ! uplo already has correct value from original call
  ! trans already has correct value from original call
  ! diag already has correct value from original call
  nsize = n
  ksize = max(0, n - 1)  ! Band width: 0 <= K <= N-1
  ! a already has correct value from original call
  lda_val = lda  ! LDA must be at least ( k + 1 )
  x = x_orig
  incx_val = 1  ! INCX 1

  ! Call the differentiated function
  call stbmv_d(uplo, trans, diag, nsize, ksize, a, a_d, lda_val, x, x_d, incx_val)

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
    a = a_orig + h * a_d_orig
    x = x_orig + h * x_d_orig
    call stbmv(uplo, trans, diag, nsize, ksize, a, lda_val, x, incx_val)
    ! Store forward perturbation results
    x_forward = x
    
    ! Backward perturbation: f(x - h)
    a = a_orig - h * a_d_orig
    x = x_orig - h * x_d_orig
    call stbmv(uplo, trans, diag, nsize, ksize, a, lda_val, x, incx_val)
    ! Store backward perturbation results
    x_backward = x
    
    ! Compute central differences and compare with AD results
    ! Check derivatives for output X
    do i = 1, min(2, n)  ! Check only first few elements
      ! Central difference: (f(x+h) - f(x-h)) / (2h)
      central_diff = (x_forward(i) - x_backward(i)) / (2.0e0 * h)
      ! AD result
      ad_result = x_d(i)
      ! Error check: |a - b| > atol + rtol * |b|
      abs_error = abs(central_diff - ad_result)
      abs_reference = abs(ad_result)
      error_bound = 2.0e-3 + 2.0e-3 * abs_reference
      if (abs_error > error_bound) then
        has_large_errors = .true.
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        write(*,*) 'Large error in output X(', i, '):'
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

end program test_stbmv