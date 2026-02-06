! Test program for DSBMV vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirsmax=4

program test_dsbmv_vector_forward
  implicit none
  include 'DIFFSIZES.inc'

  external :: dsbmv
  external :: dsbmv_dv

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, idir, band_row  ! Loop counters
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  character :: uplo
  integer :: nsize
  integer :: ksize
  real(8) :: alpha
  real(8), dimension(max_size,max_size) :: a
  integer :: lda_val
  real(8), dimension(max_size) :: x
  integer :: incx_val
  real(8) :: beta
  real(8), dimension(max_size) :: y
  integer :: incy_val

  ! Vector mode derivative variables (type-promoted)
  ! Scalars become arrays(nbdirsmax), arrays gain extra dimension
  real(8), dimension(nbdirsmax) :: alpha_dv
  real(8), dimension(nbdirsmax,max_size,max_size) :: a_dv
  real(8), dimension(nbdirsmax,max_size) :: x_dv
  real(8), dimension(nbdirsmax) :: beta_dv
  real(8), dimension(nbdirsmax,max_size) :: y_dv
  ! Declare variables for storing original values
  real(8) :: alpha_orig
  real(8), dimension(nbdirsmax) :: alpha_dv_orig
  real(8), dimension(max_size,max_size) :: a_orig
  real(8), dimension(nbdirsmax,max_size,max_size) :: a_dv_orig
  real(8), dimension(max_size) :: x_orig
  real(8), dimension(nbdirsmax,max_size) :: x_dv_orig
  real(8) :: beta_orig
  real(8), dimension(nbdirsmax) :: beta_dv_orig
  real(8), dimension(max_size) :: y_orig
  real(8), dimension(nbdirsmax,max_size) :: y_dv_orig

  ! Initialize test parameters
  nsize = n
  ksize = max(0, n - 1)  ! Band width: 0 <= K <= N-1
  lda_val = lda
  incx_val = 1
  incy_val = 1

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  seed_array = 42
  call random_seed(put=seed_array)

  uplo = 'U'
  call random_number(alpha)
  alpha = alpha * 2.0d0 - 1.0d0  ! Scale to [-1,1]
  ! Initialize a as symmetric band matrix (upper band storage)
  ! A(band_row, j) = full(i,j) with band_row = ksize+1+i-j, i = max(1,j-ksize)..j
  do j = 1, n
    do band_row = max(1, ksize+2-j), ksize+1
      call random_number(temp_real)
      a(band_row, j) = temp_real * 2.0 - 1.0  ! Scale to [-1,1]
    end do
  end do
  call random_number(x)
  x = x * 2.0d0 - 1.0d0  ! Scale to [-1,1]
  call random_number(beta)
  beta = beta * 2.0d0 - 1.0d0  ! Scale to [-1,1]
  call random_number(y)
  y = y * 2.0d0 - 1.0d0  ! Scale to [-1,1]

  ! Initialize input derivatives to random values (exactly like scalar mode)
  do idir = 1, nbdirsmax
    call random_number(temp_real)
    alpha_dv(idir) = temp_real * 2.0d0 - 1.0d0
  end do
  do idir = 1, nbdirsmax
    call random_number(a_dv(idir,:,:))
    a_dv(idir,:,:) = a_dv(idir,:,:) * 2.0d0 - 1.0d0
  end do
  do idir = 1, nbdirsmax
    call random_number(x_dv(idir,:))
    x_dv(idir,:) = x_dv(idir,:) * 2.0d0 - 1.0d0
  end do
  do idir = 1, nbdirsmax
    call random_number(temp_real)
    beta_dv(idir) = temp_real * 2.0d0 - 1.0d0
  end do
  do idir = 1, nbdirsmax
    call random_number(y_dv(idir,:))
    y_dv(idir,:) = y_dv(idir,:) * 2.0d0 - 1.0d0
  end do

  write(*,*) 'Testing DSBMV (Vector Forward Mode)'
  ! Store original values before any function calls (critical for INOUT parameters)
  alpha_orig = alpha
  alpha_dv_orig = alpha_dv
  a_orig = a
  a_dv_orig = a_dv
  x_orig = x
  x_dv_orig = x_dv
  beta_orig = beta
  beta_dv_orig = beta_dv
  y_orig = y
  y_dv_orig = y_dv

  ! Call the vector mode differentiated function

  call dsbmv_dv(uplo, nsize, ksize, alpha, alpha_dv, a, a_dv, lda_val, x, x_dv, incx_val, beta, beta_dv, y, y_dv, incy_val, nbdirsmax)

  ! Print results and compare
  write(*,*) 'Function calls completed successfully'

  ! Numerical differentiation check
  call check_derivatives_numerically()

  write(*,*) 'Vector forward mode test completed successfully'

contains

  subroutine check_derivatives_numerically()
    implicit none
    real(8), parameter :: h = 1.0e-7  ! Step size for finite differences
    real(8) :: relative_error, max_error
    real(8) :: abs_error, abs_reference, error_bound
    real(8) :: central_diff, ad_result
    integer :: i, j, idir, band_row
    logical :: has_large_errors
    real(8), dimension(max_size) :: y_forward, y_backward
    
    max_error = 0.0e0
    has_large_errors = .false.
    
    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Number of directions:', nbdirsmax
    
    ! Test each derivative direction separately
    do idir = 1, nbdirsmax
      
      ! Forward perturbation: f(x + h * direction)
      alpha = alpha_orig + h * alpha_dv_orig(idir)
      a = a_orig + h * a_dv_orig(idir,:,:)
      x = x_orig + h * x_dv_orig(idir,:)
      beta = beta_orig + h * beta_dv_orig(idir)
      y = y_orig + h * y_dv_orig(idir,:)
      call dsbmv(uplo, nsize, ksize, alpha, a, lda_val, x, incx_val, beta, y, incy_val)
      y_forward = y
      
      ! Backward perturbation: f(x - h * direction)
      alpha = alpha_orig - h * alpha_dv_orig(idir)
      a = a_orig - h * a_dv_orig(idir,:,:)
      x = x_orig - h * x_dv_orig(idir,:)
      beta = beta_orig - h * beta_dv_orig(idir)
      y = y_orig - h * y_dv_orig(idir,:)
      call dsbmv(uplo, nsize, ksize, alpha, a, lda_val, x, incx_val, beta, y, incy_val)
      y_backward = y
      
      ! Compute central differences and compare with AD results
      do i = 1, min(2, nsize)  ! Check only first few elements
        ! Central difference: (f(x+h) - f(x-h)) / (2h)
        central_diff = (y_forward(i) - y_backward(i)) / (2.0e0 * h)
        ! AD result
        ad_result = y_dv(idir,i)
        ! Error check: |a - b| > atol + rtol * |b|
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 1.0e-5 + 1.0e-5 * abs_reference
        if (abs_error > error_bound) then
          has_large_errors = .true.
          relative_error = abs_error / max(abs_reference, 1.0e-10)
          write(*,*) '  Large error in direction', idir, ' output Y(', i, '):'
          write(*,*) '    Central diff: ', central_diff
          write(*,*) '    AD result:   ', ad_result
          write(*,*) '    Absolute error:', abs_error
          write(*,*) '    Error bound:', error_bound
          write(*,*) '    Relative error:', relative_error
        end if
        ! Track max error for reporting (normalized)
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        max_error = max(max_error, relative_error)
      end do
    end do
    
    write(*,*) 'Maximum relative error across all directions:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in vector derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Vector derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_dsbmv_vector_forward