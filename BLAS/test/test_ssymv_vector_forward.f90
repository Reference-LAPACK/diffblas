! Test program for SSYMV vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirsmax=4

program test_ssymv_vector_forward
  implicit none
  include 'DIFFSIZES.inc'

  external :: ssymv
  external :: ssymv_dv

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, idir  ! Loop counters
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  character :: uplo
  integer :: nsize
  real(4) :: alpha
  real(4), dimension(max_size,max_size) :: a
  integer :: lda_val
  real(4), dimension(max_size) :: x
  integer :: incx_val
  real(4) :: beta
  real(4), dimension(max_size) :: y
  integer :: incy_val

  ! Vector mode derivative variables (type-promoted)
  ! Scalars become arrays(nbdirsmax), arrays gain extra dimension
  real(4), dimension(nbdirsmax) :: alpha_dv
  real(4), dimension(nbdirsmax,max_size,max_size) :: a_dv
  real(4), dimension(nbdirsmax,max_size) :: x_dv
  real(4), dimension(nbdirsmax) :: beta_dv
  real(4), dimension(nbdirsmax,max_size) :: y_dv
  ! Declare variables for storing original values
  real(4) :: alpha_orig
  real(4), dimension(nbdirsmax) :: alpha_dv_orig
  real(4), dimension(max_size,max_size) :: a_orig
  real(4), dimension(nbdirsmax,max_size,max_size) :: a_dv_orig
  real(4), dimension(max_size) :: x_orig
  real(4), dimension(nbdirsmax,max_size) :: x_dv_orig
  real(4) :: beta_orig
  real(4), dimension(nbdirsmax) :: beta_dv_orig
  real(4), dimension(max_size) :: y_orig
  real(4), dimension(nbdirsmax,max_size) :: y_dv_orig

  ! Initialize test parameters
  nsize = n
  lda_val = lda
  incx_val = 1
  incy_val = 1

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  seed_array = 42
  call random_seed(put=seed_array)

  uplo = 'U'
  call random_number(alpha)
  alpha = alpha * 2.0 - 1.0  ! Scale to [-1,1]
  call random_number(a)
  a = a * 2.0 - 1.0  ! Scale to [-1,1]
  call random_number(x)
  x = x * 2.0 - 1.0  ! Scale to [-1,1]
  call random_number(beta)
  beta = beta * 2.0 - 1.0  ! Scale to [-1,1]
  call random_number(y)
  y = y * 2.0 - 1.0  ! Scale to [-1,1]

  ! Initialize input derivatives to random values (exactly like scalar mode)
  do idir = 1, nbdirsmax
    call random_number(temp_real)
    alpha_dv(idir) = temp_real * 2.0 - 1.0
  end do
  do idir = 1, nbdirsmax
    call random_number(a_dv(idir,:,:))
    a_dv(idir,:,:) = a_dv(idir,:,:) * 2.0 - 1.0
  end do
  do idir = 1, nbdirsmax
    call random_number(x_dv(idir,:))
    x_dv(idir,:) = x_dv(idir,:) * 2.0 - 1.0
  end do
  do idir = 1, nbdirsmax
    call random_number(temp_real)
    beta_dv(idir) = temp_real * 2.0 - 1.0
  end do
  do idir = 1, nbdirsmax
    call random_number(y_dv(idir,:))
    y_dv(idir,:) = y_dv(idir,:) * 2.0 - 1.0
  end do

  write(*,*) 'Testing SSYMV (Vector Forward Mode)'
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

  call ssymv_dv(uplo, nsize, alpha, alpha_dv, a, a_dv, lda_val, x, x_dv, incx_val, beta, beta_dv, y, y_dv, incy_val, nbdirsmax)

  ! Print results and compare
  write(*,*) 'Function calls completed successfully'

  ! Numerical differentiation check
  call check_derivatives_numerically()

  write(*,*) 'Vector forward mode test completed successfully'

contains

  subroutine check_derivatives_numerically()
    implicit none
    real(4), parameter :: h = 1.0e-3  ! Step size for finite differences
    real(4) :: relative_error, max_error
    real(4) :: abs_error, abs_reference, error_bound
    real(4) :: central_diff, ad_result
    integer :: i, j, idir
    logical :: has_large_errors
    real(4), dimension(max_size) :: y_forward, y_backward
    
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
      call ssymv(uplo, nsize, alpha, a, lda_val, x, incx_val, beta, y, incy_val)
      y_forward = y
      
      ! Backward perturbation: f(x - h * direction)
      alpha = alpha_orig - h * alpha_dv_orig(idir)
      a = a_orig - h * a_dv_orig(idir,:,:)
      x = x_orig - h * x_dv_orig(idir,:)
      beta = beta_orig - h * beta_dv_orig(idir)
      y = y_orig - h * y_dv_orig(idir,:)
      call ssymv(uplo, nsize, alpha, a, lda_val, x, incx_val, beta, y, incy_val)
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
        error_bound = 2.0e-3 + 2.0e-3 * abs_reference
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
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in vector derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Vector derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_ssymv_vector_forward