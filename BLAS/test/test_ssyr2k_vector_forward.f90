! Test program for SSYR2K vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirsmax=4

program test_ssyr2k_vector_forward
  implicit none
  include 'DIFFSIZES.inc'

  external :: ssyr2k
  external :: ssyr2k_dv

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, idir  ! Loop counters
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  character :: uplo
  character :: trans
  integer :: nsize
  integer :: ksize
  real(4) :: alpha
  real(4), dimension(max_size,max_size) :: a
  integer :: lda_val
  real(4), dimension(max_size,max_size) :: b
  integer :: ldb_val
  real(4) :: beta
  real(4), dimension(max_size,max_size) :: c
  integer :: ldc_val

  ! Vector mode derivative variables (type-promoted)
  ! Scalars become arrays(nbdirsmax), arrays gain extra dimension
  real(4), dimension(nbdirsmax) :: alpha_dv
  real(4), dimension(nbdirsmax,max_size,max_size) :: a_dv
  real(4), dimension(nbdirsmax,max_size,max_size) :: b_dv
  real(4), dimension(nbdirsmax) :: beta_dv
  real(4), dimension(nbdirsmax,max_size,max_size) :: c_dv
  ! Declare variables for storing original values
  real(4) :: alpha_orig
  real(4), dimension(nbdirsmax) :: alpha_dv_orig
  real(4), dimension(max_size,max_size) :: a_orig
  real(4), dimension(nbdirsmax,max_size,max_size) :: a_dv_orig
  real(4), dimension(max_size,max_size) :: b_orig
  real(4), dimension(nbdirsmax,max_size,max_size) :: b_dv_orig
  real(4) :: beta_orig
  real(4), dimension(nbdirsmax) :: beta_dv_orig
  real(4), dimension(max_size,max_size) :: c_orig
  real(4), dimension(nbdirsmax,max_size,max_size) :: c_dv_orig

  ! Initialize test parameters
  nsize = n
  ksize = n
  lda_val = lda
  ldb_val = ldb
  ldc_val = ldc

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  seed_array = 42
  call random_seed(put=seed_array)

  uplo = 'U'
  trans = 'N'
  call random_number(alpha)
  alpha = alpha * 2.0 - 1.0  ! Scale to [-1,1]
  call random_number(a)
  a = a * 2.0 - 1.0  ! Scale to [-1,1]
  call random_number(b)
  b = b * 2.0 - 1.0  ! Scale to [-1,1]
  call random_number(beta)
  beta = beta * 2.0 - 1.0  ! Scale to [-1,1]
  call random_number(c)
  c = c * 2.0 - 1.0  ! Scale to [-1,1]

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
    call random_number(b_dv(idir,:,:))
    b_dv(idir,:,:) = b_dv(idir,:,:) * 2.0 - 1.0
  end do
  do idir = 1, nbdirsmax
    call random_number(temp_real)
    beta_dv(idir) = temp_real * 2.0 - 1.0
  end do
  do idir = 1, nbdirsmax
    call random_number(c_dv(idir,:,:))
    c_dv(idir,:,:) = c_dv(idir,:,:) * 2.0 - 1.0
  end do

  write(*,*) 'Testing SSYR2K (Vector Forward Mode)'
  ! Store original values before any function calls (critical for INOUT parameters)
  alpha_orig = alpha
  alpha_dv_orig = alpha_dv
  a_orig = a
  a_dv_orig = a_dv
  b_orig = b
  b_dv_orig = b_dv
  beta_orig = beta
  beta_dv_orig = beta_dv
  c_orig = c
  c_dv_orig = c_dv

  ! Call the vector mode differentiated function

  call ssyr2k_dv(uplo, trans, nsize, ksize, alpha, alpha_dv, a, a_dv, lda_val, b, b_dv, ldb_val, beta, beta_dv, c, c_dv, ldc_val, nbdirsmax)

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
    real(4), dimension(max_size,max_size) :: c_forward, c_backward
    
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
      b = b_orig + h * b_dv_orig(idir,:,:)
      beta = beta_orig + h * beta_dv_orig(idir)
      c = c_orig + h * c_dv_orig(idir,:,:)
      call ssyr2k(uplo, trans, nsize, ksize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
      c_forward = c
      
      ! Backward perturbation: f(x - h * direction)
      alpha = alpha_orig - h * alpha_dv_orig(idir)
      a = a_orig - h * a_dv_orig(idir,:,:)
      b = b_orig - h * b_dv_orig(idir,:,:)
      beta = beta_orig - h * beta_dv_orig(idir)
      c = c_orig - h * c_dv_orig(idir,:,:)
      call ssyr2k(uplo, trans, nsize, ksize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
      c_backward = c
      
      ! Compute central differences and compare with AD results
      do j = 1, min(2, nsize)  ! Check only first few elements
        do i = 1, min(2, nsize)
          ! Central difference: (f(x+h) - f(x-h)) / (2h)
          central_diff = (c_forward(i,j) - c_backward(i,j)) / (2.0e0 * h)
          ! AD result
          ad_result = c_dv(idir,i,j)
          ! Error check: |a - b| > atol + rtol * |b|
          abs_error = abs(central_diff - ad_result)
          abs_reference = abs(ad_result)
          error_bound = 2.0e-3 + 2.0e-3 * abs_reference
          if (abs_error > error_bound) then
            has_large_errors = .true.
            relative_error = abs_error / max(abs_reference, 1.0e-10)
            write(*,*) '  Large error in direction', idir, ' output C(', i, ',', j, '):'
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
    end do
    
    write(*,*) 'Maximum relative error across all directions:', max_error
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in vector derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Vector derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_ssyr2k_vector_forward