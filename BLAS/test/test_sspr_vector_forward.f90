! Test program for SSPR vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirsmax=4

program test_sspr_vector_forward
  implicit none
  include 'DIFFSIZES.inc'

  external :: sspr
  external :: sspr_dv

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
  real(4), dimension(max_size) :: x
  integer :: incx_val
  real(4), dimension((n*(n+1))/2) :: ap

  ! Vector mode derivative variables (type-promoted)
  ! Scalars become arrays(nbdirsmax), arrays gain extra dimension
  real(4), dimension(nbdirsmax) :: alpha_dv
  real(4), dimension(nbdirsmax,max_size) :: x_dv
  real(4), dimension(nbdirsmax,(n*(n+1))/2) :: ap_dv
  ! Declare variables for storing original values
  real(4) :: alpha_orig
  real(4), dimension(nbdirsmax) :: alpha_dv_orig
  real(4), dimension(max_size) :: x_orig
  real(4), dimension(nbdirsmax,max_size) :: x_dv_orig
  real(4), dimension((n*(n+1))/2) :: ap_orig
  real(4), dimension(nbdirsmax,(n*(n+1))/2) :: ap_dv_orig

  ! Initialize test parameters
  nsize = n
  incx_val = 1

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  seed_array = 42
  call random_seed(put=seed_array)

  uplo = 'U'
  call random_number(alpha)
  alpha = alpha * 2.0 - 1.0  ! Scale to [-1,1]
  call random_number(x)
  x = x * 2.0 - 1.0  ! Scale to [-1,1]
  call random_number(ap)
  ap = ap * 2.0 - 1.0  ! Scale to [-1,1]

  ! Initialize input derivatives to random values (exactly like scalar mode)
  do idir = 1, nbdirsmax
    call random_number(temp_real)
    alpha_dv(idir) = temp_real * 2.0 - 1.0
  end do
  do idir = 1, nbdirsmax
    call random_number(x_dv(idir,:))
    x_dv(idir,:) = x_dv(idir,:) * 2.0 - 1.0
  end do
  do idir = 1, nbdirsmax
    call random_number(ap_dv(idir,:))
    ap_dv(idir,:) = ap_dv(idir,:) * 2.0d0 - 1.0d0
  end do

  write(*,*) 'Testing SSPR (Vector Forward Mode)'
  ! Store original values before any function calls (critical for INOUT parameters)
  alpha_orig = alpha
  alpha_dv_orig = alpha_dv
  x_orig = x
  x_dv_orig = x_dv
  ap_orig = ap
  ap_dv_orig = ap_dv

  ! Call the vector mode differentiated function

  call sspr_dv(uplo, nsize, alpha, alpha_dv, x, x_dv, incx_val, ap, ap_dv, nbdirsmax)

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
    real(4), dimension((n*(n+1))/2) :: ap_forward, ap_backward
    
    max_error = 0.0e0
    has_large_errors = .false.
    
    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Number of directions:', nbdirsmax
    
    ! Test each derivative direction separately
    do idir = 1, nbdirsmax
      
      ! Forward perturbation: f(x + h * direction)
      alpha = alpha_orig + h * alpha_dv_orig(idir)
      x = x_orig + h * x_dv_orig(idir,:)
      ap = ap_orig + h * ap_dv_orig(idir,:)
      call sspr(uplo, nsize, alpha, x, incx_val, ap)
      ap_forward = ap
      
      ! Backward perturbation: f(x - h * direction)
      alpha = alpha_orig - h * alpha_dv_orig(idir)
      x = x_orig - h * x_dv_orig(idir,:)
      ap = ap_orig - h * ap_dv_orig(idir,:)
      call sspr(uplo, nsize, alpha, x, incx_val, ap)
      ap_backward = ap
      
      ! Compute central differences and compare with AD results
      do i = 1, min(2, nsize)  ! Check only first few elements
        ! Central difference: (f(x+h) - f(x-h)) / (2h)
        central_diff = (ap_forward(i) - ap_backward(i)) / (2.0e0 * h)
        ! AD result
        ad_result = ap_dv(idir,i)
        ! Error check: |a - b| > atol + rtol * |b|
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 2.0e-3 + 2.0e-3 * abs_reference
        if (abs_error > error_bound) then
          has_large_errors = .true.
          relative_error = abs_error / max(abs_reference, 1.0e-10)
          write(*,*) '  Large error in direction', idir, ' output AP(', i, '):'
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

end program test_sspr_vector_forward