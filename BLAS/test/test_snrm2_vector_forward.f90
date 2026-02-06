! Test program for SNRM2 vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirsmax=4

program test_snrm2_vector_forward
  use DIFFSIZES
  implicit none

  real(4), external :: snrm2
  external :: snrm2_dv

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, idir  ! Loop counters
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  integer :: nsize
  real(4), dimension(4) :: x
  integer :: incx_val

  ! Vector mode derivative variables (type-promoted)
  ! Scalars become arrays(nbdirsmax), arrays gain extra dimension
  real(4), dimension(nbdirsmax,4) :: x_dv
  ! Declare variables for storing original values
  real(4), dimension(4) :: x_orig
  real(4), dimension(nbdirsmax,4) :: x_dv_orig

  ! Function result variables
  real(4) :: snrm2_result
  real(4), dimension(nbdirsmax) :: snrm2_dv_result

  ! Initialize test parameters
  nsize = n
  incx_val = 1

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  seed_array = 42
  call random_seed(put=seed_array)

  call random_number(x)
  x = x * 2.0 - 1.0  ! Scale to [-1,1]

  ! Initialize input derivatives to random values (exactly like scalar mode)
  do idir = 1, nbdirsmax
    call random_number(x_dv(idir,:))
    x_dv(idir,:) = x_dv(idir,:) * 2.0 - 1.0
  end do

  write(*,*) 'Testing SNRM2 (Vector Forward Mode)'
  ! Store original values before any function calls (critical for INOUT parameters)
  x_orig = x
  x_dv_orig = x_dv

  ! Call the vector mode differentiated function

  call snrm2_dv(nsize, x, x_dv, incx_val, snrm2_result, snrm2_dv_result, nbdirsmax)

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
    real(4) :: snrm2_forward, snrm2_backward
    
    max_error = 0.0e0
    has_large_errors = .false.
    
    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Number of directions:', nbdirsmax
    
    ! Test each derivative direction separately
    do idir = 1, nbdirsmax
      
      ! Forward perturbation: f(x + h * direction)
      x = x_orig + h * x_dv_orig(idir,:)
      snrm2_forward = snrm2(nsize, x, incx_val)
      
      ! Backward perturbation: f(x - h * direction)
      x = x_orig - h * x_dv_orig(idir,:)
      snrm2_backward = snrm2(nsize, x, incx_val)
      
      ! Compute central differences and compare with AD results
      ! Central difference: (f(x+h) - f(x-h)) / (2h)
      central_diff = (snrm2_forward - snrm2_backward) / (2.0e0 * h)
      ! AD result
      ad_result = snrm2_dv_result(idir)
      ! Error check: |a - b| > atol + rtol * |b|
      abs_error = abs(central_diff - ad_result)
      abs_reference = abs(ad_result)
      error_bound = 2.0e-3 + 2.0e-3 * abs_reference
      if (abs_error > error_bound) then
        has_large_errors = .true.
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        write(*,*) '  Large error in direction', idir, ' output SNRM2:'
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
    
    write(*,*) 'Maximum relative error across all directions:', max_error
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in vector derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Vector derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_snrm2_vector_forward