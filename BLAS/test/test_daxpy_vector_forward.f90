! Test program for DAXPY vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirsmax=4

program test_daxpy_vector_forward
  implicit none
  include 'DIFFSIZES.inc'

  external :: daxpy
  external :: daxpy_dv

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, idir  ! Loop counters
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  integer :: nsize
  real(8) :: da
  real(8), dimension(4) :: dx
  integer :: incx_val
  real(8), dimension(max_size) :: dy
  integer :: incy_val

  ! Vector mode derivative variables (type-promoted)
  ! Scalars become arrays(nbdirsmax), arrays gain extra dimension
  real(8), dimension(nbdirsmax) :: da_dv
  real(8), dimension(nbdirsmax,4) :: dx_dv
  real(8), dimension(nbdirsmax,max_size) :: dy_dv
  ! Declare variables for storing original values
  real(8) :: da_orig
  real(8), dimension(nbdirsmax) :: da_dv_orig
  real(8), dimension(4) :: dx_orig
  real(8), dimension(nbdirsmax,4) :: dx_dv_orig
  real(8), dimension(max_size) :: dy_orig
  real(8), dimension(nbdirsmax,max_size) :: dy_dv_orig

  ! Initialize test parameters
  nsize = n
  incx_val = 1
  incy_val = 1

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  seed_array = 42
  call random_seed(put=seed_array)

  call random_number(da)
  da = da * 2.0d0 - 1.0d0  ! Scale to [-1,1]
  call random_number(dx)
  dx = dx * 2.0d0 - 1.0d0  ! Scale to [-1,1]
  call random_number(dy)
  dy = dy * 2.0d0 - 1.0d0  ! Scale to [-1,1]

  ! Initialize input derivatives to random values (exactly like scalar mode)
  do idir = 1, nbdirsmax
    call random_number(temp_real)
    da_dv(idir) = temp_real * 2.0d0 - 1.0d0
  end do
  do idir = 1, nbdirsmax
    call random_number(dx_dv(idir,:))
    dx_dv(idir,:) = dx_dv(idir,:) * 2.0d0 - 1.0d0
  end do
  do idir = 1, nbdirsmax
    call random_number(dy_dv(idir,:))
    dy_dv(idir,:) = dy_dv(idir,:) * 2.0d0 - 1.0d0
  end do

  write(*,*) 'Testing DAXPY (Vector Forward Mode)'
  ! Store original values before any function calls (critical for INOUT parameters)
  da_orig = da
  da_dv_orig = da_dv
  dx_orig = dx
  dx_dv_orig = dx_dv
  dy_orig = dy
  dy_dv_orig = dy_dv

  ! Call the vector mode differentiated function

  call daxpy_dv(nsize, da, da_dv, dx, dx_dv, incx_val, dy, dy_dv, incy_val, nbdirsmax)

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
    integer :: i, j, idir
    logical :: has_large_errors
    real(8), dimension(max_size) :: dy_forward, dy_backward
    
    max_error = 0.0e0
    has_large_errors = .false.
    
    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Number of directions:', nbdirsmax
    
    ! Test each derivative direction separately
    do idir = 1, nbdirsmax
      
      ! Forward perturbation: f(x + h * direction)
      da = da_orig + h * da_dv_orig(idir)
      dx = dx_orig + h * dx_dv_orig(idir,:)
      dy = dy_orig + h * dy_dv_orig(idir,:)
      call daxpy(nsize, da, dx, incx_val, dy, incy_val)
      dy_forward = dy
      
      ! Backward perturbation: f(x - h * direction)
      da = da_orig - h * da_dv_orig(idir)
      dx = dx_orig - h * dx_dv_orig(idir,:)
      dy = dy_orig - h * dy_dv_orig(idir,:)
      call daxpy(nsize, da, dx, incx_val, dy, incy_val)
      dy_backward = dy
      
      ! Compute central differences and compare with AD results
      do i = 1, min(2, nsize)  ! Check only first few elements
        ! Central difference: (f(x+h) - f(x-h)) / (2h)
        central_diff = (dy_forward(i) - dy_backward(i)) / (2.0e0 * h)
        ! AD result
        ad_result = dy_dv(idir,i)
        ! Error check: |a - b| > atol + rtol * |b|
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 1.0e-5 + 1.0e-5 * abs_reference
        if (abs_error > error_bound) then
          has_large_errors = .true.
          relative_error = abs_error / max(abs_reference, 1.0e-10)
          write(*,*) '  Large error in direction', idir, ' output DY(', i, '):'
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

end program test_daxpy_vector_forward