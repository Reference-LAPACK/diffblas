! Test program for CCOPY vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirsmax=4

program test_ccopy_vector_forward
  implicit none
  include 'DIFFSIZES.inc'

  external :: ccopy
  external :: ccopy_dv

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, idir  ! Loop counters
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  integer :: nsize
  complex(4), dimension(4) :: cx
  integer :: incx_val
  complex(4), dimension(max_size) :: cy
  integer :: incy_val

  ! Vector mode derivative variables (type-promoted)
  ! Scalars become arrays(nbdirsmax), arrays gain extra dimension
  complex(4), dimension(nbdirsmax,4) :: cx_dv
  complex(4), dimension(nbdirsmax,max_size) :: cy_dv
  ! Declare variables for storing original values
  complex(4), dimension(4) :: cx_orig
  complex(4), dimension(nbdirsmax,4) :: cx_dv_orig
  complex(4), dimension(max_size) :: cy_orig
  complex(4), dimension(nbdirsmax,max_size) :: cy_dv_orig

  ! Initialize test parameters
  nsize = n
  incx_val = 1
  incy_val = 1

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  seed_array = 42
  call random_seed(put=seed_array)

  do i = 1, max_size
    call random_number(temp_real)
    call random_number(temp_imag)
    cx(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  end do
  do i = 1, max_size
    call random_number(temp_real)
    call random_number(temp_imag)
    cy(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  end do

  ! Initialize input derivatives to random values (exactly like scalar mode)
  do idir = 1, nbdirsmax
    do i = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      cx_dv(idir,i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
  end do
  do idir = 1, nbdirsmax
    do i = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      cy_dv(idir,i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
  end do

  write(*,*) 'Testing CCOPY (Vector Forward Mode)'
  ! Store original values before any function calls (critical for INOUT parameters)
  cx_orig = cx
  cx_dv_orig = cx_dv
  cy_orig = cy
  cy_dv_orig = cy_dv

  ! Call the vector mode differentiated function

  ! Set ISIZE globals required by differentiated routine
  call set_ISIZE1OFCy(max_size)

  call ccopy_dv(nsize, cx, cx_dv, incx_val, cy, cy_dv, incy_val, nbdirsmax)

  ! Reset ISIZE globals to uninitialized (-1)
  call set_ISIZE1OFCy(-1)

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
    complex(4) :: central_diff, ad_result
    integer :: i, j, idir
    logical :: has_large_errors
    complex(4), dimension(max_size) :: cy_forward, cy_backward
    
    max_error = 0.0e0
    has_large_errors = .false.
    
    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Number of directions:', nbdirsmax
    
    ! Test each derivative direction separately
    do idir = 1, nbdirsmax
      
      ! Forward perturbation: f(x + h * direction)
      cx = cx_orig + cmplx(h, 0.0) * cx_dv_orig(idir,:)
      cy = cy_orig + cmplx(h, 0.0) * cy_dv_orig(idir,:)
      call ccopy(nsize, cx, incx_val, cy, incy_val)
      cy_forward = cy
      
      ! Backward perturbation: f(x - h * direction)
      cx = cx_orig - cmplx(h, 0.0) * cx_dv_orig(idir,:)
      cy = cy_orig - cmplx(h, 0.0) * cy_dv_orig(idir,:)
      call ccopy(nsize, cx, incx_val, cy, incy_val)
      cy_backward = cy
      
      ! Compute central differences and compare with AD results
      do i = 1, min(2, nsize)  ! Check only first few elements
        ! Central difference: (f(x+h) - f(x-h)) / (2h)
        central_diff = (cy_forward(i) - cy_backward(i)) / (2.0e0 * h)
        ! AD result
        ad_result = cy_dv(idir,i)
        ! Error check: |a - b| > atol + rtol * |b|
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 1.0e-3 + 1.0e-3 * abs_reference
        if (abs_error > error_bound) then
          has_large_errors = .true.
          relative_error = abs_error / max(abs_reference, 1.0e-10)
          write(*,*) '  Large error in direction', idir, ' output CY(', i, '):'
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
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in vector derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Vector derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_ccopy_vector_forward