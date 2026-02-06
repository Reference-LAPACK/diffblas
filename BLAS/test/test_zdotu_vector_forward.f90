! Test program for ZDOTU vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirsmax=4

program test_zdotu_vector_forward
  implicit none
  include 'DIFFSIZES.inc'

  complex(8), external :: zdotu
  external :: zdotu_dv

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, idir  ! Loop counters
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  integer :: nsize
  complex(8), dimension(4) :: zx
  integer :: incx_val
  complex(8), dimension(4) :: zy
  integer :: incy_val

  ! Vector mode derivative variables (type-promoted)
  ! Scalars become arrays(nbdirsmax), arrays gain extra dimension
  complex(8), dimension(nbdirsmax,4) :: zx_dv
  complex(8), dimension(nbdirsmax,4) :: zy_dv
  ! Declare variables for storing original values
  complex(8), dimension(4) :: zx_orig
  complex(8), dimension(nbdirsmax,4) :: zx_dv_orig
  complex(8), dimension(4) :: zy_orig
  complex(8), dimension(nbdirsmax,4) :: zy_dv_orig

  ! Function result variables
  complex(8) :: zdotu_result
  complex(8), dimension(nbdirsmax) :: zdotu_dv_result

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
    zx(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  end do
  do i = 1, max_size
    call random_number(temp_real)
    call random_number(temp_imag)
    zy(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  end do

  ! Initialize input derivatives to random values (exactly like scalar mode)
  do idir = 1, nbdirsmax
    do i = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      zx_dv(idir,i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
  end do
  do idir = 1, nbdirsmax
    do i = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      zy_dv(idir,i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
  end do

  write(*,*) 'Testing ZDOTU (Vector Forward Mode)'
  ! Store original values before any function calls (critical for INOUT parameters)
  zx_orig = zx
  zx_dv_orig = zx_dv
  zy_orig = zy
  zy_dv_orig = zy_dv

  ! Call the vector mode differentiated function

  call zdotu_dv(nsize, zx, zx_dv, incx_val, zy, zy_dv, incy_val, zdotu_result, zdotu_dv_result, nbdirsmax)

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
    complex(8) :: central_diff, ad_result
    integer :: i, j, idir
    logical :: has_large_errors
    complex(8) :: zdotu_forward, zdotu_backward
    
    max_error = 0.0e0
    has_large_errors = .false.
    
    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Number of directions:', nbdirsmax
    
    ! Test each derivative direction separately
    do idir = 1, nbdirsmax
      
      ! Forward perturbation: f(x + h * direction)
      zx = zx_orig + cmplx(h, 0.0) * zx_dv_orig(idir,:)
      zy = zy_orig + cmplx(h, 0.0) * zy_dv_orig(idir,:)
      zdotu_forward = zdotu(nsize, zx, incx_val, zy, incy_val)
      
      ! Backward perturbation: f(x - h * direction)
      zx = zx_orig - cmplx(h, 0.0) * zx_dv_orig(idir,:)
      zy = zy_orig - cmplx(h, 0.0) * zy_dv_orig(idir,:)
      zdotu_backward = zdotu(nsize, zx, incx_val, zy, incy_val)
      
      ! Compute central differences and compare with AD results
      ! Central difference: (f(x+h) - f(x-h)) / (2h)
      central_diff = (zdotu_forward - zdotu_backward) / (2.0e0 * h)
      ! AD result
      ad_result = zdotu_dv_result(idir)
      ! Error check: |a - b| > atol + rtol * |b|
      abs_error = abs(central_diff - ad_result)
      abs_reference = abs(ad_result)
      error_bound = 1.0e-5 + 1.0e-5 * abs_reference
      if (abs_error > error_bound) then
        has_large_errors = .true.
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        write(*,*) '  Large error in direction', idir, ' output ZDOTU:'
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
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in vector derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Vector derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_zdotu_vector_forward