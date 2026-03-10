! Test program for CDOTU vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=4

program test_cdotu_vector_forward
  implicit none
  integer, parameter :: nbdirs = 4

  complex(4), external :: cdotu
  external :: cdotu_dv

  ! Test parameters
  integer :: n  ! Current size (set in loop)
  integer, parameter :: max_size = 100  ! Maximum array dimension (multi-size: 1,4,40,100)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, idir  ! Loop counters
  integer :: test_sizes(1), itest
  logical :: passed, all_passed
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  integer :: nsize
  complex(4), dimension(max_size) :: cx
  integer :: incx_val
  complex(4), dimension(max_size) :: cy
  integer :: incy_val

  ! Vector mode derivative variables (type-promoted)
  ! Scalars become arrays(nbdirs), arrays gain extra dimension
  complex(4), dimension(nbdirs,max_size) :: cx_dv
  complex(4), dimension(nbdirs,max_size) :: cy_dv
  ! Declare variables for storing original values
  complex(4), dimension(max_size) :: cx_orig
  complex(4), dimension(nbdirs,max_size) :: cx_dv_orig
  complex(4), dimension(max_size) :: cy_orig
  complex(4), dimension(nbdirs,max_size) :: cy_dv_orig

  ! Function result variables
  complex(4) :: cdotu_result
  complex(4), dimension(nbdirs) :: cdotu_dv_result

  test_sizes = (/ 4 /)
  write(*,*) 'Testing CDOTU (Vector Forward, multi-size: n = 4)'
  all_passed = .true.
  do itest = 1, 1
    n = test_sizes(itest)
    write(*,*) 'Testing CDOTU (Vector Forward, n =', n, ')'

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
  do idir = 1, nbdirs
    do i = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      cx_dv(idir,i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
  end do
  do idir = 1, nbdirs
    do i = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      cy_dv(idir,i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
  end do

  write(*,*) 'Testing CDOTU (Vector Forward Mode)'
  ! Store original values before any function calls (critical for INOUT parameters)
  cx_orig = cx
  cx_dv_orig = cx_dv
  cy_orig = cy
  cy_dv_orig = cy_dv

  ! Call the vector mode differentiated function

  call cdotu_dv(nsize, cx, cx_dv, incx_val, cy, cy_dv, incy_val, cdotu_result, cdotu_dv_result, nbdirs)

  ! Print results and compare
  write(*,*) 'Function calls completed successfully'

  ! Numerical differentiation check
  call check_derivatives_numerically(passed)
  all_passed = all_passed .and. passed
  end do
  if (all_passed) then
    write(*,*) 'PASS: Vector forward mode - all sizes completed successfully'
  else
    write(*,*) 'FAIL: Vector forward mode - one or more sizes had derivative errors'
  end if

contains

  subroutine check_derivatives_numerically(passed)
    implicit none
    logical, intent(out) :: passed
    real(4), parameter :: h = 1.0e-3  ! Step size for finite differences
    real(4) :: relative_error, max_error
    real(4) :: abs_error, abs_reference, error_bound
    complex(4) :: central_diff, ad_result
    integer :: i, j, idir
    logical :: has_large_errors
    complex(4) :: cdotu_forward, cdotu_backward
    
    max_error = 0.0e0
    has_large_errors = .false.
    
    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Number of directions:', nbdirs
    
    ! Test each derivative direction separately
    do idir = 1, nbdirs
      
      ! Forward perturbation: f(x + h * direction)
      cx = cx_orig + cmplx(h, 0.0) * cx_dv_orig(idir,:)
      cy = cy_orig + cmplx(h, 0.0) * cy_dv_orig(idir,:)
      cdotu_forward = cdotu(nsize, cx, incx_val, cy, incy_val)
      
      ! Backward perturbation: f(x - h * direction)
      cx = cx_orig - cmplx(h, 0.0) * cx_dv_orig(idir,:)
      cy = cy_orig - cmplx(h, 0.0) * cy_dv_orig(idir,:)
      cdotu_backward = cdotu(nsize, cx, incx_val, cy, incy_val)
      
      ! Compute central differences and compare with AD results
      ! Central difference: (f(x+h) - f(x-h)) / (2h)
      central_diff = (cdotu_forward - cdotu_backward) / (2.0e0 * h)
      ! AD result
      ad_result = cdotu_dv_result(idir)
      ! Error check: |a - b| > atol + rtol * |b|
      abs_error = abs(central_diff - ad_result)
      abs_reference = abs(ad_result)
      error_bound = 1.0e-3 + 1.0e-3 * abs_reference
      if (abs_error > error_bound) then
        has_large_errors = .true.
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        write(*,*) '  Large error in direction', idir, ' output CDOTU:'
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
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in vector derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Vector derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_cdotu_vector_forward