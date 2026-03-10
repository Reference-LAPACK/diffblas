! Test program for ZDOTC vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirs=4

program test_zdotc_vector_forward
  implicit none
  integer, parameter :: nbdirs = 4

  complex(8), external :: zdotc
  external :: zdotc_dv

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
  complex(8), dimension(max_size) :: zx
  integer :: incx_val
  complex(8), dimension(max_size) :: zy
  integer :: incy_val

  ! Vector mode derivative variables (type-promoted)
  ! Scalars become arrays(nbdirs), arrays gain extra dimension
  complex(8), dimension(nbdirs,max_size) :: zx_dv
  complex(8), dimension(nbdirs,max_size) :: zy_dv
  ! Declare variables for storing original values
  complex(8), dimension(max_size) :: zx_orig
  complex(8), dimension(nbdirs,max_size) :: zx_dv_orig
  complex(8), dimension(max_size) :: zy_orig
  complex(8), dimension(nbdirs,max_size) :: zy_dv_orig

  ! Function result variables
  complex(8) :: zdotc_result
  complex(8), dimension(nbdirs) :: zdotc_dv_result

  test_sizes = (/ 4 /)
  write(*,*) 'Testing ZDOTC (Vector Forward, multi-size: n = 4)'
  all_passed = .true.
  do itest = 1, 1
    n = test_sizes(itest)
    write(*,*) 'Testing ZDOTC (Vector Forward, n =', n, ')'

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
  do idir = 1, nbdirs
    do i = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      zx_dv(idir,i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
  end do
  do idir = 1, nbdirs
    do i = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      zy_dv(idir,i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
  end do

  write(*,*) 'Testing ZDOTC (Vector Forward Mode)'
  ! Store original values before any function calls (critical for INOUT parameters)
  zx_orig = zx
  zx_dv_orig = zx_dv
  zy_orig = zy
  zy_dv_orig = zy_dv

  ! Call the vector mode differentiated function

  call zdotc_dv(nsize, zx, zx_dv, incx_val, zy, zy_dv, incy_val, zdotc_result, zdotc_dv_result, nbdirs)

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
    real(8), parameter :: h = 1.0e-7  ! Step size for finite differences
    real(8) :: relative_error, max_error
    real(8) :: abs_error, abs_reference, error_bound
    complex(8) :: central_diff, ad_result
    integer :: i, j, idir
    logical :: has_large_errors
    complex(8) :: zdotc_forward, zdotc_backward
    
    max_error = 0.0e0
    has_large_errors = .false.
    
    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Number of directions:', nbdirs
    
    ! Test each derivative direction separately
    do idir = 1, nbdirs
      
      ! Forward perturbation: f(x + h * direction)
      zx = zx_orig + cmplx(h, 0.0) * zx_dv_orig(idir,:)
      zy = zy_orig + cmplx(h, 0.0) * zy_dv_orig(idir,:)
      zdotc_forward = zdotc(nsize, zx, incx_val, zy, incy_val)
      
      ! Backward perturbation: f(x - h * direction)
      zx = zx_orig - cmplx(h, 0.0) * zx_dv_orig(idir,:)
      zy = zy_orig - cmplx(h, 0.0) * zy_dv_orig(idir,:)
      zdotc_backward = zdotc(nsize, zx, incx_val, zy, incy_val)
      
      ! Compute central differences and compare with AD results
      ! Central difference: (f(x+h) - f(x-h)) / (2h)
      central_diff = (zdotc_forward - zdotc_backward) / (2.0e0 * h)
      ! AD result
      ad_result = zdotc_dv_result(idir)
      ! Error check: |a - b| > atol + rtol * |b|
      abs_error = abs(central_diff - ad_result)
      abs_reference = abs(ad_result)
      error_bound = 1.0e-5 + 1.0e-5 * abs_reference
      if (abs_error > error_bound) then
        has_large_errors = .true.
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        write(*,*) '  Large error in direction', idir, ' output ZDOTC:'
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
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in vector derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Vector derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_zdotc_vector_forward