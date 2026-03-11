! Test program for STRSV vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=4

program test_strsv_vector_forward
  implicit none
  integer, parameter :: nbdirs = 4

  external :: strsv
  external :: strsv_dv

  ! Test parameters
  integer :: n  ! Current size (set in loop)
  integer, parameter :: max_size = 100  ! Maximum array dimension (multi-size: 1,4,40,100)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, idir  ! Loop counters
  integer :: test_sizes(1), itest
  logical :: passed, all_passed
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  character :: uplo
  character :: trans
  character :: diag
  integer :: nsize
  real(4), dimension(max_size,max_size) :: a
  integer :: lda_val
  real(4), dimension(max_size) :: x
  integer :: incx_val

  ! Vector mode derivative variables (type-promoted)
  ! Scalars become arrays(nbdirs), arrays gain extra dimension
  real(4), dimension(nbdirs,max_size,max_size) :: a_dv
  real(4), dimension(nbdirs,max_size) :: x_dv
  ! Declare variables for storing original values
  real(4), dimension(max_size,max_size) :: a_orig
  real(4), dimension(nbdirs,max_size,max_size) :: a_dv_orig
  real(4), dimension(max_size) :: x_orig
  real(4), dimension(nbdirs,max_size) :: x_dv_orig

  test_sizes = (/ 4 /)
  write(*,*) 'Testing STRSV (Vector Forward, multi-size: n = 4)'
  all_passed = .true.
  do itest = 1, 1
    n = test_sizes(itest)
    write(*,*) 'Testing STRSV (Vector Forward, n =', n, ')'

    call run_test_for_size(n, passed)
  all_passed = all_passed .and. passed
  end do
  if (all_passed) then
    write(*,*) 'PASS: Vector forward mode - all sizes completed successfully'
  else
    write(*,*) 'FAIL: Vector forward mode - one or more sizes had derivative errors'
  end if

contains

  subroutine run_test_for_size(n, passed)
    implicit none
    integer, intent(in) :: n
    logical, intent(out) :: passed

    ! Initialize test parameters
    nsize = n
    lda_val = lda
    incx_val = 1
    
    ! Initialize test data with random numbers
    ! Initialize random seed for reproducible results
    seed_array = 42
    call random_seed(put=seed_array)
    
    uplo = 'U'
    trans = 'N'
    diag = 'N'
    call random_number(a)
    a = a * 2.0 - 1.0  ! Scale to [-1,1]
    call random_number(x)
    x = x * 2.0 - 1.0  ! Scale to [-1,1]
    
    ! Initialize input derivatives to random values (exactly like scalar mode)
    do idir = 1, nbdirs
      call random_number(a_dv(idir,:,:))
      a_dv(idir,:,:) = a_dv(idir,:,:) * 2.0 - 1.0
    end do
    do idir = 1, nbdirs
      call random_number(x_dv(idir,:))
      x_dv(idir,:) = x_dv(idir,:) * 2.0 - 1.0
    end do
    
    write(*,*) 'Testing STRSV (Vector Forward Mode)'
    ! Store original values before any function calls (critical for INOUT parameters)
    a_orig = a
    a_dv_orig = a_dv
    x_orig = x
    x_dv_orig = x_dv
    
    ! Call the vector mode differentiated function
    
    call strsv_dv(uplo, trans, diag, nsize, a, a_dv, lda_val, x, x_dv, incx_val, nbdirs)
    
    ! Print results and compare
    write(*,*) 'Function calls completed successfully'
    
    ! Numerical differentiation check
    call check_derivatives_numerically(passed)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(passed)
    implicit none
    logical, intent(out) :: passed
    real(4), parameter :: h = 1.0e-3  ! Step size for finite differences
    real(4) :: relative_error, max_error
    real(4) :: abs_error, abs_reference, error_bound
    real(4) :: central_diff, ad_result
    integer :: i, j, idir
    logical :: has_large_errors
    real(4), dimension(max_size) :: x_forward, x_backward
    
    max_error = 0.0e0
    has_large_errors = .false.
    
    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Number of directions:', nbdirs
    
    ! Test each derivative direction separately
    do idir = 1, nbdirs
      
      ! Forward perturbation: f(x + h * direction)
      a = a_orig + h * a_dv_orig(idir,:,:)
      x = x_orig + h * x_dv_orig(idir,:)
      call strsv(uplo, trans, diag, nsize, a, lda_val, x, incx_val)
      x_forward = x
      
      ! Backward perturbation: f(x - h * direction)
      a = a_orig - h * a_dv_orig(idir,:,:)
      x = x_orig - h * x_dv_orig(idir,:)
      call strsv(uplo, trans, diag, nsize, a, lda_val, x, incx_val)
      x_backward = x
      
      ! Compute central differences and compare with AD results
      do i = 1, min(2, nsize)  ! Check only first few elements
        ! Central difference: (f(x+h) - f(x-h)) / (2h)
        central_diff = (x_forward(i) - x_backward(i)) / (2.0e0 * h)
        ! AD result
        ad_result = x_dv(idir,i)
        ! Error check: |a - b| > atol + rtol * |b|
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 2.0e-3 + 2.0e-3 * abs_reference
        if (abs_error > error_bound) then
          has_large_errors = .true.
          relative_error = abs_error / max(abs_reference, 1.0e-10)
          write(*,*) '  Large error in direction', idir, ' output X(', i, '):'
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
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in vector derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Vector derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_strsv_vector_forward