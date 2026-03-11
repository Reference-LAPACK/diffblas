! Test program for CTRSV vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=4

program test_ctrsv_vector_forward
  implicit none
  integer, parameter :: nbdirs = 4

  external :: ctrsv
  external :: ctrsv_dv

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
  complex(4), dimension(max_size,max_size) :: a
  integer :: lda_val
  complex(4), dimension(max_size) :: x
  integer :: incx_val

  ! Vector mode derivative variables (type-promoted)
  ! Scalars become arrays(nbdirs), arrays gain extra dimension
  complex(4), dimension(nbdirs,max_size,max_size) :: a_dv
  complex(4), dimension(nbdirs,max_size) :: x_dv
  ! Declare variables for storing original values
  complex(4), dimension(max_size,max_size) :: a_orig
  complex(4), dimension(nbdirs,max_size,max_size) :: a_dv_orig
  complex(4), dimension(max_size) :: x_orig
  complex(4), dimension(nbdirs,max_size) :: x_dv_orig

  test_sizes = (/ 4 /)
  write(*,*) 'Testing CTRSV (Vector Forward, multi-size: n = 4)'
  all_passed = .true.
  do itest = 1, 1
    n = test_sizes(itest)
    write(*,*) 'Testing CTRSV (Vector Forward, n =', n, ')'

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
    do i = 1, max_size
      do j = 1, max_size
        call random_number(temp_real)
        call random_number(temp_imag)
        a(i,j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
      end do
    end do
    do i = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      x(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
    
    ! Initialize input derivatives to random values (exactly like scalar mode)
    do idir = 1, nbdirs
      do i = 1, max_size
        do j = 1, max_size
          call random_number(temp_real)
          call random_number(temp_imag)
          a_dv(idir,i,j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
        end do
      end do
    end do
    do idir = 1, nbdirs
      do i = 1, max_size
        call random_number(temp_real)
        call random_number(temp_imag)
        x_dv(idir,i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
      end do
    end do
    
    write(*,*) 'Testing CTRSV (Vector Forward Mode)'
    ! Store original values before any function calls (critical for INOUT parameters)
    a_orig = a
    a_dv_orig = a_dv
    x_orig = x
    x_dv_orig = x_dv
    
    ! Call the vector mode differentiated function
    
    call ctrsv_dv(uplo, trans, diag, nsize, a, a_dv, lda_val, x, x_dv, incx_val, nbdirs)
    
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
    complex(4) :: central_diff, ad_result
    integer :: i, j, idir
    logical :: has_large_errors
    complex(4), dimension(max_size) :: x_forward, x_backward
    
    max_error = 0.0e0
    has_large_errors = .false.
    
    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Number of directions:', nbdirs
    
    ! Test each derivative direction separately
    do idir = 1, nbdirs
      
      ! Forward perturbation: f(x + h * direction)
      a = a_orig + cmplx(h, 0.0) * a_dv_orig(idir,:,:)
      x = x_orig + cmplx(h, 0.0) * x_dv_orig(idir,:)
      call ctrsv(uplo, trans, diag, nsize, a, lda_val, x, incx_val)
      x_forward = x
      
      ! Backward perturbation: f(x - h * direction)
      a = a_orig - cmplx(h, 0.0) * a_dv_orig(idir,:,:)
      x = x_orig - cmplx(h, 0.0) * x_dv_orig(idir,:)
      call ctrsv(uplo, trans, diag, nsize, a, lda_val, x, incx_val)
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
        error_bound = 1.0e-3 + 1.0e-3 * abs_reference
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
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in vector derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Vector derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_ctrsv_vector_forward