! Test program for SGBMV vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=4

program test_sgbmv_vector_forward
  implicit none
  integer, parameter :: nbdirs = 4

  external :: sgbmv
  external :: sgbmv_dv

  ! Test parameters
  integer :: n  ! Current size (set in loop)
  integer, parameter :: max_size = 100  ! Maximum array dimension (multi-size: 1,4,40,100)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, idir, band_row  ! Loop counters
  integer :: test_sizes(1), itest
  logical :: passed, all_passed
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  character :: trans
  integer :: msize
  integer :: nsize
  integer :: kl
  integer :: ku
  real(4) :: alpha
  real(4), dimension(max_size,max_size) :: a
  integer :: lda_val
  real(4), dimension(max_size) :: x
  integer :: incx_val
  real(4) :: beta
  real(4), dimension(max_size) :: y
  integer :: incy_val

  ! Vector mode derivative variables (type-promoted)
  ! Scalars become arrays(nbdirs), arrays gain extra dimension
  real(4), dimension(nbdirs) :: alpha_dv
  real(4), dimension(nbdirs,max_size,max_size) :: a_dv
  real(4), dimension(nbdirs,max_size) :: x_dv
  real(4), dimension(nbdirs) :: beta_dv
  real(4), dimension(nbdirs,max_size) :: y_dv
  ! Declare variables for storing original values
  real(4) :: alpha_orig
  real(4), dimension(nbdirs) :: alpha_dv_orig
  real(4), dimension(max_size,max_size) :: a_orig
  real(4), dimension(nbdirs,max_size,max_size) :: a_dv_orig
  real(4), dimension(max_size) :: x_orig
  real(4), dimension(nbdirs,max_size) :: x_dv_orig
  real(4) :: beta_orig
  real(4), dimension(nbdirs) :: beta_dv_orig
  real(4), dimension(max_size) :: y_orig
  real(4), dimension(nbdirs,max_size) :: y_dv_orig

  test_sizes = (/ 4 /)
  write(*,*) 'Testing SGBMV (Vector Forward, multi-size: n = 4)'
  all_passed = .true.
  do itest = 1, 1
    n = test_sizes(itest)
    write(*,*) 'Testing SGBMV (Vector Forward, n =', n, ')'

  ! Initialize test parameters
  msize = n
  nsize = n
  kl = 1
  ku = 1
  lda_val = lda
  incx_val = 1
  incy_val = 1

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  seed_array = 42
  call random_seed(put=seed_array)

  trans = 'N'
  call random_number(alpha)
  alpha = alpha * 2.0 - 1.0  ! Scale to [-1,1]
  ! Initialize a as general band matrix (kl, ku band storage)
  do j = 1, n
    do band_row = max(1, ku+2-j), min(kl+ku+1, ku+msize-j+1)
      call random_number(temp_real)
      a(band_row, j) = temp_real * 2.0 - 1.0
    end do
  end do
  call random_number(x)
  x = x * 2.0 - 1.0  ! Scale to [-1,1]
  call random_number(beta)
  beta = beta * 2.0 - 1.0  ! Scale to [-1,1]
  call random_number(y)
  y = y * 2.0 - 1.0  ! Scale to [-1,1]

  ! Initialize input derivatives to random values (exactly like scalar mode)
  do idir = 1, nbdirs
    call random_number(temp_real)
    alpha_dv(idir) = temp_real * 2.0 - 1.0
  end do
  do idir = 1, nbdirs
    call random_number(a_dv(idir,:,:))
    a_dv(idir,:,:) = a_dv(idir,:,:) * 2.0 - 1.0
  end do
  do idir = 1, nbdirs
    call random_number(x_dv(idir,:))
    x_dv(idir,:) = x_dv(idir,:) * 2.0 - 1.0
  end do
  do idir = 1, nbdirs
    call random_number(temp_real)
    beta_dv(idir) = temp_real * 2.0 - 1.0
  end do
  do idir = 1, nbdirs
    call random_number(y_dv(idir,:))
    y_dv(idir,:) = y_dv(idir,:) * 2.0 - 1.0
  end do

  write(*,*) 'Testing SGBMV (Vector Forward Mode)'
  ! Store original values before any function calls (critical for INOUT parameters)
  alpha_orig = alpha
  alpha_dv_orig = alpha_dv
  a_orig = a
  a_dv_orig = a_dv
  x_orig = x
  x_dv_orig = x_dv
  beta_orig = beta
  beta_dv_orig = beta_dv
  y_orig = y
  y_dv_orig = y_dv

  ! Call the vector mode differentiated function

  call sgbmv_dv(trans, msize, nsize, kl, ku, alpha, alpha_dv, a, a_dv, lda_val, x, x_dv, incx_val, beta, beta_dv, y, y_dv, incy_val, nbdirs)

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
    real(4) :: central_diff, ad_result
    integer :: i, j, idir, band_row
    logical :: has_large_errors
    real(4), dimension(max_size) :: y_forward, y_backward
    
    max_error = 0.0e0
    has_large_errors = .false.
    
    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Number of directions:', nbdirs
    
    ! Test each derivative direction separately
    do idir = 1, nbdirs
      
      ! Forward perturbation: f(x + h * direction)
      alpha = alpha_orig + h * alpha_dv_orig(idir)
      a = a_orig + h * a_dv_orig(idir,:,:)
      x = x_orig + h * x_dv_orig(idir,:)
      beta = beta_orig + h * beta_dv_orig(idir)
      y = y_orig + h * y_dv_orig(idir,:)
      call sgbmv(trans, msize, nsize, kl, ku, alpha, a, lda_val, x, incx_val, beta, y, incy_val)
      y_forward = y
      
      ! Backward perturbation: f(x - h * direction)
      alpha = alpha_orig - h * alpha_dv_orig(idir)
      a = a_orig - h * a_dv_orig(idir,:,:)
      x = x_orig - h * x_dv_orig(idir,:)
      beta = beta_orig - h * beta_dv_orig(idir)
      y = y_orig - h * y_dv_orig(idir,:)
      call sgbmv(trans, msize, nsize, kl, ku, alpha, a, lda_val, x, incx_val, beta, y, incy_val)
      y_backward = y
      
      ! Compute central differences and compare with AD results
      do i = 1, min(2, nsize)  ! Check only first few elements
        ! Central difference: (f(x+h) - f(x-h)) / (2h)
        central_diff = (y_forward(i) - y_backward(i)) / (2.0e0 * h)
        ! AD result
        ad_result = y_dv(idir,i)
        ! Error check: |a - b| > atol + rtol * |b|
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 2.0e-3 + 2.0e-3 * abs_reference
        if (abs_error > error_bound) then
          has_large_errors = .true.
          relative_error = abs_error / max(abs_reference, 1.0e-10)
          write(*,*) '  Large error in direction', idir, ' output Y(', i, '):'
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

end program test_sgbmv_vector_forward