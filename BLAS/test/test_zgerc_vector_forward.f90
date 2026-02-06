! Test program for ZGERC vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirsmax=4

program test_zgerc_vector_forward
  implicit none
  include 'DIFFSIZES.inc'

  external :: zgerc
  external :: zgerc_dv

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, idir  ! Loop counters
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  integer :: msize
  integer :: nsize
  complex(8) :: alpha
  complex(8), dimension(max_size) :: x
  integer :: incx_val
  complex(8), dimension(max_size) :: y
  integer :: incy_val
  complex(8), dimension(max_size,max_size) :: a
  integer :: lda_val

  ! Vector mode derivative variables (type-promoted)
  ! Scalars become arrays(nbdirsmax), arrays gain extra dimension
  complex(8), dimension(nbdirsmax) :: alpha_dv
  complex(8), dimension(nbdirsmax,max_size) :: x_dv
  complex(8), dimension(nbdirsmax,max_size) :: y_dv
  complex(8), dimension(nbdirsmax,max_size,max_size) :: a_dv
  ! Declare variables for storing original values
  complex(8) :: alpha_orig
  complex(8), dimension(nbdirsmax) :: alpha_dv_orig
  complex(8), dimension(max_size) :: x_orig
  complex(8), dimension(nbdirsmax,max_size) :: x_dv_orig
  complex(8), dimension(max_size) :: y_orig
  complex(8), dimension(nbdirsmax,max_size) :: y_dv_orig
  complex(8), dimension(max_size,max_size) :: a_orig
  complex(8), dimension(nbdirsmax,max_size,max_size) :: a_dv_orig

  ! Initialize test parameters
  msize = n
  nsize = n
  incx_val = 1
  incy_val = 1
  lda_val = lda

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  seed_array = 42
  call random_seed(put=seed_array)

  call random_number(temp_real)
  call random_number(temp_imag)
  alpha = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  do i = 1, max_size
    call random_number(temp_real)
    call random_number(temp_imag)
    x(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  end do
  do i = 1, max_size
    call random_number(temp_real)
    call random_number(temp_imag)
    y(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  end do
  do i = 1, max_size
    do j = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      a(i,j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
  end do

  ! Initialize input derivatives to random values (exactly like scalar mode)
  do idir = 1, nbdirsmax
    call random_number(temp_real)
    call random_number(temp_imag)
    alpha_dv(idir) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  end do
  do idir = 1, nbdirsmax
    do i = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      x_dv(idir,i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
  end do
  do idir = 1, nbdirsmax
    do i = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      y_dv(idir,i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
  end do
  do idir = 1, nbdirsmax
    do i = 1, max_size
      do j = 1, max_size
        call random_number(temp_real)
        call random_number(temp_imag)
        a_dv(idir,i,j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
      end do
    end do
  end do

  write(*,*) 'Testing ZGERC (Vector Forward Mode)'
  ! Store original values before any function calls (critical for INOUT parameters)
  alpha_orig = alpha
  alpha_dv_orig = alpha_dv
  x_orig = x
  x_dv_orig = x_dv
  y_orig = y
  y_dv_orig = y_dv
  a_orig = a
  a_dv_orig = a_dv

  ! Call the vector mode differentiated function

  call zgerc_dv(msize, nsize, alpha, alpha_dv, x, x_dv, incx_val, y, y_dv, incy_val, a, a_dv, lda_val, nbdirsmax)

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
    complex(8), dimension(max_size,max_size) :: a_forward, a_backward
    
    max_error = 0.0e0
    has_large_errors = .false.
    
    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Number of directions:', nbdirsmax
    
    ! Test each derivative direction separately
    do idir = 1, nbdirsmax
      
      ! Forward perturbation: f(x + h * direction)
      alpha = alpha_orig + cmplx(h, 0.0) * alpha_dv_orig(idir)
      x = x_orig + cmplx(h, 0.0) * x_dv_orig(idir,:)
      y = y_orig + cmplx(h, 0.0) * y_dv_orig(idir,:)
      a = a_orig + cmplx(h, 0.0) * a_dv_orig(idir,:,:)
      call zgerc(msize, nsize, alpha, x, incx_val, y, incy_val, a, lda_val)
      a_forward = a
      
      ! Backward perturbation: f(x - h * direction)
      alpha = alpha_orig - cmplx(h, 0.0) * alpha_dv_orig(idir)
      x = x_orig - cmplx(h, 0.0) * x_dv_orig(idir,:)
      y = y_orig - cmplx(h, 0.0) * y_dv_orig(idir,:)
      a = a_orig - cmplx(h, 0.0) * a_dv_orig(idir,:,:)
      call zgerc(msize, nsize, alpha, x, incx_val, y, incy_val, a, lda_val)
      a_backward = a
      
      ! Compute central differences and compare with AD results
      do j = 1, min(2, nsize)  ! Check only first few elements
        do i = 1, min(2, nsize)
          ! Central difference: (f(x+h) - f(x-h)) / (2h)
          central_diff = (a_forward(i,j) - a_backward(i,j)) / (2.0e0 * h)
          ! AD result
          ad_result = a_dv(idir,i,j)
          ! Error check: |a - b| > atol + rtol * |b|
          abs_error = abs(central_diff - ad_result)
          abs_reference = abs(ad_result)
          error_bound = 1.0e-5 + 1.0e-5 * abs_reference
          if (abs_error > error_bound) then
            has_large_errors = .true.
            relative_error = abs_error / max(abs_reference, 1.0e-10)
            write(*,*) '  Large error in direction', idir, ' output A(', i, ',', j, '):'
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
    end do
    
    write(*,*) 'Maximum relative error across all directions:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in vector derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Vector derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_zgerc_vector_forward