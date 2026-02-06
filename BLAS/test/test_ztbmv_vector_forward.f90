! Test program for ZTBMV vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirsmax=4

program test_ztbmv_vector_forward
  implicit none
  include 'DIFFSIZES.inc'

  external :: ztbmv
  external :: ztbmv_dv

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, idir, band_row  ! Loop counters
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  character :: uplo
  character :: trans
  character :: diag
  integer :: nsize
  integer :: ksize
  complex(8), dimension(max_size,max_size) :: a
  integer :: lda_val
  complex(8), dimension(max_size) :: x
  integer :: incx_val

  ! Vector mode derivative variables (type-promoted)
  ! Scalars become arrays(nbdirsmax), arrays gain extra dimension
  complex(8), dimension(nbdirsmax,max_size,max_size) :: a_dv
  complex(8), dimension(nbdirsmax,max_size) :: x_dv
  ! Declare variables for storing original values
  complex(8), dimension(max_size,max_size) :: a_orig
  complex(8), dimension(nbdirsmax,max_size,max_size) :: a_dv_orig
  complex(8), dimension(max_size) :: x_orig
  complex(8), dimension(nbdirsmax,max_size) :: x_dv_orig

  ! Initialize test parameters
  nsize = n
  ksize = max(0, n - 1)  ! Band width: 0 <= K <= N-1
  lda_val = lda
  incx_val = 1

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  seed_array = 42
  call random_seed(put=seed_array)

  uplo = 'U'
  trans = 'N'
  diag = 'N'
  ! Initialize a as triangular band matrix (upper band storage)
  do j = 1, n
    do band_row = max(1, ksize+2-j), ksize+1
      call random_number(temp_real)
      call random_number(temp_imag)
      a(band_row, j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
  end do
  do i = 1, max_size
    call random_number(temp_real)
    call random_number(temp_imag)
    x(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  end do

  ! Initialize input derivatives to random values (exactly like scalar mode)
  do idir = 1, nbdirsmax
    do i = 1, max_size
      do j = 1, max_size
        call random_number(temp_real)
        call random_number(temp_imag)
        a_dv(idir,i,j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
      end do
    end do
  end do
  do idir = 1, nbdirsmax
    do i = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      x_dv(idir,i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
  end do

  write(*,*) 'Testing ZTBMV (Vector Forward Mode)'
  ! Store original values before any function calls (critical for INOUT parameters)
  a_orig = a
  a_dv_orig = a_dv
  x_orig = x
  x_dv_orig = x_dv

  ! Call the vector mode differentiated function

  call ztbmv_dv(uplo, trans, diag, nsize, ksize, a, a_dv, lda_val, x, x_dv, incx_val, nbdirsmax)

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
    integer :: i, j, idir, band_row
    logical :: has_large_errors
    complex(8), dimension(max_size) :: x_forward, x_backward
    
    max_error = 0.0e0
    has_large_errors = .false.
    
    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Number of directions:', nbdirsmax
    
    ! Test each derivative direction separately
    do idir = 1, nbdirsmax
      
      ! Forward perturbation: f(x + h * direction)
      a = a_orig + cmplx(h, 0.0) * a_dv_orig(idir,:,:)
      x = x_orig + cmplx(h, 0.0) * x_dv_orig(idir,:)
      call ztbmv(uplo, trans, diag, nsize, ksize, a, lda_val, x, incx_val)
      x_forward = x
      
      ! Backward perturbation: f(x - h * direction)
      a = a_orig - cmplx(h, 0.0) * a_dv_orig(idir,:,:)
      x = x_orig - cmplx(h, 0.0) * x_dv_orig(idir,:)
      call ztbmv(uplo, trans, diag, nsize, ksize, a, lda_val, x, incx_val)
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
        error_bound = 1.0e-5 + 1.0e-5 * abs_reference
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
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in vector derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Vector derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_ztbmv_vector_forward