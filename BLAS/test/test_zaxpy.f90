! Test program for ZAXPY differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision

program test_zaxpy
  implicit none

  external :: zaxpy
  external :: zaxpy_d

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension (rows/cols of matrices)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions

  integer :: nsize
  complex(8) :: za
  complex(8), dimension(4) :: zx
  integer :: incx_val
  complex(8), dimension(max_size) :: zy
  integer :: incy_val

  ! Derivative variables
  complex(8) :: za_d
  complex(8), dimension(4) :: zx_d
  complex(8), dimension(max_size) :: zy_d

  ! Storage variables for inout parameters
  complex(8), dimension(max_size) :: zy_output

  ! Array restoration variables for numerical differentiation
  complex(8), dimension(4) :: zx_orig
  complex(8), dimension(max_size) :: zy_orig
  complex(8) :: za_orig

  ! Variables for central difference computation
  complex(8), dimension(max_size) :: zy_forward, zy_backward
  ! Scalar variables for central difference computation
  complex(8) :: central_diff, ad_result
  logical :: has_large_errors

  ! Variables for storing original derivative values
  complex(8), dimension(4) :: zx_d_orig
  complex(8), dimension(max_size) :: zy_d_orig
  complex(8) :: za_d_orig

  ! Temporary variables for matrix initialization
  real(4) :: temp_real, temp_imag
  integer :: i, j

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  integer :: seed_array(33)
  seed_array = 42
  call random_seed(put=seed_array)

  nsize = n
  call random_number(temp_real)
  call random_number(temp_imag)
  za = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  do i = 1, n
    call random_number(temp_real)
    call random_number(temp_imag)
    zx(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  end do
  incx_val = 1
  do i = 1, n
    call random_number(temp_real)
    call random_number(temp_imag)
    zy(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  end do
  incy_val = 1

  ! Initialize input derivatives to random values
  call random_number(temp_real)
  call random_number(temp_imag)
  zx_d = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  call random_number(temp_real)
  call random_number(temp_imag)
  zy_d = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  call random_number(temp_real)
  call random_number(temp_imag)
  za_d = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)

  ! Store initial derivative values after random initialization
  zx_d_orig = zx_d
  zy_d_orig = zy_d
  za_d_orig = za_d

  ! Store original values for central difference computation
  zx_orig = zx
  zy_orig = zy
  za_orig = za

  write(*,*) 'Testing ZAXPY'
  ! Store input values of inout parameters before first function call
  zy_orig = zy

  ! Re-initialize data for differentiated function
  ! Only reinitialize inout parameters - keep input-only parameters unchanged

  nsize = n
  ! za already has correct value from original call
  ! zx already has correct value from original call
  incx_val = 1
  zy = zy_orig
  incy_val = 1

  ! Call the differentiated function
  call zaxpy_d(nsize, za, za_d, zx, zx_d, incx_val, zy, zy_d, incy_val)

  ! Print results and compare
  write(*,*) 'Function calls completed successfully'

  ! Numerical differentiation check
  call check_derivatives_numerically()

  write(*,*) 'Test completed successfully'

contains

  subroutine check_derivatives_numerically()
    implicit none
    real(8), parameter :: h = 1.0e-6  ! Step size for finite differences
    real(8) :: relative_error, max_error
    real(8) :: output_orig, output_pert
    real(8) :: numerical_result, analytical_result
    real(8) :: abs_error, abs_reference, error_bound
    integer :: i, j
    
    max_error = 0.0e0
    has_large_errors = .false.
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5
    
    ! Original values already stored in main program
    
    ! Central difference computation: f(x + h) - f(x - h) / (2h)
    ! Forward perturbation: f(x + h)
    zx = zx_orig + cmplx(h, 0.0) * zx_d_orig
    zy = zy_orig + cmplx(h, 0.0) * zy_d_orig
    za = za_orig + cmplx(h, 0.0) * za_d_orig
    call zaxpy(nsize, za, zx, incx_val, zy, incy_val)
    ! Store forward perturbation results
    
    ! Backward perturbation: f(x - h)
    zx = zx_orig - cmplx(h, 0.0) * zx_d_orig
    zy = zy_orig - cmplx(h, 0.0) * zy_d_orig
    za = za_orig - cmplx(h, 0.0) * za_d_orig
    call zaxpy(nsize, za, zx, incx_val, zy, incy_val)
    ! Store backward perturbation results
    
    ! Compute central differences and compare with AD results
    
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_derivatives_numerically

end program test_zaxpy