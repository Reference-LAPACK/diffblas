! Test program for ZGEMM vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirsmax=4

program test_zgemm_vector_forward
  implicit none
  include 'DIFFSIZES.inc'

  external :: zgemm
  external :: zgemm_dv

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, idir  ! Loop counters
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  character :: transa
  character :: transb
  integer :: msize
  integer :: nsize
  integer :: ksize
  complex(8) :: alpha
  complex(8), dimension(max_size,max_size) :: a
  integer :: lda_val
  complex(8), dimension(max_size,max_size) :: b
  integer :: ldb_val
  complex(8) :: beta
  complex(8), dimension(max_size,max_size) :: c
  integer :: ldc_val

  ! Vector mode derivative variables (type-promoted)
  ! Scalars become arrays(nbdirsmax), arrays gain extra dimension
  complex(8), dimension(nbdirsmax) :: alpha_dv
  complex(8), dimension(nbdirsmax,max_size,max_size) :: a_dv
  complex(8), dimension(nbdirsmax,max_size,max_size) :: b_dv
  complex(8), dimension(nbdirsmax) :: beta_dv
  complex(8), dimension(nbdirsmax,max_size,max_size) :: c_dv
  ! Declare variables for storing original values
  complex(8) :: alpha_orig
  complex(8), dimension(nbdirsmax) :: alpha_dv_orig
  complex(8), dimension(max_size,max_size) :: a_orig
  complex(8), dimension(nbdirsmax,max_size,max_size) :: a_dv_orig
  complex(8), dimension(max_size,max_size) :: b_orig
  complex(8), dimension(nbdirsmax,max_size,max_size) :: b_dv_orig
  complex(8) :: beta_orig
  complex(8), dimension(nbdirsmax) :: beta_dv_orig
  complex(8), dimension(max_size,max_size) :: c_orig
  complex(8), dimension(nbdirsmax,max_size,max_size) :: c_dv_orig

  ! Initialize test parameters
  msize = n
  nsize = n
  ksize = n
  lda_val = lda
  ldb_val = ldb
  ldc_val = ldc

  ! Initialize test data with random numbers
  ! Initialize random seed for reproducible results
  seed_array = 42
  call random_seed(put=seed_array)

  transa = 'N'
  transb = 'N'
  call random_number(temp_real)
  call random_number(temp_imag)
  alpha = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  do i = 1, max_size
    do j = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      a(i,j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
  end do
  do i = 1, max_size
    do j = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      b(i,j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
  end do
  call random_number(temp_real)
  call random_number(temp_imag)
  beta = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  do i = 1, max_size
    do j = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      c(i,j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
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
      do j = 1, max_size
        call random_number(temp_real)
        call random_number(temp_imag)
        a_dv(idir,i,j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
      end do
    end do
  end do
  do idir = 1, nbdirsmax
    do i = 1, max_size
      do j = 1, max_size
        call random_number(temp_real)
        call random_number(temp_imag)
        b_dv(idir,i,j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
      end do
    end do
  end do
  do idir = 1, nbdirsmax
    call random_number(temp_real)
    call random_number(temp_imag)
    beta_dv(idir) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
  end do
  do idir = 1, nbdirsmax
    do i = 1, max_size
      do j = 1, max_size
        call random_number(temp_real)
        call random_number(temp_imag)
        c_dv(idir,i,j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
      end do
    end do
  end do

  write(*,*) 'Testing ZGEMM (Vector Forward Mode)'
  ! Store original values before any function calls (critical for INOUT parameters)
  alpha_orig = alpha
  alpha_dv_orig = alpha_dv
  a_orig = a
  a_dv_orig = a_dv
  b_orig = b
  b_dv_orig = b_dv
  beta_orig = beta
  beta_dv_orig = beta_dv
  c_orig = c
  c_dv_orig = c_dv

  ! Call the vector mode differentiated function

  call zgemm_dv(transa, transb, msize, nsize, ksize, alpha, alpha_dv, a, a_dv, lda_val, b, b_dv, ldb_val, beta, beta_dv, c, c_dv, ldc_val, nbdirsmax)

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
    complex(8), dimension(max_size,max_size) :: c_forward, c_backward
    
    max_error = 0.0e0
    has_large_errors = .false.
    
    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Number of directions:', nbdirsmax
    
    ! Test each derivative direction separately
    do idir = 1, nbdirsmax
      
      ! Forward perturbation: f(x + h * direction)
      alpha = alpha_orig + cmplx(h, 0.0) * alpha_dv_orig(idir)
      a = a_orig + cmplx(h, 0.0) * a_dv_orig(idir,:,:)
      b = b_orig + cmplx(h, 0.0) * b_dv_orig(idir,:,:)
      beta = beta_orig + cmplx(h, 0.0) * beta_dv_orig(idir)
      c = c_orig + cmplx(h, 0.0) * c_dv_orig(idir,:,:)
      call zgemm(transa, transb, msize, nsize, ksize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
      c_forward = c
      
      ! Backward perturbation: f(x - h * direction)
      alpha = alpha_orig - cmplx(h, 0.0) * alpha_dv_orig(idir)
      a = a_orig - cmplx(h, 0.0) * a_dv_orig(idir,:,:)
      b = b_orig - cmplx(h, 0.0) * b_dv_orig(idir,:,:)
      beta = beta_orig - cmplx(h, 0.0) * beta_dv_orig(idir)
      c = c_orig - cmplx(h, 0.0) * c_dv_orig(idir,:,:)
      call zgemm(transa, transb, msize, nsize, ksize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
      c_backward = c
      
      ! Compute central differences and compare with AD results
      do j = 1, min(2, nsize)  ! Check only first few elements
        do i = 1, min(2, nsize)
          ! Central difference: (f(x+h) - f(x-h)) / (2h)
          central_diff = (c_forward(i,j) - c_backward(i,j)) / (2.0e0 * h)
          ! AD result
          ad_result = c_dv(idir,i,j)
          ! Error check: |a - b| > atol + rtol * |b|
          abs_error = abs(central_diff - ad_result)
          abs_reference = abs(ad_result)
          error_bound = 1.0e-5 + 1.0e-5 * abs_reference
          if (abs_error > error_bound) then
            has_large_errors = .true.
            relative_error = abs_error / max(abs_reference, 1.0e-10)
            write(*,*) '  Large error in direction', idir, ' output C(', i, ',', j, '):'
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

end program test_zgemm_vector_forward