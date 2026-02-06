! Test program for ZGERU vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirsmax=4

program test_zgeru_vector_reverse
  implicit none
  include 'DIFFSIZES.inc'

  external :: zgeru
  external :: zgeru_bv

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, k  ! Loop counters
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

  ! Adjoint variables (reverse vector mode)
  ! In reverse mode: output adjoints are INPUT (cotangents/seeds)
  !                  input adjoints are OUTPUT (computed gradients)
  complex(8), dimension(nbdirsmax) :: alphab
  complex(8), dimension(nbdirsmax,max_size) :: xb
  complex(8), dimension(nbdirsmax,max_size) :: yb
  complex(8), dimension(nbdirsmax,max_size,max_size) :: ab

  ! Storage for original cotangents (for INOUT parameters in VJP verification)
  complex(8), dimension(nbdirsmax,max_size,max_size) :: ab_orig

  ! Storage for original values (for VJP verification)
  complex(8) :: alpha_orig
  complex(8), dimension(max_size) :: x_orig
  complex(8), dimension(max_size) :: y_orig
  complex(8), dimension(max_size,max_size) :: a_orig

  ! Variables for VJP verification via finite differences
  real(8), parameter :: h = 1.0e-7
  real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
  logical :: has_large_errors
  real(8), dimension(max_size*max_size) :: temp_products  ! For sorted summation
  integer :: n_products

  ! Initialize random seed for reproducibility
  seed_array = 42
  call random_seed(put=seed_array)

  ! Initialize primal values
  msize = n
  nsize = n
  call random_number(temp_real)
  call random_number(temp_imag)
  alpha = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
  do i = 1, n
    call random_number(temp_real)
    call random_number(temp_imag)
    x(i) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
  end do
  incx_val = 1
  do i = 1, n
    call random_number(temp_real)
    call random_number(temp_imag)
    y(i) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
  end do
  incy_val = 1
  do j = 1, n
    do i = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      a(i,j) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
    end do
  end do
  lda_val = lda

  ! Store original primal values
  alpha_orig = alpha
  x_orig = x
  y_orig = y
  a_orig = a

  ! Initialize output adjoints (cotangents) with random values for each direction
  ! These are the 'seeds' for reverse mode
  do k = 1, nbdirsmax
    do j = 1, n
      do i = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        ab(k,i,j) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
      end do
    end do
  end do

  ! Initialize input adjoints to zero (they will be computed)
  ! Note: Inout parameters are skipped - they already have output adjoints initialized
  alphab = 0.0
  xb = 0.0
  yb = 0.0

  ! Save original cotangent seeds for OUTPUT/INOUT parameters (before function call)
  ab_orig = ab

  ! Set ISIZE globals required by differentiated routine (dimension 2 of arrays).
  ! Differentiated code checks they are set via check_ISIZE*_initialized.
  call set_ISIZE1OFX(max_size)
  call set_ISIZE1OFY(max_size)

  ! Call reverse vector mode differentiated function
  call zgeru_bv(msize, nsize, alpha, alphab, x, xb, incx_val, y, yb, incy_val, a, ab, lda_val, nbdirsmax)

  ! Reset ISIZE globals to uninitialized (-1) for completeness
  call set_ISIZE1OFX(-1)
  call set_ISIZE1OFY(-1)

  ! VJP Verification using finite differences
  call check_vjp_numerically()

  write(*,*) ''
  write(*,*) 'Test completed successfully'

contains

  subroutine check_vjp_numerically()
    implicit none
    
    ! Direction vectors for VJP testing
    complex(8) :: alpha_dir
    complex(8), dimension(max_size) :: x_dir
    complex(8), dimension(max_size) :: y_dir
    complex(8), dimension(max_size,max_size) :: a_dir
    complex(8), dimension(max_size,max_size) :: a_plus, a_minus, a_central_diff
    
    max_error = 0.0d0
    has_large_errors = .false.
    
    write(*,*) 'Function calls completed successfully'
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Test each differentiation direction separately
    do k = 1, nbdirsmax
      
      ! Initialize random direction vectors for all inputs
      call random_number(temp_real)
      call random_number(temp_imag)
      alpha_dir = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
      do i = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        x_dir(i) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
      end do
      do i = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        y_dir(i) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
      end do
      do j = 1, n
        do i = 1, n
          call random_number(temp_real)
          call random_number(temp_imag)
          a_dir(i,j) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
        end do
      end do
      
      ! Forward perturbation: f(x + h*dir)
      alpha = alpha_orig + cmplx(h, 0.0) * alpha_dir
      x = x_orig + cmplx(h, 0.0) * x_dir
      y = y_orig + cmplx(h, 0.0) * y_dir
      a = a_orig + cmplx(h, 0.0) * a_dir
      call zgeru(msize, nsize, alpha, x, incx_val, y, incy_val, a, lda_val)
      a_plus = a
      
      ! Backward perturbation: f(x - h*dir)
      alpha = alpha_orig - cmplx(h, 0.0) * alpha_dir
      x = x_orig - cmplx(h, 0.0) * x_dir
      y = y_orig - cmplx(h, 0.0) * y_dir
      a = a_orig - cmplx(h, 0.0) * a_dir
      call zgeru(msize, nsize, alpha, x, incx_val, y, incy_val, a, lda_val)
      a_minus = a
      
      ! Compute central differences and VJP verification
      ! VJP check: direction^T @ adjoint should equal finite difference
      
      ! Compute central differences: (f(x+h*dir) - f(x-h*dir)) / (2h)
      a_central_diff = (a_plus - a_minus) / (2.0d0 * h)
      
      ! VJP verification:
      ! cotangent^T @ central_diff should equal direction^T @ computed_adjoint
      ! Left side: cotangent^T @ Jacobian @ direction (via finite differences, with sorted summation)
      vjp_fd = 0.0d0
      ! Compute and sort products for a (FD)
      n_products = 0
      do j = 1, n
        do i = 1, n
          n_products = n_products + 1
          temp_products(n_products) = real(conjg(ab_orig(k,i,j)) * a_central_diff(i,j))
        end do
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_fd = vjp_fd + temp_products(i)
      end do
      
      ! Right side: direction^T @ computed_adjoint (with sorted summation)
      ! For INOUT parameters: use cb directly (it contains the computed input adjoint after reverse pass)
      ! For pure inputs: use adjoint directly
      vjp_ad = 0.0d0
      vjp_ad = vjp_ad + real(conjg(alpha_dir) * alphab(k))
      ! Compute and sort products for a
      n_products = 0
      do j = 1, n
        do i = 1, n
          n_products = n_products + 1
          temp_products(n_products) = real(conjg(a_dir(i,j)) * ab(k,i,j))
        end do
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_ad = vjp_ad + temp_products(i)
      end do
      ! Compute and sort products for x
      n_products = n
      do i = 1, n
        temp_products(i) = real(conjg(x_dir(i)) * xb(k,i))
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_ad = vjp_ad + temp_products(i)
      end do
      ! Compute and sort products for y
      n_products = n
      do i = 1, n
        temp_products(i) = real(conjg(y_dir(i)) * yb(k,i))
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_ad = vjp_ad + temp_products(i)
      end do
      
      ! Error check: |vjp_fd - vjp_ad| > atol + rtol * |vjp_ad|
      abs_error = abs(vjp_fd - vjp_ad)
      abs_reference = abs(vjp_ad)
      error_bound = 1.0e-5 + 1.0e-5 * abs_reference
      if (abs_error > error_bound) then
        has_large_errors = .true.
      end if
      
      ! Compute relative error for reporting
      if (abs_reference > 1.0e-10) then
        relative_error = abs_error / abs_reference
      else
        relative_error = abs_error
      end if
      if (relative_error > max_error) max_error = relative_error
    end do
    
    write(*,*) ''
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_vjp_numerically

  subroutine sort_array(arr, n)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(inout) :: arr
    integer :: i, j, min_idx
    real(8) :: temp
    
    ! Simple selection sort
    do i = 1, n-1
      min_idx = i
      do j = i+1, n
        if (abs(arr(j)) < abs(arr(min_idx))) then
          min_idx = j
        end if
      end do
      if (min_idx /= i) then
        temp = arr(i)
        arr(i) = arr(min_idx)
        arr(min_idx) = temp
      end if
    end do
  end subroutine sort_array

end program test_zgeru_vector_reverse