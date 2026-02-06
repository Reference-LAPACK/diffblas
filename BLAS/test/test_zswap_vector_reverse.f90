! Test program for ZSWAP vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirsmax=4

program test_zswap_vector_reverse
  implicit none
  include 'DIFFSIZES.inc'

  external :: zswap
  external :: zswap_bv

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, k  ! Loop counters
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  integer :: nsize
  complex(8), dimension(max_size) :: zx
  integer :: incx_val
  complex(8), dimension(max_size) :: zy
  integer :: incy_val

  ! Adjoint variables (reverse vector mode)
  ! In reverse mode: output adjoints are INPUT (cotangents/seeds)
  !                  input adjoints are OUTPUT (computed gradients)
  complex(8), dimension(nbdirsmax,max_size) :: zxb
  complex(8), dimension(nbdirsmax,max_size) :: zyb

  ! Storage for original cotangents (for INOUT parameters in VJP verification)
  complex(8), dimension(nbdirsmax,max_size) :: zxb_orig
  complex(8), dimension(nbdirsmax,max_size) :: zyb_orig

  ! Storage for original values (for VJP verification)
  complex(8), dimension(max_size) :: zx_orig
  complex(8), dimension(max_size) :: zy_orig

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
  nsize = n
  do i = 1, n
    call random_number(temp_real)
    call random_number(temp_imag)
    zx(i) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
  end do
  incx_val = 1
  do i = 1, n
    call random_number(temp_real)
    call random_number(temp_imag)
    zy(i) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
  end do
  incy_val = 1

  ! Store original primal values
  zx_orig = zx
  zy_orig = zy

  ! Initialize output adjoints (cotangents) with random values for each direction
  ! These are the 'seeds' for reverse mode
  do k = 1, nbdirsmax
    do i = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      zxb(k,i) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
    end do
  end do
  do k = 1, nbdirsmax
    do i = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      zyb(k,i) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
    end do
  end do

  ! Initialize input adjoints to zero (they will be computed)
  ! Note: Inout parameters are skipped - they already have output adjoints initialized

  ! Save original cotangent seeds for OUTPUT/INOUT parameters (before function call)
  zxb_orig = zxb
  zyb_orig = zyb

  ! Call reverse vector mode differentiated function
  call zswap_bv(nsize, zx, zxb, incx_val, zy, zyb, incy_val, nbdirsmax)

  ! VJP Verification using finite differences
  call check_vjp_numerically()

  write(*,*) ''
  write(*,*) 'Test completed successfully'

contains

  subroutine check_vjp_numerically()
    implicit none
    
    ! Direction vectors for VJP testing
    complex(8), dimension(max_size) :: zx_dir
    complex(8), dimension(max_size) :: zy_dir
    complex(8), dimension(max_size) :: zx_plus, zx_minus, zx_central_diff
    complex(8), dimension(max_size) :: zy_plus, zy_minus, zy_central_diff
    
    max_error = 0.0d0
    has_large_errors = .false.
    
    write(*,*) 'Function calls completed successfully'
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Test each differentiation direction separately
    do k = 1, nbdirsmax
      
      ! Initialize random direction vectors for all inputs
      do i = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        zx_dir(i) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
      end do
      do i = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        zy_dir(i) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
      end do
      
      ! Forward perturbation: f(x + h*dir)
      zx = zx_orig + cmplx(h, 0.0) * zx_dir
      zy = zy_orig + cmplx(h, 0.0) * zy_dir
      call zswap(nsize, zx, incx_val, zy, incy_val)
      zx_plus = zx
      zy_plus = zy
      
      ! Backward perturbation: f(x - h*dir)
      zx = zx_orig - cmplx(h, 0.0) * zx_dir
      zy = zy_orig - cmplx(h, 0.0) * zy_dir
      call zswap(nsize, zx, incx_val, zy, incy_val)
      zx_minus = zx
      zy_minus = zy
      
      ! Compute central differences and VJP verification
      ! VJP check: direction^T @ adjoint should equal finite difference
      
      ! Compute central differences: (f(x+h*dir) - f(x-h*dir)) / (2h)
      zx_central_diff = (zx_plus - zx_minus) / (2.0d0 * h)
      zy_central_diff = (zy_plus - zy_minus) / (2.0d0 * h)
      
      ! VJP verification:
      ! cotangent^T @ central_diff should equal direction^T @ computed_adjoint
      ! Left side: cotangent^T @ Jacobian @ direction (via finite differences, with sorted summation)
      vjp_fd = 0.0d0
      ! Compute and sort products for zx (FD)
      n_products = n
      do i = 1, n
        temp_products(i) = real(conjg(zxb_orig(k,i)) * zx_central_diff(i))
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_fd = vjp_fd + temp_products(i)
      end do
      ! Compute and sort products for zy (FD)
      n_products = n
      do i = 1, n
        temp_products(i) = real(conjg(zyb_orig(k,i)) * zy_central_diff(i))
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_fd = vjp_fd + temp_products(i)
      end do
      
      ! Right side: direction^T @ computed_adjoint (with sorted summation)
      ! For INOUT parameters: use cb directly (it contains the computed input adjoint after reverse pass)
      ! For pure inputs: use adjoint directly
      vjp_ad = 0.0d0
      ! Compute and sort products for zx
      n_products = n
      do i = 1, n
        temp_products(i) = real(conjg(zx_dir(i)) * zxb(k,i))
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_ad = vjp_ad + temp_products(i)
      end do
      ! Compute and sort products for zy
      n_products = n
      do i = 1, n
        temp_products(i) = real(conjg(zy_dir(i)) * zyb(k,i))
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

end program test_zswap_vector_reverse