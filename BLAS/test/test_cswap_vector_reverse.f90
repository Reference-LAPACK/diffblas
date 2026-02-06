! Test program for CSWAP vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirsmax=4

program test_cswap_vector_reverse
  implicit none
  include 'DIFFSIZES.inc'

  external :: cswap
  external :: cswap_bv

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, k  ! Loop counters
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  integer :: nsize
  complex(4), dimension(max_size) :: cx
  integer :: incx_val
  complex(4), dimension(max_size) :: cy
  integer :: incy_val

  ! Adjoint variables (reverse vector mode)
  ! In reverse mode: output adjoints are INPUT (cotangents/seeds)
  !                  input adjoints are OUTPUT (computed gradients)
  complex(4), dimension(nbdirsmax,max_size) :: cxb
  complex(4), dimension(nbdirsmax,max_size) :: cyb

  ! Storage for original cotangents (for INOUT parameters in VJP verification)
  complex(4), dimension(nbdirsmax,max_size) :: cxb_orig
  complex(4), dimension(nbdirsmax,max_size) :: cyb_orig

  ! Storage for original values (for VJP verification)
  complex(4), dimension(max_size) :: cx_orig
  complex(4), dimension(max_size) :: cy_orig

  ! Variables for VJP verification via finite differences
  real(4), parameter :: h = 1.0e-3
  real(4) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
  logical :: has_large_errors
  real(4), dimension(max_size*max_size) :: temp_products  ! For sorted summation
  integer :: n_products

  ! Initialize random seed for reproducibility
  seed_array = 42
  call random_seed(put=seed_array)

  ! Initialize primal values
  nsize = n
  do i = 1, n
    call random_number(temp_real)
    call random_number(temp_imag)
    cx(i) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
  end do
  incx_val = 1
  do i = 1, n
    call random_number(temp_real)
    call random_number(temp_imag)
    cy(i) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
  end do
  incy_val = 1

  ! Store original primal values
  cx_orig = cx
  cy_orig = cy

  ! Initialize output adjoints (cotangents) with random values for each direction
  ! These are the 'seeds' for reverse mode
  do k = 1, nbdirsmax
    do i = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      cxb(k,i) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
    end do
  end do
  do k = 1, nbdirsmax
    do i = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      cyb(k,i) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
    end do
  end do

  ! Initialize input adjoints to zero (they will be computed)
  ! Note: Inout parameters are skipped - they already have output adjoints initialized

  ! Save original cotangent seeds for OUTPUT/INOUT parameters (before function call)
  cxb_orig = cxb
  cyb_orig = cyb

  ! Call reverse vector mode differentiated function
  call cswap_bv(nsize, cx, cxb, incx_val, cy, cyb, incy_val, nbdirsmax)

  ! VJP Verification using finite differences
  call check_vjp_numerically()

  write(*,*) ''
  write(*,*) 'Test completed successfully'

contains

  subroutine check_vjp_numerically()
    implicit none
    
    ! Direction vectors for VJP testing
    complex(4), dimension(max_size) :: cx_dir
    complex(4), dimension(max_size) :: cy_dir
    complex(4), dimension(max_size) :: cx_plus, cx_minus, cx_central_diff
    complex(4), dimension(max_size) :: cy_plus, cy_minus, cy_central_diff
    
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
        cx_dir(i) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
      end do
      do i = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        cy_dir(i) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
      end do
      
      ! Forward perturbation: f(x + h*dir)
      cx = cx_orig + cmplx(h, 0.0) * cx_dir
      cy = cy_orig + cmplx(h, 0.0) * cy_dir
      call cswap(nsize, cx, incx_val, cy, incy_val)
      cx_plus = cx
      cy_plus = cy
      
      ! Backward perturbation: f(x - h*dir)
      cx = cx_orig - cmplx(h, 0.0) * cx_dir
      cy = cy_orig - cmplx(h, 0.0) * cy_dir
      call cswap(nsize, cx, incx_val, cy, incy_val)
      cx_minus = cx
      cy_minus = cy
      
      ! Compute central differences and VJP verification
      ! VJP check: direction^T @ adjoint should equal finite difference
      
      ! Compute central differences: (f(x+h*dir) - f(x-h*dir)) / (2h)
      cx_central_diff = (cx_plus - cx_minus) / (2.0 * h)
      cy_central_diff = (cy_plus - cy_minus) / (2.0 * h)
      
      ! VJP verification:
      ! cotangent^T @ central_diff should equal direction^T @ computed_adjoint
      ! Left side: cotangent^T @ Jacobian @ direction (via finite differences, with sorted summation)
      vjp_fd = 0.0
      ! Compute and sort products for cx (FD)
      n_products = n
      do i = 1, n
        temp_products(i) = real(conjg(cxb_orig(k,i)) * cx_central_diff(i))
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_fd = vjp_fd + temp_products(i)
      end do
      ! Compute and sort products for cy (FD)
      n_products = n
      do i = 1, n
        temp_products(i) = real(conjg(cyb_orig(k,i)) * cy_central_diff(i))
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_fd = vjp_fd + temp_products(i)
      end do
      
      ! Right side: direction^T @ computed_adjoint (with sorted summation)
      ! For INOUT parameters: use cb directly (it contains the computed input adjoint after reverse pass)
      ! For pure inputs: use adjoint directly
      vjp_ad = 0.0
      ! Compute and sort products for cx
      n_products = n
      do i = 1, n
        temp_products(i) = real(conjg(cx_dir(i)) * cxb(k,i))
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_ad = vjp_ad + temp_products(i)
      end do
      ! Compute and sort products for cy
      n_products = n
      do i = 1, n
        temp_products(i) = real(conjg(cy_dir(i)) * cyb(k,i))
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_ad = vjp_ad + temp_products(i)
      end do
      
      ! Error check: |vjp_fd - vjp_ad| > atol + rtol * |vjp_ad|
      abs_error = abs(vjp_fd - vjp_ad)
      abs_reference = abs(vjp_ad)
      error_bound = 1.0e-3 + 1.0e-3 * abs_reference
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
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
    
  end subroutine check_vjp_numerically

  subroutine sort_array(arr, n)
    implicit none
    integer, intent(in) :: n
    real(4), dimension(n), intent(inout) :: arr
    integer :: i, j, min_idx
    real(4) :: temp
    
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

end program test_cswap_vector_reverse