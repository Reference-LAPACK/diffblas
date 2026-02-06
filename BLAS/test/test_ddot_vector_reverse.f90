! Test program for DDOT vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirsmax=4

program test_ddot_vector_reverse
  implicit none
  include 'DIFFSIZES.inc'

  real(8), external :: ddot
  external :: ddot_bv

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, k  ! Loop counters
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  integer :: nsize
  real(8), dimension(4) :: dx
  integer :: incx_val
  real(8), dimension(4) :: dy
  integer :: incy_val

  ! Adjoint variables (reverse vector mode)
  ! In reverse mode: output adjoints are INPUT (cotangents/seeds)
  !                  input adjoints are OUTPUT (computed gradients)
  real(8), dimension(nbdirsmax,4) :: dxb
  real(8), dimension(nbdirsmax,4) :: dyb
  real(8), dimension(nbdirsmax) :: ddotb

  ! Storage for original cotangents (for INOUT parameters in VJP verification)
  real(8), dimension(nbdirsmax) :: ddotb_orig

  ! Storage for original values (for VJP verification)
  real(8), dimension(4) :: dx_orig
  real(8), dimension(4) :: dy_orig

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
  call random_number(dx)
  dx = dx * 2.0 - 1.0
  incx_val = 1
  call random_number(dy)
  dy = dy * 2.0 - 1.0
  incy_val = 1

  ! Store original primal values
  dx_orig = dx
  dy_orig = dy

  ! Initialize output adjoints (cotangents) with random values for each direction
  ! These are the 'seeds' for reverse mode
  ! Initialize function result adjoint (output cotangent)
  do k = 1, nbdirsmax
    call random_number(ddotb(k))
    ddotb(k) = ddotb(k) * 2.0 - 1.0
  end do

  ! Initialize input adjoints to zero (they will be computed)
  ! Note: Inout parameters are skipped - they already have output adjoints initialized
  dxb = 0.0
  dyb = 0.0

  ! Save original cotangent seeds for OUTPUT/INOUT parameters (before function call)
  ddotb_orig = ddotb

  ! Set ISIZE globals required by differentiated routine (dimension 2 of arrays).
  ! Differentiated code checks they are set via check_ISIZE*_initialized.
  call set_ISIZE1OFDx(max_size)
  call set_ISIZE1OFDy(max_size)

  ! Call reverse vector mode differentiated function
  call ddot_bv(nsize, dx, dxb, incx_val, dy, dyb, incy_val, ddotb, nbdirsmax)

  ! Reset ISIZE globals to uninitialized (-1) for completeness
  call set_ISIZE1OFDx(-1)
  call set_ISIZE1OFDy(-1)

  ! VJP Verification using finite differences
  call check_vjp_numerically()

  write(*,*) ''
  write(*,*) 'Test completed successfully'

contains

  subroutine check_vjp_numerically()
    implicit none
    
    ! Direction vectors for VJP testing
    real(8), dimension(4) :: dx_dir
    real(8), dimension(4) :: dy_dir
    real(8) :: ddot_plus, ddot_minus
    
    max_error = 0.0d0
    has_large_errors = .false.
    
    write(*,*) 'Function calls completed successfully'
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Test each differentiation direction separately
    do k = 1, nbdirsmax
      
      ! Initialize random direction vectors for all inputs
      call random_number(dx_dir)
      dx_dir = dx_dir * 2.0 - 1.0
      call random_number(dy_dir)
      dy_dir = dy_dir * 2.0 - 1.0
      
      ! Forward perturbation: f(x + h*dir)
      dx = dx_orig + h * dx_dir
      dy = dy_orig + h * dy_dir
      ddot_plus = ddot(nsize, dx, incx_val, dy, incy_val)
      
      ! Backward perturbation: f(x - h*dir)
      dx = dx_orig - h * dx_dir
      dy = dy_orig - h * dy_dir
      ddot_minus = ddot(nsize, dx, incx_val, dy, incy_val)
      
      ! Compute central differences and VJP verification
      ! VJP check: direction^T @ adjoint should equal finite difference
      
      ! Compute finite difference VJP (central difference)
      ! For functions: vjp_fd = adjoint * central_diff
      vjp_fd = ddotb(k) * (ddot_plus - ddot_minus) / (2.0d0 * h)
      
      ! Right side: direction^T @ computed_adjoint (with sorted summation)
      ! For INOUT parameters: use cb directly (it contains the computed input adjoint after reverse pass)
      ! For pure inputs: use adjoint directly
      vjp_ad = 0.0d0
      ! Compute and sort products for dy
      n_products = n
      do i = 1, n
        temp_products(i) = dy_dir(i) * dyb(k,i)
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_ad = vjp_ad + temp_products(i)
      end do
      ! Compute and sort products for dx
      n_products = n
      do i = 1, n
        temp_products(i) = dx_dir(i) * dxb(k,i)
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

end program test_ddot_vector_reverse