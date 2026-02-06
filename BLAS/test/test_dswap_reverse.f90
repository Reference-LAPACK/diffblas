! Test program for DSWAP reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Verification uses VJP methodology with finite differences

program test_dswap_reverse
  implicit none

  external :: dswap
  external :: dswap_b

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension (rows/cols of matrices)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions

  integer :: nsize
  real(8), dimension(max_size) :: dx
  integer :: incx_val
  real(8), dimension(max_size) :: dy
  integer :: incy_val

  ! Adjoint variables (reverse mode)
  ! In reverse mode: output adjoints are INPUT (cotangents/seeds)
  !                  input adjoints are OUTPUT (computed gradients)
  real(8), dimension(max_size) :: dxb
  real(8), dimension(max_size) :: dyb

  ! Storage for original values (for VJP verification)
  real(8), dimension(max_size) :: dx_orig
  real(8), dimension(max_size) :: dy_orig

  ! Variables for VJP verification via finite differences
  real(8), dimension(max_size) :: dy_plus, dy_minus
  real(8), dimension(max_size) :: dx_plus, dx_minus

  ! Saved cotangents (output adjoints) for VJP verification
  real(8), dimension(max_size) :: dyb_orig
  real(8), dimension(max_size) :: dxb_orig
  real(8), parameter :: h = 1.0e-7
  real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
  logical :: has_large_errors
  integer :: i, j
  real(8), dimension(max_size*max_size) :: temp_products  ! For sorted summation
  integer :: n_products

  ! Initialize random seed for reproducibility
  integer :: seed_array(33)
  seed_array = 42
  call random_seed(put=seed_array)

  ! Initialize primal values
  nsize = n
  call random_number(dx)
  dx = dx * 2.0d0 - 1.0d0
  incx_val = 1
  call random_number(dy)
  dy = dy * 2.0d0 - 1.0d0
  incy_val = 1

  ! Store original primal values
  dx_orig = dx
  dy_orig = dy

  write(*,*) 'Testing DSWAP'

  ! Initialize output adjoints (cotangents) with random values
  ! These are the 'seeds' for reverse mode
  call random_number(dyb)
  dyb = dyb * 2.0d0 - 1.0d0
  call random_number(dxb)
  dxb = dxb * 2.0d0 - 1.0d0

  ! Save output adjoints (cotangents) for VJP verification
  ! Note: output adjoints may be modified by reverse mode function
  dyb_orig = dyb
  dxb_orig = dxb

  ! Initialize input adjoints to zero (they will be computed)

  ! Call reverse mode differentiated function
  call dswap_b(nsize, dx, dxb, incx_val, dy, dyb, incy_val)

  ! VJP Verification using finite differences
  ! For reverse mode, we verify: cotangent^T @ J @ direction = direction^T @ adjoint
  ! Equivalently: cotangent^T @ (f(x+h*dir) - f(x-h*dir))/(2h) should equal dir^T @ computed_adjoint
  call check_vjp_numerically()

  write(*,*) ''
  write(*,*) 'Test completed successfully'

contains

  subroutine check_vjp_numerically()
    implicit none
    
    ! Direction vectors for VJP testing (like tangents in forward mode)
    real(8), dimension(max_size) :: dx_dir
    real(8), dimension(max_size) :: dy_dir
    
    real(8), dimension(max_size) :: dy_central_diff
    real(8), dimension(max_size) :: dx_central_diff
    
    max_error = 0.0d0
    has_large_errors = .false.
    
    write(*,*) 'Function calls completed successfully'
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Initialize random direction vectors for all inputs
    call random_number(dx_dir)
    dx_dir = dx_dir * 2.0d0 - 1.0d0
    call random_number(dy_dir)
    dy_dir = dy_dir * 2.0d0 - 1.0d0
    
    ! Forward perturbation: f(x + h*dir)
    dx = dx_orig + h * dx_dir
    dy = dy_orig + h * dy_dir
    call dswap(nsize, dx, incx_val, dy, incy_val)
    dy_plus = dy
    dx_plus = dx
    
    ! Backward perturbation: f(x - h*dir)
    dx = dx_orig - h * dx_dir
    dy = dy_orig - h * dy_dir
    call dswap(nsize, dx, incx_val, dy, incy_val)
    dy_minus = dy
    dx_minus = dx
    
    ! Compute central differences: (f(x+h*dir) - f(x-h*dir)) / (2h)
    dy_central_diff = (dy_plus - dy_minus) / (2.0d0 * h)
    dx_central_diff = (dx_plus - dx_minus) / (2.0d0 * h)
    
    ! VJP verification:
    ! cotangent^T @ central_diff should equal direction^T @ computed_adjoint
    ! Left side: cotangent^T @ Jacobian @ direction (via finite differences, with sorted summation)
    vjp_fd = 0.0d0
    ! Compute and sort products for dy (FD)
    n_products = n
    do i = 1, n
      temp_products(i) = dyb_orig(i) * dy_central_diff(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do
    ! Compute and sort products for dx (FD)
    n_products = n
    do i = 1, n
      temp_products(i) = dxb_orig(i) * dx_central_diff(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do
    
    ! Right side: direction^T @ computed_adjoint (with sorted summation)
    ! For INOUT parameters: use cb directly (it contains the computed input adjoint after reverse pass)
    ! For pure inputs: use adjoint directly
    vjp_ad = 0.0d0
    ! Compute and sort products for dx
    n_products = n
    do i = 1, n
      temp_products(i) = dx_dir(i) * dxb(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    ! Compute and sort products for dy
    n_products = n
    do i = 1, n
      temp_products(i) = dy_dir(i) * dyb(i)
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
    
    
    if (abs_reference > 1.0e-10) then
      relative_error = abs_error / abs_reference
    else
      relative_error = abs_error
    end if
    max_error = relative_error
    
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

end program test_dswap_reverse