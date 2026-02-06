! Test program for SSWAP reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Verification uses VJP methodology with finite differences

program test_sswap_reverse
  implicit none

  external :: sswap
  external :: sswap_b

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension (rows/cols of matrices)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions

  integer :: nsize
  real(4), dimension(max_size) :: sx
  integer :: incx_val
  real(4), dimension(max_size) :: sy
  integer :: incy_val

  ! Adjoint variables (reverse mode)
  ! In reverse mode: output adjoints are INPUT (cotangents/seeds)
  !                  input adjoints are OUTPUT (computed gradients)
  real(4), dimension(max_size) :: sxb
  real(4), dimension(max_size) :: syb

  ! Storage for original values (for VJP verification)
  real(4), dimension(max_size) :: sx_orig
  real(4), dimension(max_size) :: sy_orig

  ! Variables for VJP verification via finite differences
  real(4), dimension(max_size) :: sy_plus, sy_minus
  real(4), dimension(max_size) :: sx_plus, sx_minus

  ! Saved cotangents (output adjoints) for VJP verification
  real(4), dimension(max_size) :: syb_orig
  real(4), dimension(max_size) :: sxb_orig
  real(4), parameter :: h = 1.0e-3
  real(4) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
  logical :: has_large_errors
  integer :: i, j
  real(4), dimension(max_size*max_size) :: temp_products  ! For sorted summation
  integer :: n_products

  ! Initialize random seed for reproducibility
  integer :: seed_array(33)
  seed_array = 42
  call random_seed(put=seed_array)

  ! Initialize primal values
  nsize = n
  call random_number(sx)
  sx = sx * 2.0 - 1.0
  incx_val = 1
  call random_number(sy)
  sy = sy * 2.0 - 1.0
  incy_val = 1

  ! Store original primal values
  sx_orig = sx
  sy_orig = sy

  write(*,*) 'Testing SSWAP'

  ! Initialize output adjoints (cotangents) with random values
  ! These are the 'seeds' for reverse mode
  call random_number(syb)
  syb = syb * 2.0 - 1.0
  call random_number(sxb)
  sxb = sxb * 2.0 - 1.0

  ! Save output adjoints (cotangents) for VJP verification
  ! Note: output adjoints may be modified by reverse mode function
  syb_orig = syb
  sxb_orig = sxb

  ! Initialize input adjoints to zero (they will be computed)

  ! Call reverse mode differentiated function
  call sswap_b(nsize, sx, sxb, incx_val, sy, syb, incy_val)

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
    real(4), dimension(max_size) :: sx_dir
    real(4), dimension(max_size) :: sy_dir
    
    real(4), dimension(max_size) :: sy_central_diff
    real(4), dimension(max_size) :: sx_central_diff
    
    max_error = 0.0
    has_large_errors = .false.
    
    write(*,*) 'Function calls completed successfully'
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Initialize random direction vectors for all inputs
    call random_number(sx_dir)
    sx_dir = sx_dir * 2.0 - 1.0
    call random_number(sy_dir)
    sy_dir = sy_dir * 2.0 - 1.0
    
    ! Forward perturbation: f(x + h*dir)
    sx = sx_orig + h * sx_dir
    sy = sy_orig + h * sy_dir
    call sswap(nsize, sx, incx_val, sy, incy_val)
    sy_plus = sy
    sx_plus = sx
    
    ! Backward perturbation: f(x - h*dir)
    sx = sx_orig - h * sx_dir
    sy = sy_orig - h * sy_dir
    call sswap(nsize, sx, incx_val, sy, incy_val)
    sy_minus = sy
    sx_minus = sx
    
    ! Compute central differences: (f(x+h*dir) - f(x-h*dir)) / (2h)
    sy_central_diff = (sy_plus - sy_minus) / (2.0d0 * h)
    sx_central_diff = (sx_plus - sx_minus) / (2.0d0 * h)
    
    ! VJP verification:
    ! cotangent^T @ central_diff should equal direction^T @ computed_adjoint
    ! Left side: cotangent^T @ Jacobian @ direction (via finite differences, with sorted summation)
    vjp_fd = 0.0
    ! Compute and sort products for sy (FD)
    n_products = n
    do i = 1, n
      temp_products(i) = syb_orig(i) * sy_central_diff(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do
    ! Compute and sort products for sx (FD)
    n_products = n
    do i = 1, n
      temp_products(i) = sxb_orig(i) * sx_central_diff(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do
    
    ! Right side: direction^T @ computed_adjoint (with sorted summation)
    ! For INOUT parameters: use cb directly (it contains the computed input adjoint after reverse pass)
    ! For pure inputs: use adjoint directly
    vjp_ad = 0.0
    ! Compute and sort products for sx
    n_products = n
    do i = 1, n
      temp_products(i) = sx_dir(i) * sxb(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    ! Compute and sort products for sy
    n_products = n
    do i = 1, n
      temp_products(i) = sy_dir(i) * syb(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    
    ! Error check: |vjp_fd - vjp_ad| > atol + rtol * |vjp_ad|
    abs_error = abs(vjp_fd - vjp_ad)
    abs_reference = abs(vjp_ad)
    error_bound = 2.0e-3 + 2.0e-3 * abs_reference
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
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
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

end program test_sswap_reverse