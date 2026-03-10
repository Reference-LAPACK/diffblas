! Test program for SSWAP vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=4

program test_sswap_vector_reverse
  implicit none
  integer, parameter :: nbdirs = 4

  external :: sswap
  external :: sswap_bv

  ! Test parameters
  integer :: n  ! Current size (set in loop)
  integer, parameter :: max_size = 100  ! Maximum array dimension (multi-size: 1,4,40,100)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, k  ! Loop counters
  integer :: test_sizes(1), itest
  logical :: passed, all_passed
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  integer :: nsize
  real(4), dimension(max_size) :: sx
  integer :: incx_val
  real(4), dimension(max_size) :: sy
  integer :: incy_val

  ! Adjoint variables (reverse vector mode)
  ! In reverse mode: output adjoints are INPUT (cotangents/seeds)
  !                  input adjoints are OUTPUT (computed gradients)
  real(4), dimension(nbdirs,max_size) :: sxb
  real(4), dimension(nbdirs,max_size) :: syb

  ! Storage for original cotangents (for INOUT parameters in VJP verification)
  real(4), dimension(nbdirs,max_size) :: sxb_orig
  real(4), dimension(nbdirs,max_size) :: syb_orig

  ! Storage for original values (for VJP verification)
  real(4), dimension(max_size) :: sx_orig
  real(4), dimension(max_size) :: sy_orig

  ! Variables for VJP verification via finite differences
  real(4), parameter :: h = 1.0e-3
  real(4) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
  logical :: has_large_errors
  real(4), dimension(max_size*max_size) :: temp_products  ! For sorted summation
  integer :: n_products

  ! Initialize random seed for reproducibility
  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing SSWAP (Vector Reverse, multi-size: n = 4)'
  all_passed = .true.
  do itest = 1, 1
    n = test_sizes(itest)
    write(*,*) 'Testing SSWAP (Vector Reverse, n =', n, ')'

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

  ! Initialize output adjoints (cotangents) with random values for each direction
  ! These are the 'seeds' for reverse mode
  do k = 1, nbdirs
    call random_number(sxb(k,:))
    sxb(k,:) = sxb(k,:) * 2.0 - 1.0
  end do
  do k = 1, nbdirs
    call random_number(syb(k,:))
    syb(k,:) = syb(k,:) * 2.0 - 1.0
  end do

  ! Initialize input adjoints to zero (they will be computed)
  ! Note: Inout parameters are skipped - they already have output adjoints initialized

  ! Save original cotangent seeds for OUTPUT/INOUT parameters (before function call)
  sxb_orig = sxb
  syb_orig = syb

  ! Call reverse vector mode differentiated function
  call sswap_bv(nsize, sx, sxb, incx_val, sy, syb, incy_val, nbdirs)

  ! VJP Verification using finite differences
  call check_vjp_numerically(passed)
  all_passed = all_passed .and. passed
  end do
  if (all_passed) then
    write(*,*) 'PASS: Vector reverse mode - all sizes completed successfully'
  else
    write(*,*) 'FAIL: Vector reverse mode - one or more sizes had derivative errors'
  end if

contains

  subroutine check_vjp_numerically(passed)
    implicit none
    logical, intent(out) :: passed
    
    ! Direction vectors for VJP testing
    real(4), dimension(max_size) :: sx_dir
    real(4), dimension(max_size) :: sy_dir
    real(4), dimension(max_size) :: sx_plus, sx_minus, sx_central_diff
    real(4), dimension(max_size) :: sy_plus, sy_minus, sy_central_diff
    
    max_error = 0.0d0
    has_large_errors = .false.
    
    write(*,*) 'Function calls completed successfully'
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Test each differentiation direction separately
    do k = 1, nbdirs
      
      ! Initialize random direction vectors for all inputs
      call random_number(sx_dir)
      sx_dir = sx_dir * 2.0 - 1.0
      call random_number(sy_dir)
      sy_dir = sy_dir * 2.0 - 1.0
      
      ! Forward perturbation: f(x + h*dir)
      sx = sx_orig + h * sx_dir
      sy = sy_orig + h * sy_dir
      call sswap(nsize, sx, incx_val, sy, incy_val)
      sx_plus = sx
      sy_plus = sy
      
      ! Backward perturbation: f(x - h*dir)
      sx = sx_orig - h * sx_dir
      sy = sy_orig - h * sy_dir
      call sswap(nsize, sx, incx_val, sy, incy_val)
      sx_minus = sx
      sy_minus = sy
      
      ! Compute central differences and VJP verification
      ! VJP check: direction^T @ adjoint should equal finite difference
      
      ! Compute central differences: (f(x+h*dir) - f(x-h*dir)) / (2h)
      sx_central_diff = (sx_plus - sx_minus) / (2.0 * h)
      sy_central_diff = (sy_plus - sy_minus) / (2.0 * h)
      
      ! VJP verification:
      ! cotangent^T @ central_diff should equal direction^T @ computed_adjoint
      ! Left side: cotangent^T @ Jacobian @ direction (via finite differences, with sorted summation)
      vjp_fd = 0.0
      ! Compute and sort products for sx (FD)
      n_products = n
      do i = 1, n
        temp_products(i) = sxb_orig(k,i) * sx_central_diff(i)
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_fd = vjp_fd + temp_products(i)
      end do
      ! Compute and sort products for sy (FD)
      n_products = n
      do i = 1, n
        temp_products(i) = syb_orig(k,i) * sy_central_diff(i)
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
        temp_products(i) = sx_dir(i) * sxb(k,i)
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_ad = vjp_ad + temp_products(i)
      end do
      ! Compute and sort products for sy
      n_products = n
      do i = 1, n
        temp_products(i) = sy_dir(i) * syb(k,i)
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
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
    passed = .not. has_large_errors
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

end program test_sswap_vector_reverse