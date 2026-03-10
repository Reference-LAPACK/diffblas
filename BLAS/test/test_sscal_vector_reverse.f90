! Test program for SSCAL vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=4

program test_sscal_vector_reverse
  implicit none
  integer, parameter :: nbdirs = 4

  external :: sscal
  external :: sscal_bv

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
  real(4) :: sa
  real(4), dimension(max_size) :: sx
  integer :: incx_val

  ! Adjoint variables (reverse vector mode)
  ! In reverse mode: output adjoints are INPUT (cotangents/seeds)
  !                  input adjoints are OUTPUT (computed gradients)
  real(4), dimension(nbdirs) :: sab
  real(4), dimension(nbdirs,max_size) :: sxb

  ! Storage for original cotangents (for INOUT parameters in VJP verification)
  real(4), dimension(nbdirs,max_size) :: sxb_orig

  ! Storage for original values (for VJP verification)
  real(4) :: sa_orig
  real(4), dimension(max_size) :: sx_orig

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
  write(*,*) 'Testing SSCAL (Vector Reverse, multi-size: n = 4)'
  all_passed = .true.
  do itest = 1, 1
    n = test_sizes(itest)
    write(*,*) 'Testing SSCAL (Vector Reverse, n =', n, ')'

  ! Initialize primal values
  nsize = n
  call random_number(sa)
  sa = sa * 2.0 - 1.0
  call random_number(sx)
  sx = sx * 2.0 - 1.0
  incx_val = 1

  ! Store original primal values
  sa_orig = sa
  sx_orig = sx

  ! Initialize output adjoints (cotangents) with random values for each direction
  ! These are the 'seeds' for reverse mode
  do k = 1, nbdirs
    call random_number(sxb(k,:))
    sxb(k,:) = sxb(k,:) * 2.0 - 1.0
  end do

  ! Initialize input adjoints to zero (they will be computed)
  ! Note: Inout parameters are skipped - they already have output adjoints initialized
  sab = 0.0

  ! Save original cotangent seeds for OUTPUT/INOUT parameters (before function call)
  sxb_orig = sxb

  ! Call reverse vector mode differentiated function
  call sscal_bv(nsize, sa, sab, sx, sxb, incx_val, nbdirs)

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
    real(4) :: sa_dir
    real(4), dimension(max_size) :: sx_dir
    real(4), dimension(max_size) :: sx_plus, sx_minus, sx_central_diff
    
    max_error = 0.0d0
    has_large_errors = .false.
    
    write(*,*) 'Function calls completed successfully'
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Test each differentiation direction separately
    do k = 1, nbdirs
      
      ! Initialize random direction vectors for all inputs
      call random_number(sa_dir)
      sa_dir = sa_dir * 2.0 - 1.0
      call random_number(sx_dir)
      sx_dir = sx_dir * 2.0 - 1.0
      
      ! Forward perturbation: f(x + h*dir)
      sa = sa_orig + h * sa_dir
      sx = sx_orig + h * sx_dir
      call sscal(nsize, sa, sx, incx_val)
      sx_plus = sx
      
      ! Backward perturbation: f(x - h*dir)
      sa = sa_orig - h * sa_dir
      sx = sx_orig - h * sx_dir
      call sscal(nsize, sa, sx, incx_val)
      sx_minus = sx
      
      ! Compute central differences and VJP verification
      ! VJP check: direction^T @ adjoint should equal finite difference
      
      ! Compute central differences: (f(x+h*dir) - f(x-h*dir)) / (2h)
      sx_central_diff = (sx_plus - sx_minus) / (2.0 * h)
      
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
      vjp_ad = vjp_ad + sa_dir * sab(k)
      
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

end program test_sscal_vector_reverse