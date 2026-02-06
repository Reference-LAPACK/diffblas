! Test program for SSPMV vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirsmax=4

program test_sspmv_vector_reverse
  implicit none
  include 'DIFFSIZES.inc'

  external :: sspmv
  external :: sspmv_bv

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, k  ! Loop counters
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  character :: uplo
  integer :: nsize
  real(4) :: alpha
  real(4), dimension((n*(n+1))/2) :: ap
  real(4), dimension(max_size) :: x
  integer :: incx_val
  real(4) :: beta
  real(4), dimension(max_size) :: y
  integer :: incy_val

  ! Adjoint variables (reverse vector mode)
  ! In reverse mode: output adjoints are INPUT (cotangents/seeds)
  !                  input adjoints are OUTPUT (computed gradients)
  real(4), dimension(nbdirsmax) :: alphab
  real(4), dimension(nbdirsmax,(n*(n+1))/2) :: apb
  real(4), dimension(nbdirsmax,max_size) :: xb
  real(4), dimension(nbdirsmax) :: betab
  real(4), dimension(nbdirsmax,max_size) :: yb

  ! Storage for original cotangents (for INOUT parameters in VJP verification)
  real(4), dimension(nbdirsmax,max_size) :: yb_orig

  ! Storage for original values (for VJP verification)
  real(4) :: alpha_orig
  real(4), dimension((n*(n+1))/2) :: ap_orig
  real(4), dimension(max_size) :: x_orig
  real(4) :: beta_orig
  real(4), dimension(max_size) :: y_orig

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
  uplo = 'U'
  nsize = n
  call random_number(alpha)
  alpha = alpha * 2.0 - 1.0
  call random_number(x)
  x = x * 2.0 - 1.0
  incx_val = 1
  call random_number(beta)
  beta = beta * 2.0 - 1.0
  call random_number(y)
  y = y * 2.0 - 1.0
  incy_val = 1

  ! Store original primal values
  alpha_orig = alpha
  ap_orig = ap
  x_orig = x
  beta_orig = beta
  y_orig = y

  ! Initialize output adjoints (cotangents) with random values for each direction
  ! These are the 'seeds' for reverse mode
  do k = 1, nbdirsmax
    call random_number(yb(k,:))
    yb(k,:) = yb(k,:) * 2.0 - 1.0
  end do

  ! Initialize input adjoints to zero (they will be computed)
  ! Note: Inout parameters are skipped - they already have output adjoints initialized
  alphab = 0.0
  apb = 0.0
  xb = 0.0
  betab = 0.0

  ! Save original cotangent seeds for OUTPUT/INOUT parameters (before function call)
  yb_orig = yb

  ! Set ISIZE globals required by differentiated routine (dimension 2 of arrays).
  ! Differentiated code checks they are set via check_ISIZE*_initialized.
  call set_ISIZE1OFAp(max_size)
  call set_ISIZE1OFX(max_size)

  ! Call reverse vector mode differentiated function
  call sspmv_bv(uplo, nsize, alpha, alphab, ap, apb, x, xb, incx_val, beta, betab, y, yb, incy_val, nbdirsmax)

  ! Reset ISIZE globals to uninitialized (-1) for completeness
  call set_ISIZE1OFAp(-1)
  call set_ISIZE1OFX(-1)

  ! VJP Verification using finite differences
  call check_vjp_numerically()

  write(*,*) ''
  write(*,*) 'Test completed successfully'

contains

  subroutine check_vjp_numerically()
    implicit none
    
    ! Direction vectors for VJP testing
    real(4) :: alpha_dir
    real(4), dimension((n*(n+1))/2) :: ap_dir
    real(4), dimension(max_size) :: x_dir
    real(4) :: beta_dir
    real(4), dimension(max_size) :: y_dir
    real(4), dimension(max_size) :: y_plus, y_minus, y_central_diff
    
    max_error = 0.0d0
    has_large_errors = .false.
    
    write(*,*) 'Function calls completed successfully'
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Test each differentiation direction separately
    do k = 1, nbdirsmax
      
      ! Initialize random direction vectors for all inputs
      call random_number(alpha_dir)
      alpha_dir = alpha_dir * 2.0 - 1.0
      call random_number(ap_dir)
      ap_dir = ap_dir * 2.0 - 1.0
      call random_number(x_dir)
      x_dir = x_dir * 2.0 - 1.0
      call random_number(beta_dir)
      beta_dir = beta_dir * 2.0 - 1.0
      call random_number(y_dir)
      y_dir = y_dir * 2.0 - 1.0
      
      ! Forward perturbation: f(x + h*dir)
      alpha = alpha_orig + h * alpha_dir
      ap = ap_orig + h * ap_dir
      x = x_orig + h * x_dir
      beta = beta_orig + h * beta_dir
      y = y_orig + h * y_dir
      call sspmv(uplo, nsize, alpha, ap, x, incx_val, beta, y, incy_val)
      y_plus = y
      
      ! Backward perturbation: f(x - h*dir)
      alpha = alpha_orig - h * alpha_dir
      ap = ap_orig - h * ap_dir
      x = x_orig - h * x_dir
      beta = beta_orig - h * beta_dir
      y = y_orig - h * y_dir
      call sspmv(uplo, nsize, alpha, ap, x, incx_val, beta, y, incy_val)
      y_minus = y
      
      ! Compute central differences and VJP verification
      ! VJP check: direction^T @ adjoint should equal finite difference
      
      ! Compute central differences: (f(x+h*dir) - f(x-h*dir)) / (2h)
      y_central_diff = (y_plus - y_minus) / (2.0 * h)
      
      ! VJP verification:
      ! cotangent^T @ central_diff should equal direction^T @ computed_adjoint
      ! Left side: cotangent^T @ Jacobian @ direction (via finite differences, with sorted summation)
      vjp_fd = 0.0
      ! Compute and sort products for y (FD)
      n_products = n
      do i = 1, n
        temp_products(i) = yb_orig(k,i) * y_central_diff(i)
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_fd = vjp_fd + temp_products(i)
      end do
      
      ! Right side: direction^T @ computed_adjoint (with sorted summation)
      ! For INOUT parameters: use cb directly (it contains the computed input adjoint after reverse pass)
      ! For pure inputs: use adjoint directly
      vjp_ad = 0.0
      vjp_ad = vjp_ad + alpha_dir * alphab(k)
      vjp_ad = vjp_ad + beta_dir * betab(k)
      ! Compute and sort products for ap
      n_products = (n*(n+1))/2
      do i = 1, (n*(n+1))/2
        temp_products(i) = ap_dir(i) * apb(k,i)
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_ad = vjp_ad + temp_products(i)
      end do
      ! Compute and sort products for x
      n_products = n
      do i = 1, n
        temp_products(i) = x_dir(i) * xb(k,i)
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_ad = vjp_ad + temp_products(i)
      end do
      ! Compute and sort products for y
      n_products = n
      do i = 1, n
        temp_products(i) = y_dir(i) * yb(k,i)
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

end program test_sspmv_vector_reverse