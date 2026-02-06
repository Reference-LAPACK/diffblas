! Test program for DSYR2 reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Verification uses VJP methodology with finite differences

program test_dsyr2_reverse
  implicit none

  external :: dsyr2
  external :: dsyr2_b

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension (rows/cols of matrices)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions

  character :: uplo
  integer :: nsize
  real(8) :: alpha
  real(8), dimension(max_size) :: x
  integer :: incx_val
  real(8), dimension(max_size) :: y
  integer :: incy_val
  real(8), dimension(max_size,max_size) :: a
  integer :: lda_val

  ! Adjoint variables (reverse mode)
  ! In reverse mode: output adjoints are INPUT (cotangents/seeds)
  !                  input adjoints are OUTPUT (computed gradients)
  real(8) :: alphab
  real(8), dimension(max_size) :: xb
  real(8), dimension(max_size) :: yb
  real(8), dimension(max_size,max_size) :: ab

  ! Storage for original values (for VJP verification)
  real(8) :: alpha_orig
  real(8), dimension(max_size) :: x_orig
  real(8), dimension(max_size) :: y_orig
  real(8), dimension(max_size,max_size) :: a_orig

  ! Variables for VJP verification via finite differences
  real(8), dimension(max_size,max_size) :: a_plus, a_minus

  ! Saved cotangents (output adjoints) for VJP verification
  real(8), dimension(max_size,max_size) :: ab_orig
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
  uplo = 'U'
  nsize = n
  call random_number(alpha)
  alpha = alpha * 2.0d0 - 1.0d0
  call random_number(x)
  x = x * 2.0d0 - 1.0d0
  incx_val = 1
  call random_number(y)
  y = y * 2.0d0 - 1.0d0
  incy_val = 1
  call random_number(a)
  a = a * 2.0d0 - 1.0d0
  lda_val = lda

  ! Store original primal values
  alpha_orig = alpha
  x_orig = x
  y_orig = y
  a_orig = a

  write(*,*) 'Testing DSYR2'

  ! Initialize output adjoints (cotangents) with random values
  ! These are the 'seeds' for reverse mode
  call random_number(ab)
  ab = ab * 2.0d0 - 1.0d0

  ! Save output adjoints (cotangents) for VJP verification
  ! Note: output adjoints may be modified by reverse mode function
  ab_orig = ab

  ! Initialize input adjoints to zero (they will be computed)
  alphab = 0.0d0
  xb = 0.0d0
  yb = 0.0d0

  ! Set ISIZE globals required by differentiated routine (dimension 2 of arrays).
  ! Differentiated code checks they are set via check_ISIZE*_initialized.
  call set_ISIZE1OFX(max_size)
  call set_ISIZE1OFY(max_size)

  ! Call reverse mode differentiated function
  call dsyr2_b(uplo, nsize, alpha, alphab, x, xb, incx_val, y, yb, incy_val, a, ab, lda_val)

  ! Reset ISIZE globals to uninitialized (-1) for completeness
  call set_ISIZE1OFX(-1)
  call set_ISIZE1OFY(-1)

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
    real(8) :: alpha_dir
    real(8), dimension(max_size) :: x_dir
    real(8), dimension(max_size) :: y_dir
    real(8), dimension(max_size,max_size) :: a_dir
    
    real(8), dimension(max_size,max_size) :: a_central_diff
    
    max_error = 0.0d0
    has_large_errors = .false.
    
    write(*,*) 'Function calls completed successfully'
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Initialize random direction vectors for all inputs
    call random_number(alpha_dir)
    alpha_dir = alpha_dir * 2.0d0 - 1.0d0
    call random_number(x_dir)
    x_dir = x_dir * 2.0d0 - 1.0d0
    call random_number(y_dir)
    y_dir = y_dir * 2.0d0 - 1.0d0
    call random_number(a_dir)
    a_dir = a_dir * 2.0d0 - 1.0d0
    
    ! Forward perturbation: f(x + h*dir)
    alpha = alpha_orig + h * alpha_dir
    x = x_orig + h * x_dir
    y = y_orig + h * y_dir
    a = a_orig + h * a_dir
    call dsyr2(uplo, nsize, alpha, x, incx_val, y, incy_val, a, lda_val)
    a_plus = a
    
    ! Backward perturbation: f(x - h*dir)
    alpha = alpha_orig - h * alpha_dir
    x = x_orig - h * x_dir
    y = y_orig - h * y_dir
    a = a_orig - h * a_dir
    call dsyr2(uplo, nsize, alpha, x, incx_val, y, incy_val, a, lda_val)
    a_minus = a
    
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
        temp_products(n_products) = ab_orig(i,j) * a_central_diff(i,j)
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
    vjp_ad = vjp_ad + alpha_dir * alphab
    ! Compute and sort products for x
    n_products = n
    do i = 1, n
      temp_products(i) = x_dir(i) * xb(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    ! Compute and sort products for y
    n_products = n
    do i = 1, n
      temp_products(i) = y_dir(i) * yb(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    ! Compute and sort products for a
    n_products = 0
    do j = 1, n
      do i = 1, n
        n_products = n_products + 1
        temp_products(n_products) = a_dir(i,j) * ab(i,j)
      end do
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

end program test_dsyr2_reverse