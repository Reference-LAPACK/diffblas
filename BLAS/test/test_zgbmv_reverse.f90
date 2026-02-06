! Test program for ZGBMV reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Verification uses VJP methodology with finite differences

program test_zgbmv_reverse
  implicit none

  external :: zgbmv
  external :: zgbmv_b

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension (rows/cols of matrices)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions

  character :: trans
  integer :: msize
  integer :: nsize
  integer :: kl
  integer :: ku
  complex(8) :: alpha
  complex(8), dimension(max_size,max_size) :: a
  integer :: lda_val
  complex(8), dimension(max_size) :: x
  integer :: incx_val
  complex(8) :: beta
  complex(8), dimension(max_size) :: y
  integer :: incy_val

  ! Adjoint variables (reverse mode)
  ! In reverse mode: output adjoints are INPUT (cotangents/seeds)
  !                  input adjoints are OUTPUT (computed gradients)
  complex(8) :: alphab
  complex(8), dimension(max_size,max_size) :: ab
  complex(8), dimension(max_size) :: xb
  complex(8) :: betab
  complex(8), dimension(max_size) :: yb

  ! Storage for original values (for VJP verification)
  complex(8) :: alpha_orig
  complex(8), dimension(max_size,max_size) :: a_orig
  complex(8), dimension(max_size) :: x_orig
  complex(8) :: beta_orig
  complex(8), dimension(max_size) :: y_orig

  ! Variables for VJP verification via finite differences
  complex(8), dimension(max_size) :: y_plus, y_minus

  ! Saved cotangents (output adjoints) for VJP verification
  complex(8), dimension(max_size) :: yb_orig
  real(8), parameter :: h = 1.0e-7
  real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
  logical :: has_large_errors
  integer :: i, j
  real(8), dimension(max_size*max_size) :: temp_products  ! For sorted summation
  integer :: n_products

  ! Temporary variables for complex random initialization
  real(4) :: temp_real_init, temp_imag_init

  ! Initialize random seed for reproducibility
  integer :: seed_array(33)
  seed_array = 42
  call random_seed(put=seed_array)

  ! Initialize primal values
  trans = 'N'
  msize = n
  nsize = n
  kl = 1
  ku = 1
  call random_number(temp_real_init)
  call random_number(temp_imag_init)
  alpha = cmplx(temp_real_init, temp_imag_init) * (2.0,2.0) - (1.0,1.0)
  do j = 1, max_size
    do i = 1, max_size
      call random_number(temp_real_init)
      call random_number(temp_imag_init)
      a(i,j) = cmplx(temp_real_init, temp_imag_init) * (2.0,2.0) - (1.0,1.0)
    end do
  end do
  lda_val = lda
  do i = 1, max_size
    call random_number(temp_real_init)
    call random_number(temp_imag_init)
    x(i) = cmplx(temp_real_init, temp_imag_init) * (2.0,2.0) - (1.0,1.0)
  end do
  incx_val = 1
  call random_number(temp_real_init)
  call random_number(temp_imag_init)
  beta = cmplx(temp_real_init, temp_imag_init) * (2.0,2.0) - (1.0,1.0)
  do i = 1, max_size
    call random_number(temp_real_init)
    call random_number(temp_imag_init)
    y(i) = cmplx(temp_real_init, temp_imag_init) * (2.0,2.0) - (1.0,1.0)
  end do
  incy_val = 1

  ! Store original primal values
  alpha_orig = alpha
  a_orig = a
  x_orig = x
  beta_orig = beta
  y_orig = y

  write(*,*) 'Testing ZGBMV'

  ! Initialize output adjoints (cotangents) with random values
  ! These are the 'seeds' for reverse mode
  do i = 1, max_size
    call random_number(temp_real_init)
    call random_number(temp_imag_init)
    yb(i) = cmplx(temp_real_init, temp_imag_init) * (2.0,2.0) - (1.0,1.0)
  end do

  ! Save output adjoints (cotangents) for VJP verification
  ! Note: output adjoints may be modified by reverse mode function
  yb_orig = yb

  ! Initialize input adjoints to zero (they will be computed)
  alphab = 0.0d0
  betab = 0.0d0
  ab = 0.0d0
  xb = 0.0d0

  ! Set ISIZE globals required by differentiated routine (dimension 2 of arrays).
  ! Differentiated code checks they are set via check_ISIZE*_initialized.
  call set_ISIZE1OFX(max_size)
  call set_ISIZE2OFA(max_size)

  ! Call reverse mode differentiated function
  call zgbmv_b(trans, msize, nsize, kl, ku, alpha, alphab, a, ab, lda_val, x, xb, incx_val, beta, betab, y, yb, incy_val)

  ! Reset ISIZE globals to uninitialized (-1) for completeness
  call set_ISIZE1OFX(-1)
  call set_ISIZE2OFA(-1)

  ! VJP Verification using finite differences
  ! For reverse mode, we verify: cotangent^T @ J @ direction = direction^T @ adjoint
  ! Equivalently: cotangent^T @ (f(x+h*dir) - f(x-h*dir))/(2h) should equal dir^T @ computed_adjoint
  call check_vjp_numerically()

  write(*,*) ''
  write(*,*) 'Test completed successfully'

contains

  subroutine check_vjp_numerically()
    implicit none
    
    ! Temporary variables for complex random number generation
    real(4) :: temp_real, temp_imag
    
    ! Direction vectors for VJP testing (like tangents in forward mode)
    complex(8) :: alpha_dir
    complex(8), dimension(max_size,max_size) :: a_dir
    complex(8), dimension(max_size) :: x_dir
    complex(8) :: beta_dir
    complex(8), dimension(max_size) :: y_dir
    
    complex(8), dimension(max_size) :: y_central_diff
    
    max_error = 0.0d0
    has_large_errors = .false.
    
    write(*,*) 'Function calls completed successfully'
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Initialize random direction vectors for all inputs
    call random_number(temp_real)
    call random_number(temp_imag)
    alpha_dir = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    do j = 1, max_size
      do i = 1, max_size
        call random_number(temp_real)
        call random_number(temp_imag)
        a_dir(i,j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
      end do
    end do
    do i = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      x_dir(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
    call random_number(temp_real)
    call random_number(temp_imag)
    beta_dir = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    do i = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      y_dir(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
    
    ! Forward perturbation: f(x + h*dir)
    alpha = alpha_orig + cmplx(h, 0.0) * alpha_dir
    a = a_orig + cmplx(h, 0.0) * a_dir
    x = x_orig + cmplx(h, 0.0) * x_dir
    beta = beta_orig + cmplx(h, 0.0) * beta_dir
    y = y_orig + cmplx(h, 0.0) * y_dir
    call zgbmv(trans, msize, nsize, kl, ku, alpha, a, lda_val, x, incx_val, beta, y, incy_val)
    y_plus = y
    
    ! Backward perturbation: f(x - h*dir)
    alpha = alpha_orig - cmplx(h, 0.0) * alpha_dir
    a = a_orig - cmplx(h, 0.0) * a_dir
    x = x_orig - cmplx(h, 0.0) * x_dir
    beta = beta_orig - cmplx(h, 0.0) * beta_dir
    y = y_orig - cmplx(h, 0.0) * y_dir
    call zgbmv(trans, msize, nsize, kl, ku, alpha, a, lda_val, x, incx_val, beta, y, incy_val)
    y_minus = y
    
    ! Compute central differences: (f(x+h*dir) - f(x-h*dir)) / (2h)
    y_central_diff = (y_plus - y_minus) / (2.0d0 * h)
    
    ! VJP verification:
    ! cotangent^T @ central_diff should equal direction^T @ computed_adjoint
    ! Left side: cotangent^T @ Jacobian @ direction (via finite differences, with sorted summation)
    vjp_fd = 0.0d0
    ! Compute and sort products for y (FD)
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(yb_orig(i)) * y_central_diff(i))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do
    
    ! Right side: direction^T @ computed_adjoint (with sorted summation)
    ! For INOUT parameters: use cb directly (it contains the computed input adjoint after reverse pass)
    ! For pure inputs: use adjoint directly
    vjp_ad = 0.0d0
    vjp_ad = vjp_ad + real(conjg(alpha_dir) * alphab)
    ! Compute and sort products for a
    n_products = 0
    do j = 1, n
      do i = 1, n
        n_products = n_products + 1
        temp_products(n_products) = real(conjg(a_dir(i,j)) * ab(i,j))
      end do
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    ! Compute and sort products for x
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(x_dir(i)) * xb(i))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    vjp_ad = vjp_ad + real(conjg(beta_dir) * betab)
    ! Compute and sort products for y
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(y_dir(i)) * yb(i))
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

end program test_zgbmv_reverse