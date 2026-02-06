! Test program for CTPMV reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Verification uses VJP methodology with finite differences

program test_ctpmv_reverse
  implicit none

  external :: ctpmv
  external :: ctpmv_b

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension (rows/cols of matrices)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions

  character :: uplo
  character :: trans
  character :: diag
  integer :: nsize
  complex(4), dimension((n*(n+1))/2) :: ap
  complex(4), dimension(max_size) :: x
  integer :: incx_val

  ! Adjoint variables (reverse mode)
  ! In reverse mode: output adjoints are INPUT (cotangents/seeds)
  !                  input adjoints are OUTPUT (computed gradients)
  complex(4), dimension((n*(n+1))/2) :: apb
  complex(4), dimension(max_size) :: xb

  ! Storage for original values (for VJP verification)
  complex(4), dimension((n*(n+1))/2) :: ap_orig
  complex(4), dimension(max_size) :: x_orig

  ! Variables for VJP verification via finite differences
  complex(4), dimension(max_size) :: x_plus, x_minus

  ! Saved cotangents (output adjoints) for VJP verification
  complex(4), dimension(max_size) :: xb_orig
  real(4), parameter :: h = 1.0e-3
  real(4) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
  logical :: has_large_errors
  integer :: i, j
  real(4), dimension(max_size*max_size) :: temp_products  ! For sorted summation
  integer :: n_products

  ! Temporary variables for complex random initialization
  real(4) :: temp_real_init, temp_imag_init

  ! Initialize random seed for reproducibility
  integer :: seed_array(33)
  seed_array = 42
  call random_seed(put=seed_array)

  ! Initialize primal values
  uplo = 'U'
  trans = 'N'
  diag = 'N'
  nsize = n
  do i = 1, (n*(n+1))/2
    call random_number(temp_real_init)
    call random_number(temp_imag_init)
    ap(i) = cmplx(temp_real_init, temp_imag_init) * (2.0,2.0) - (1.0,1.0)
  end do
  do i = 1, max_size
    call random_number(temp_real_init)
    call random_number(temp_imag_init)
    x(i) = cmplx(temp_real_init, temp_imag_init) * (2.0,2.0) - (1.0,1.0)
  end do
  incx_val = 1

  ! Store original primal values
  ap_orig = ap
  x_orig = x

  write(*,*) 'Testing CTPMV'

  ! Initialize output adjoints (cotangents) with random values
  ! These are the 'seeds' for reverse mode
  do i = 1, max_size
    call random_number(temp_real_init)
    call random_number(temp_imag_init)
    xb(i) = cmplx(temp_real_init, temp_imag_init) * (2.0,2.0) - (1.0,1.0)
  end do

  ! Save output adjoints (cotangents) for VJP verification
  ! Note: output adjoints may be modified by reverse mode function
  xb_orig = xb

  ! Initialize input adjoints to zero (they will be computed)
  apb = 0.0

  ! Set ISIZE globals required by differentiated routine (dimension 2 of arrays).
  ! Differentiated code checks they are set via check_ISIZE*_initialized.
  call set_ISIZE1OFAp(max_size)

  ! Call reverse mode differentiated function
  call ctpmv_b(uplo, trans, diag, nsize, ap, apb, x, xb, incx_val)

  ! Reset ISIZE globals to uninitialized (-1) for completeness
  call set_ISIZE1OFAp(-1)

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
    complex(4), dimension(max_size*(max_size+1)/2) :: ap_dir
    complex(4), dimension(max_size) :: x_dir
    
    complex(4), dimension(max_size) :: x_central_diff
    
    max_error = 0.0
    has_large_errors = .false.
    
    write(*,*) 'Function calls completed successfully'
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Initialize random direction vectors for all inputs
    do i = 1, max_size*(max_size+1)/2
      call random_number(temp_real)
      call random_number(temp_imag)
      ap_dir(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
    do i = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      x_dir(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
    
    ! Forward perturbation: f(x + h*dir)
    ap = ap_orig + cmplx(h, 0.0) * ap_dir
    x = x_orig + cmplx(h, 0.0) * x_dir
    call ctpmv(uplo, trans, diag, nsize, ap, x, incx_val)
    x_plus = x
    
    ! Backward perturbation: f(x - h*dir)
    ap = ap_orig - cmplx(h, 0.0) * ap_dir
    x = x_orig - cmplx(h, 0.0) * x_dir
    call ctpmv(uplo, trans, diag, nsize, ap, x, incx_val)
    x_minus = x
    
    ! Compute central differences: (f(x+h*dir) - f(x-h*dir)) / (2h)
    x_central_diff = (x_plus - x_minus) / (2.0d0 * h)
    
    ! VJP verification:
    ! cotangent^T @ central_diff should equal direction^T @ computed_adjoint
    ! Left side: cotangent^T @ Jacobian @ direction (via finite differences, with sorted summation)
    vjp_fd = 0.0
    ! Compute and sort products for x (FD)
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(xb_orig(i)) * x_central_diff(i))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do
    
    ! Right side: direction^T @ computed_adjoint (with sorted summation)
    ! For INOUT parameters: use cb directly (it contains the computed input adjoint after reverse pass)
    ! For pure inputs: use adjoint directly
    vjp_ad = 0.0
    ! Compute and sort products for ap
    n_products = n*(n+1)/2
    do i = 1, n_products
      temp_products(i) = real(conjg(ap_dir(i)) * apb(i))
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
    
    ! Error check: |vjp_fd - vjp_ad| > atol + rtol * |vjp_ad|
    abs_error = abs(vjp_fd - vjp_ad)
    abs_reference = abs(vjp_ad)
    error_bound = 1.0e-3 + 1.0e-3 * abs_reference
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

end program test_ctpmv_reverse