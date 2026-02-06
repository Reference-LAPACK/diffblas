! Test program for ZCOPY reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Verification uses VJP methodology with finite differences

program test_zcopy_reverse
  implicit none

  external :: zcopy
  external :: zcopy_b

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension (rows/cols of matrices)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions

  integer :: nsize
  complex(8), dimension(max_size) :: zx
  integer :: incx_val
  complex(8), dimension(max_size) :: zy
  integer :: incy_val

  ! Adjoint variables (reverse mode)
  ! In reverse mode: output adjoints are INPUT (cotangents/seeds)
  !                  input adjoints are OUTPUT (computed gradients)
  complex(8), dimension(max_size) :: zxb
  complex(8), dimension(max_size) :: zyb

  ! Storage for original values (for VJP verification)
  complex(8), dimension(max_size) :: zx_orig
  complex(8), dimension(max_size) :: zy_orig

  ! Variables for VJP verification via finite differences
  complex(8), dimension(max_size) :: zy_plus, zy_minus

  ! Saved cotangents (output adjoints) for VJP verification
  complex(8), dimension(max_size) :: zyb_orig
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
  nsize = n
  do i = 1, max_size
    call random_number(temp_real_init)
    call random_number(temp_imag_init)
    zx(i) = cmplx(temp_real_init, temp_imag_init) * (2.0,2.0) - (1.0,1.0)
  end do
  incx_val = 1
  do i = 1, max_size
    call random_number(temp_real_init)
    call random_number(temp_imag_init)
    zy(i) = cmplx(temp_real_init, temp_imag_init) * (2.0,2.0) - (1.0,1.0)
  end do
  incy_val = 1

  ! Store original primal values
  zx_orig = zx
  zy_orig = zy

  write(*,*) 'Testing ZCOPY'

  ! Initialize output adjoints (cotangents) with random values
  ! These are the 'seeds' for reverse mode
  do i = 1, max_size
    call random_number(temp_real_init)
    call random_number(temp_imag_init)
    zyb(i) = cmplx(temp_real_init, temp_imag_init) * (2.0,2.0) - (1.0,1.0)
  end do

  ! Save output adjoints (cotangents) for VJP verification
  ! Note: output adjoints may be modified by reverse mode function
  zyb_orig = zyb

  ! Initialize input adjoints to zero (they will be computed)
  zxb = 0.0d0

  ! Set ISIZE globals required by differentiated routine (dimension 2 of arrays).
  ! Differentiated code checks they are set via check_ISIZE*_initialized.
  call set_ISIZE1OFZx(max_size)

  ! Call reverse mode differentiated function
  call zcopy_b(nsize, zx, zxb, incx_val, zy, zyb, incy_val)

  ! Reset ISIZE globals to uninitialized (-1) for completeness
  call set_ISIZE1OFZx(-1)

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
    complex(8), dimension(max_size) :: zx_dir
    complex(8), dimension(max_size) :: zy_dir
    
    complex(8), dimension(max_size) :: zy_central_diff
    
    max_error = 0.0d0
    has_large_errors = .false.
    
    write(*,*) 'Function calls completed successfully'
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Initialize random direction vectors for all inputs
    do i = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      zx_dir(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
    do i = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      zy_dir(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
    
    ! Forward perturbation: f(x + h*dir)
    zx = zx_orig + cmplx(h, 0.0) * zx_dir
    zy = zy_orig + cmplx(h, 0.0) * zy_dir
    call zcopy(nsize, zx, incx_val, zy, incy_val)
    zy_plus = zy
    
    ! Backward perturbation: f(x - h*dir)
    zx = zx_orig - cmplx(h, 0.0) * zx_dir
    zy = zy_orig - cmplx(h, 0.0) * zy_dir
    call zcopy(nsize, zx, incx_val, zy, incy_val)
    zy_minus = zy
    
    ! Compute central differences: (f(x+h*dir) - f(x-h*dir)) / (2h)
    zy_central_diff = (zy_plus - zy_minus) / (2.0d0 * h)
    
    ! VJP verification:
    ! cotangent^T @ central_diff should equal direction^T @ computed_adjoint
    ! Left side: cotangent^T @ Jacobian @ direction (via finite differences, with sorted summation)
    vjp_fd = 0.0d0
    ! Compute and sort products for zy (FD)
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(zyb_orig(i)) * zy_central_diff(i))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do
    
    ! Right side: direction^T @ computed_adjoint (with sorted summation)
    ! For INOUT parameters: use cb directly (it contains the computed input adjoint after reverse pass)
    ! For pure inputs: use adjoint directly
    vjp_ad = 0.0d0
    ! Compute and sort products for zx
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(zx_dir(i)) * zxb(i))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    ! Compute and sort products for zy
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(zy_dir(i)) * zyb(i))
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

end program test_zcopy_reverse