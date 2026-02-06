! Test program for ZTBMV reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Verification uses VJP methodology with finite differences

program test_ztbmv_reverse
  implicit none

  external :: ztbmv
  external :: ztbmv_b

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = 5  ! Maximum array dimension (adjusted for LD constraints)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions

  character :: uplo
  character :: trans
  character :: diag
  integer :: nsize
  integer :: ksize
  complex(8), dimension(max_size,max_size) :: a  ! Band storage (k+1) x n
  integer :: lda_val
  complex(8), dimension(max_size) :: x
  integer :: incx_val

  ! Adjoint variables (reverse mode)
  ! In reverse mode: output adjoints are INPUT (cotangents/seeds)
  !                  input adjoints are OUTPUT (computed gradients)
  complex(8), dimension(max_size,max_size) :: ab  ! Band storage
  complex(8), dimension(max_size) :: xb

  ! Storage for original values (for VJP verification)
  complex(8), dimension(max_size,max_size) :: a_orig  ! Band storage
  complex(8), dimension(max_size) :: x_orig

  ! Variables for VJP verification via finite differences
  complex(8), dimension(max_size) :: x_plus, x_minus

  ! Saved cotangents (output adjoints) for VJP verification
  complex(8), dimension(max_size) :: xb_orig
  real(8), parameter :: h = 1.0e-7
  real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
  logical :: has_large_errors
  integer :: i, j, band_row
  real(4) :: temp_real, temp_imag  ! For band matrix initialization
  real(8), dimension(max_size*max_size) :: temp_products  ! For sorted summation
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
  ksize = max(0, n - 1)  ! Band width: 0 <= K <= N-1
  ! Initialize a as triangular band matrix (upper band storage)
  do j = 1, n
    do band_row = max(1, ksize+2-j), ksize+1
      call random_number(temp_real)
      call random_number(temp_imag)
      a(band_row, j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
  end do
  lda_val = lda
  do i = 1, max_size
    call random_number(temp_real_init)
    call random_number(temp_imag_init)
    x(i) = cmplx(temp_real_init, temp_imag_init) * (2.0,2.0) - (1.0,1.0)
  end do
  incx_val = 1

  ! Store original primal values
  a_orig = a
  x_orig = x

  write(*,*) 'Testing ZTBMV'

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
  ab = 0.0d0

  ! Set ISIZE globals required by differentiated routine (dimension 2 of arrays).
  ! Differentiated code checks they are set via check_ISIZE*_initialized.
  call set_ISIZE2OFA(max_size)

  ! Call reverse mode differentiated function
  call ztbmv_b(uplo, trans, diag, nsize, ksize, a, ab, lda_val, x, xb, incx_val)

  ! Reset ISIZE globals to uninitialized (-1) for completeness
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
    
    integer :: band_row  ! Loop variable for band storage
    ! Temporary variables for complex random number generation
    real(4) :: temp_real, temp_imag
    
    ! Direction vectors for VJP testing (like tangents in forward mode)
    complex(8), dimension(max_size,max_size) :: a_dir  ! Band storage
    complex(8), dimension(max_size) :: x_dir
    
    complex(8), dimension(max_size) :: x_central_diff
    
    max_error = 0.0d0
    has_large_errors = .false.
    
    write(*,*) 'Function calls completed successfully'
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Initialize random direction vectors for all inputs
      ! Keep direction consistent with triangular band: only band entries used
      do j = 1, n
        do band_row = max(1, ksize+2-j), ksize+1
          call random_number(temp_real)
          call random_number(temp_imag)
          a_dir(band_row, j) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
        end do
      end do
    do i = 1, max_size
      call random_number(temp_real)
      call random_number(temp_imag)
      x_dir(i) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
    
    ! Forward perturbation: f(x + h*dir)
    a = a_orig + cmplx(h, 0.0) * a_dir
    x = x_orig + cmplx(h, 0.0) * x_dir
    call ztbmv(uplo, trans, diag, nsize, ksize, a, lda_val, x, incx_val)
    x_plus = x
    
    ! Backward perturbation: f(x - h*dir)
    a = a_orig - cmplx(h, 0.0) * a_dir
    x = x_orig - cmplx(h, 0.0) * x_dir
    call ztbmv(uplo, trans, diag, nsize, ksize, a, lda_val, x, incx_val)
    x_minus = x
    
    ! Compute central differences: (f(x+h*dir) - f(x-h*dir)) / (2h)
    x_central_diff = (x_plus - x_minus) / (2.0d0 * h)
    
    ! VJP verification:
    ! cotangent^T @ central_diff should equal direction^T @ computed_adjoint
    ! Left side: cotangent^T @ Jacobian @ direction (via finite differences, with sorted summation)
    vjp_fd = 0.0d0
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
    vjp_ad = 0.0d0
    ! Compute and sort products for a (band storage)
    n_products = 0
    do j = 1, n
      do band_row = max(1, ksize+2-j), ksize+1
        n_products = n_products + 1
        temp_products(n_products) = real(conjg(a_dir(band_row,j)) * ab(band_row,j))
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

end program test_ztbmv_reverse