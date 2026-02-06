! Test program for DTRMM vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirsmax=4

program test_dtrmm_vector_reverse
  implicit none
  include 'DIFFSIZES.inc'

  external :: dtrmm
  external :: dtrmm_bv

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, k  ! Loop counters
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  character :: side
  character :: uplo
  character :: transa
  character :: diag
  integer :: msize
  integer :: nsize
  real(8) :: alpha
  real(8), dimension(max_size,max_size) :: a
  integer :: lda_val
  real(8), dimension(max_size,max_size) :: b
  integer :: ldb_val

  ! Adjoint variables (reverse vector mode)
  ! In reverse mode: output adjoints are INPUT (cotangents/seeds)
  !                  input adjoints are OUTPUT (computed gradients)
  real(8), dimension(nbdirsmax) :: alphab
  real(8), dimension(nbdirsmax,max_size,max_size) :: ab
  real(8), dimension(nbdirsmax,max_size,max_size) :: bb

  ! Storage for original cotangents (for INOUT parameters in VJP verification)
  real(8), dimension(nbdirsmax,max_size,max_size) :: bb_orig

  ! Storage for original values (for VJP verification)
  real(8) :: alpha_orig
  real(8), dimension(max_size,max_size) :: a_orig
  real(8), dimension(max_size,max_size) :: b_orig

  ! Variables for VJP verification via finite differences
  real(8), parameter :: h = 1.0e-7
  real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
  logical :: has_large_errors
  real(8), dimension(max_size*max_size) :: temp_products  ! For sorted summation
  integer :: n_products

  ! Initialize random seed for reproducibility
  seed_array = 42
  call random_seed(put=seed_array)

  ! Initialize primal values
  side = 'L'
  uplo = 'U'
  transa = 'N'
  diag = 'N'
  msize = n
  nsize = n
  call random_number(alpha)
  alpha = alpha * 2.0 - 1.0
  call random_number(a)
  a = a * 2.0 - 1.0
  lda_val = lda
  call random_number(b)
  b = b * 2.0 - 1.0
  ldb_val = ldb

  ! Store original primal values
  alpha_orig = alpha
  a_orig = a
  b_orig = b

  ! Initialize output adjoints (cotangents) with random values for each direction
  ! These are the 'seeds' for reverse mode
  do k = 1, nbdirsmax
    call random_number(bb(k,:,:))
    bb(k,:,:) = bb(k,:,:) * 2.0 - 1.0
  end do

  ! Initialize input adjoints to zero (they will be computed)
  ! Note: Inout parameters are skipped - they already have output adjoints initialized
  alphab = 0.0
  ab = 0.0

  ! Save original cotangent seeds for OUTPUT/INOUT parameters (before function call)
  bb_orig = bb

  ! Set ISIZE globals required by differentiated routine (dimension 2 of arrays).
  ! Differentiated code checks they are set via check_ISIZE*_initialized.
  call set_ISIZE2OFA(max_size)

  ! Call reverse vector mode differentiated function
  call dtrmm_bv(side, uplo, transa, diag, msize, nsize, alpha, alphab, a, ab, lda_val, b, bb, ldb_val, nbdirsmax)

  ! Reset ISIZE globals to uninitialized (-1) for completeness
  call set_ISIZE2OFA(-1)

  ! VJP Verification using finite differences
  call check_vjp_numerically()

  write(*,*) ''
  write(*,*) 'Test completed successfully'

contains

  subroutine check_vjp_numerically()
    implicit none
    
    ! Direction vectors for VJP testing
    real(8) :: alpha_dir
    real(8), dimension(max_size,max_size) :: a_dir
    real(8), dimension(max_size,max_size) :: b_dir
    real(8), dimension(max_size,max_size) :: b_plus, b_minus, b_central_diff
    
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
      call random_number(a_dir)
      a_dir = a_dir * 2.0 - 1.0
      call random_number(b_dir)
      b_dir = b_dir * 2.0 - 1.0
      
      ! Forward perturbation: f(x + h*dir)
      alpha = alpha_orig + h * alpha_dir
      a = a_orig + h * a_dir
      b = b_orig + h * b_dir
      call dtrmm(side, uplo, transa, diag, msize, nsize, alpha, a, lda_val, b, ldb_val)
      b_plus = b
      
      ! Backward perturbation: f(x - h*dir)
      alpha = alpha_orig - h * alpha_dir
      a = a_orig - h * a_dir
      b = b_orig - h * b_dir
      call dtrmm(side, uplo, transa, diag, msize, nsize, alpha, a, lda_val, b, ldb_val)
      b_minus = b
      
      ! Compute central differences and VJP verification
      ! VJP check: direction^T @ adjoint should equal finite difference
      
      ! Compute central differences: (f(x+h*dir) - f(x-h*dir)) / (2h)
      b_central_diff = (b_plus - b_minus) / (2.0d0 * h)
      
      ! VJP verification:
      ! cotangent^T @ central_diff should equal direction^T @ computed_adjoint
      ! Left side: cotangent^T @ Jacobian @ direction (via finite differences, with sorted summation)
      vjp_fd = 0.0d0
      ! Compute and sort products for b (FD)
      n_products = 0
      do j = 1, n
        do i = 1, n
          n_products = n_products + 1
          temp_products(n_products) = bb_orig(k,i,j) * b_central_diff(i,j)
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
      vjp_ad = vjp_ad + alpha_dir * alphab(k)
      ! Compute and sort products for a
      n_products = 0
      do j = 1, n
        do i = 1, n
          n_products = n_products + 1
          temp_products(n_products) = a_dir(i,j) * ab(k,i,j)
        end do
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_ad = vjp_ad + temp_products(i)
      end do
      ! Compute and sort products for b
      n_products = 0
      do j = 1, n
        do i = 1, n
          n_products = n_products + 1
          temp_products(n_products) = b_dir(i,j) * bb(k,i,j)
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

end program test_dtrmm_vector_reverse