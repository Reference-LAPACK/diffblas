! Test program for DTRMV vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirs=4

program test_dtrmv_vector_reverse
  implicit none
  integer, parameter :: nbdirs = 4

  external :: dtrmv
  external :: dtrmv_bv

  ! Test parameters
  integer :: n  ! Current size (set in loop)
  integer, parameter :: max_size = 100  ! Maximum array dimension (multi-size: 1,4,40,100)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, k  ! Loop counters
  integer :: test_sizes(1), itest
  logical :: passed, all_passed
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  character :: uplo
  character :: trans
  character :: diag
  integer :: nsize
  real(8), dimension(max_size,max_size) :: a
  integer :: lda_val
  real(8), dimension(max_size) :: x
  integer :: incx_val

  ! Adjoint variables (reverse vector mode)
  ! In reverse mode: output adjoints are INPUT (cotangents/seeds)
  !                  input adjoints are OUTPUT (computed gradients)
  real(8), dimension(nbdirs,max_size,max_size) :: ab
  real(8), dimension(nbdirs,max_size) :: xb

  ! Storage for original cotangents (for INOUT parameters in VJP verification)
  real(8), dimension(nbdirs,max_size) :: xb_orig

  ! Storage for original values (for VJP verification)
  real(8), dimension(max_size,max_size) :: a_orig
  real(8), dimension(max_size) :: x_orig

  ! Variables for VJP verification via finite differences
  real(8), parameter :: h = 1.0e-7
  real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
  logical :: has_large_errors
  real(8), dimension(max_size*max_size) :: temp_products  ! For sorted summation
  integer :: n_products

  ! Initialize random seed for reproducibility
  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing DTRMV (Vector Reverse, multi-size: n = 4)'
  all_passed = .true.
  do itest = 1, 1
    n = test_sizes(itest)
    write(*,*) 'Testing DTRMV (Vector Reverse, n =', n, ')'

    call run_test_for_size(n, passed)
  all_passed = all_passed .and. passed
  end do
  if (all_passed) then
    write(*,*) 'PASS: Vector reverse mode - all sizes completed successfully'
  else
    write(*,*) 'FAIL: Vector reverse mode - one or more sizes had derivative errors'
  end if

contains

  subroutine run_test_for_size(n, passed)
    implicit none
    integer, intent(in) :: n
    logical, intent(out) :: passed

    ! Initialize primal values
    uplo = 'U'
    trans = 'N'
    diag = 'N'
    nsize = n
    call random_number(a)
    a = a * 2.0 - 1.0
    lda_val = lda
    call random_number(x)
    x = x * 2.0 - 1.0
    incx_val = 1
    
    ! Store original primal values
    a_orig = a
    x_orig = x
    
    ! Initialize output adjoints (cotangents) with random values for each direction
    ! These are the 'seeds' for reverse mode
    do k = 1, nbdirs
      call random_number(xb(k,:))
      xb(k,:) = xb(k,:) * 2.0 - 1.0
    end do
    
    ! Initialize input adjoints to zero (they will be computed)
    ! Note: Inout parameters are skipped - they already have output adjoints initialized
    ab = 0.0
    
    ! Save original cotangent seeds for OUTPUT/INOUT parameters (before function call)
    xb_orig = xb
    
    ! Set ISIZE globals required by differentiated routine (dimension 2 of arrays).
    ! ISIZE1OF* (vectors): use n to match adjoint array size; ISIZE2OF* (matrices): use max_size.
    call set_ISIZE2OFA(max_size)
    
    ! Call reverse vector mode differentiated function
    call dtrmv_bv(uplo, trans, diag, nsize, a, ab, lda_val, x, xb, incx_val, nbdirs)
    
    ! Reset ISIZE globals to uninitialized (-1) for completeness
    call set_ISIZE2OFA(-1)
    
    ! VJP Verification using finite differences
    call check_vjp_numerically(passed)
  end subroutine run_test_for_size

  subroutine check_vjp_numerically(passed)
    implicit none
    logical, intent(out) :: passed
    
    ! Direction vectors for VJP testing
    real(8), dimension(max_size,max_size) :: a_dir
    real(8), dimension(max_size) :: x_dir
    real(8), dimension(max_size) :: x_plus, x_minus, x_central_diff
    
    max_error = 0.0d0
    has_large_errors = .false.
    
    write(*,*) 'Function calls completed successfully'
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Test each differentiation direction separately
    do k = 1, nbdirs
      
      ! Initialize random direction vectors for all inputs
      call random_number(a_dir)
      a_dir = a_dir * 2.0 - 1.0
      call random_number(x_dir)
      x_dir = x_dir * 2.0 - 1.0
      
      ! Forward perturbation: f(x + h*dir)
      a = a_orig + h * a_dir
      x = x_orig + h * x_dir
      call dtrmv(uplo, trans, diag, nsize, a, lda_val, x, incx_val)
      x_plus = x
      
      ! Backward perturbation: f(x - h*dir)
      a = a_orig - h * a_dir
      x = x_orig - h * x_dir
      call dtrmv(uplo, trans, diag, nsize, a, lda_val, x, incx_val)
      x_minus = x
      
      ! Compute central differences and VJP verification
      ! VJP check: direction^T @ adjoint should equal finite difference
      
      ! Compute central differences: (f(x+h*dir) - f(x-h*dir)) / (2h)
      x_central_diff = (x_plus - x_minus) / (2.0d0 * h)
      
      ! VJP verification:
      ! cotangent^T @ central_diff should equal direction^T @ computed_adjoint
      ! Left side: cotangent^T @ Jacobian @ direction (via finite differences, with sorted summation)
      vjp_fd = 0.0d0
      ! Compute and sort products for x (FD)
      n_products = n
      do i = 1, n
        temp_products(i) = xb_orig(k,i) * x_central_diff(i)
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_fd = vjp_fd + temp_products(i)
      end do
      
      ! Right side: direction^T @ computed_adjoint (with sorted summation)
      ! For INOUT parameters: use cb directly (it contains the computed input adjoint after reverse pass)
      ! For pure inputs: use adjoint directly
      vjp_ad = 0.0d0
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
      ! Compute and sort products for x
      n_products = n
      do i = 1, n
        temp_products(i) = x_dir(i) * xb(k,i)
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

end program test_dtrmv_vector_reverse