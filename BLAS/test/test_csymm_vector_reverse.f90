! Test program for CSYMM vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=4

program test_csymm_vector_reverse
  implicit none
  integer, parameter :: nbdirs = 4

  external :: csymm
  external :: csymm_bv

  ! Test parameters
  integer :: n  ! Current size (set in loop)
  integer, parameter :: max_size = 100  ! Maximum array dimension (multi-size: 1,4,40,100)
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, k  ! Loop counters
  integer :: test_sizes(1), itest
  logical :: passed, all_passed
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  character :: side
  character :: uplo
  integer :: msize
  integer :: nsize
  complex(4) :: alpha
  complex(4), dimension(max_size,max_size) :: a
  integer :: lda_val
  complex(4), dimension(max_size,max_size) :: b
  integer :: ldb_val
  complex(4) :: beta
  complex(4), dimension(max_size,max_size) :: c
  integer :: ldc_val

  ! Adjoint variables (reverse vector mode)
  ! In reverse mode: output adjoints are INPUT (cotangents/seeds)
  !                  input adjoints are OUTPUT (computed gradients)
  complex(4), dimension(nbdirs) :: alphab
  complex(4), dimension(nbdirs,max_size,max_size) :: ab
  complex(4), dimension(nbdirs,max_size,max_size) :: bb
  complex(4), dimension(nbdirs) :: betab
  complex(4), dimension(nbdirs,max_size,max_size) :: cb

  ! Storage for original cotangents (for INOUT parameters in VJP verification)
  complex(4), dimension(nbdirs,max_size,max_size) :: cb_orig

  ! Storage for original values (for VJP verification)
  complex(4) :: alpha_orig
  complex(4), dimension(max_size,max_size) :: a_orig
  complex(4), dimension(max_size,max_size) :: b_orig
  complex(4) :: beta_orig
  complex(4), dimension(max_size,max_size) :: c_orig

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
  write(*,*) 'Testing CSYMM (Vector Reverse, multi-size: n = 4)'
  all_passed = .true.
  do itest = 1, 1
    n = test_sizes(itest)
    write(*,*) 'Testing CSYMM (Vector Reverse, n =', n, ')'

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
    side = 'L'
    uplo = 'U'
    msize = n
    nsize = n
    call random_number(temp_real)
    call random_number(temp_imag)
    alpha = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
    do j = 1, n
      do i = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        a(i,j) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
      end do
    end do
    lda_val = lda
    do j = 1, n
      do i = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        b(i,j) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
      end do
    end do
    ldb_val = ldb
    call random_number(temp_real)
    call random_number(temp_imag)
    beta = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
    do j = 1, n
      do i = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        c(i,j) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
      end do
    end do
    ldc_val = ldc
    
    ! Store original primal values
    alpha_orig = alpha
    a_orig = a
    b_orig = b
    beta_orig = beta
    c_orig = c
    
    ! Initialize output adjoints (cotangents) with random values for each direction
    ! These are the 'seeds' for reverse mode
    do k = 1, nbdirs
      do j = 1, n
        do i = 1, n
          call random_number(temp_real)
          call random_number(temp_imag)
          cb(k,i,j) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
        end do
      end do
    end do
    
    ! Initialize input adjoints to zero (they will be computed)
    ! Note: Inout parameters are skipped - they already have output adjoints initialized
    alphab = 0.0
    ab = 0.0
    bb = 0.0
    betab = 0.0
    
    ! Save original cotangent seeds for OUTPUT/INOUT parameters (before function call)
    cb_orig = cb
    
    ! Set ISIZE globals required by differentiated routine (dimension 2 of arrays).
    ! ISIZE1OF* (vectors): use n to match adjoint array size; ISIZE2OF* (matrices): use max_size.
    call set_ISIZE2OFA(max_size)
    call set_ISIZE2OFB(max_size)
    
    ! Call reverse vector mode differentiated function
    call csymm_bv(side, uplo, msize, nsize, alpha, alphab, a, ab, lda_val, b, bb, ldb_val, beta, betab, c, cb, ldc_val, nbdirs)
    
    ! Reset ISIZE globals to uninitialized (-1) for completeness
    call set_ISIZE2OFA(-1)
    call set_ISIZE2OFB(-1)
    
    ! VJP Verification using finite differences
    call check_vjp_numerically(passed)
  end subroutine run_test_for_size

  subroutine check_vjp_numerically(passed)
    implicit none
    logical, intent(out) :: passed
    
    ! Direction vectors for VJP testing
    complex(4) :: alpha_dir
    complex(4), dimension(max_size,max_size) :: a_dir
    complex(4), dimension(max_size,max_size) :: b_dir
    complex(4) :: beta_dir
    complex(4), dimension(max_size,max_size) :: c_dir
    complex(4), dimension(max_size,max_size) :: c_plus, c_minus, c_central_diff
    
    max_error = 0.0d0
    has_large_errors = .false.
    
    write(*,*) 'Function calls completed successfully'
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Test each differentiation direction separately
    do k = 1, nbdirs
      
      ! Initialize random direction vectors for all inputs
      call random_number(temp_real)
      call random_number(temp_imag)
      alpha_dir = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
      do j = 1, n
        do i = 1, n
          call random_number(temp_real)
          call random_number(temp_imag)
          a_dir(i,j) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
        end do
      end do
      do j = 1, n
        do i = 1, n
          call random_number(temp_real)
          call random_number(temp_imag)
          b_dir(i,j) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
        end do
      end do
      call random_number(temp_real)
      call random_number(temp_imag)
      beta_dir = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
      do j = 1, n
        do i = 1, n
          call random_number(temp_real)
          call random_number(temp_imag)
          c_dir(i,j) = cmplx(temp_real * 2.0 - 1.0, temp_imag * 2.0 - 1.0)
        end do
      end do
      
      ! Forward perturbation: f(x + h*dir)
      alpha = alpha_orig + cmplx(h, 0.0) * alpha_dir
      a = a_orig + cmplx(h, 0.0) * a_dir
      b = b_orig + cmplx(h, 0.0) * b_dir
      beta = beta_orig + cmplx(h, 0.0) * beta_dir
      c = c_orig + cmplx(h, 0.0) * c_dir
      call csymm(side, uplo, msize, nsize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
      c_plus = c
      
      ! Backward perturbation: f(x - h*dir)
      alpha = alpha_orig - cmplx(h, 0.0) * alpha_dir
      a = a_orig - cmplx(h, 0.0) * a_dir
      b = b_orig - cmplx(h, 0.0) * b_dir
      beta = beta_orig - cmplx(h, 0.0) * beta_dir
      c = c_orig - cmplx(h, 0.0) * c_dir
      call csymm(side, uplo, msize, nsize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
      c_minus = c
      
      ! Compute central differences and VJP verification
      ! VJP check: direction^T @ adjoint should equal finite difference
      
      ! Compute central differences: (f(x+h*dir) - f(x-h*dir)) / (2h)
      c_central_diff = (c_plus - c_minus) / (2.0 * h)
      
      ! VJP verification:
      ! cotangent^T @ central_diff should equal direction^T @ computed_adjoint
      ! Left side: cotangent^T @ Jacobian @ direction (via finite differences, with sorted summation)
      vjp_fd = 0.0
      ! Compute and sort products for c (FD)
      n_products = 0
      do j = 1, n
        do i = 1, n
          n_products = n_products + 1
          temp_products(n_products) = real(conjg(cb_orig(k,i,j)) * c_central_diff(i,j))
        end do
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_fd = vjp_fd + temp_products(i)
      end do
      
      ! Right side: direction^T @ computed_adjoint (with sorted summation)
      ! For INOUT parameters: use cb directly (it contains the computed input adjoint after reverse pass)
      ! For pure inputs: use adjoint directly
      vjp_ad = 0.0
      ! Compute and sort products for a
      n_products = 0
      do j = 1, n
        do i = 1, n
          n_products = n_products + 1
          temp_products(n_products) = real(conjg(a_dir(i,j)) * ab(k,i,j))
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
          temp_products(n_products) = real(conjg(b_dir(i,j)) * bb(k,i,j))
        end do
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_ad = vjp_ad + temp_products(i)
      end do
      vjp_ad = vjp_ad + real(conjg(alpha_dir) * alphab(k))
      ! Compute and sort products for c
      n_products = 0
      do j = 1, n
        do i = 1, n
          n_products = n_products + 1
          temp_products(n_products) = real(conjg(c_dir(i,j)) * cb(k,i,j))
        end do
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_ad = vjp_ad + temp_products(i)
      end do
      vjp_ad = vjp_ad + real(conjg(beta_dir) * betab(k))
      
      ! Error check: |vjp_fd - vjp_ad| > atol + rtol * |vjp_ad|
      abs_error = abs(vjp_fd - vjp_ad)
      abs_reference = abs(vjp_ad)
      error_bound = 1.0e-3 + 1.0e-3 * abs_reference
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
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
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

end program test_csymm_vector_reverse