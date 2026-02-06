! Test program for SNRM2 vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirsmax=4

program test_snrm2_vector_reverse
  use DIFFSIZES
  implicit none

  real(4), external :: snrm2
  external :: snrm2_bv

  ! Test parameters
  integer, parameter :: n = 4  ! Matrix/vector size for test
  integer, parameter :: max_size = n  ! Maximum array dimension
  integer, parameter :: lda = max_size, ldb = max_size, ldc = max_size  ! Leading dimensions
  integer :: i, j, k  ! Loop counters
  integer :: seed_array(33)  ! Random seed
  real(4) :: temp_real, temp_imag  ! Temporary variables for complex initialization

  integer :: nsize
  real(4), dimension(4) :: x
  integer :: incx_val

  ! Adjoint variables (reverse vector mode)
  ! In reverse mode: output adjoints are INPUT (cotangents/seeds)
  !                  input adjoints are OUTPUT (computed gradients)
  real(4), dimension(nbdirsmax,4) :: xb
  real(4), dimension(nbdirsmax) :: snrm2b

  ! Storage for original cotangents (for INOUT parameters in VJP verification)
  real(4), dimension(nbdirsmax) :: snrm2b_orig

  ! Storage for original values (for VJP verification)
  real(4), dimension(4) :: x_orig

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
  nsize = n
  call random_number(x)
  x = x * 2.0 - 1.0
  incx_val = 1

  ! Store original primal values
  x_orig = x

  ! Initialize output adjoints (cotangents) with random values for each direction
  ! These are the 'seeds' for reverse mode
  ! Initialize function result adjoint (output cotangent)
  do k = 1, nbdirsmax
    call random_number(snrm2b(k))
    snrm2b(k) = snrm2b(k) * 2.0 - 1.0
  end do

  ! Initialize input adjoints to zero (they will be computed)
  ! Note: Inout parameters are skipped - they already have output adjoints initialized
  xb = 0.0

  ! Save original cotangent seeds for OUTPUT/INOUT parameters (before function call)
  snrm2b_orig = snrm2b

  ! Call reverse vector mode differentiated function
  call snrm2_bv(nsize, x, xb, incx_val, snrm2b, nbdirsmax)

  ! VJP Verification using finite differences
  call check_vjp_numerically()

  write(*,*) ''
  write(*,*) 'Test completed successfully'

contains

  subroutine check_vjp_numerically()
    implicit none
    
    ! Direction vectors for VJP testing
    real(4), dimension(4) :: x_dir
    real(4) :: snrm2_plus, snrm2_minus
    
    max_error = 0.0d0
    has_large_errors = .false.
    
    write(*,*) 'Function calls completed successfully'
    
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    
    ! Test each differentiation direction separately
    do k = 1, nbdirsmax
      
      ! Initialize random direction vectors for all inputs
      call random_number(x_dir)
      x_dir = x_dir * 2.0 - 1.0
      
      ! Forward perturbation: f(x + h*dir)
      x = x_orig + h * x_dir
      snrm2_plus = snrm2(nsize, x, incx_val)
      
      ! Backward perturbation: f(x - h*dir)
      x = x_orig - h * x_dir
      snrm2_minus = snrm2(nsize, x, incx_val)
      
      ! Compute central differences and VJP verification
      ! VJP check: direction^T @ adjoint should equal finite difference
      
      ! Compute finite difference VJP (central difference)
      ! For functions: vjp_fd = adjoint * central_diff
      vjp_fd = snrm2b(k) * (snrm2_plus - snrm2_minus) / (2.0 * h)
      
      ! Right side: direction^T @ computed_adjoint (with sorted summation)
      ! For INOUT parameters: use cb directly (it contains the computed input adjoint after reverse pass)
      ! For pure inputs: use adjoint directly
      vjp_ad = 0.0
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

end program test_snrm2_vector_reverse