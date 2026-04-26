! Test program for DGEMM reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_dgemm_reverse
  implicit none

  external :: dgemm
  external :: dgemm_b

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(3)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing DGEMM (multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 3
    n_test = test_sizes(i)
    call run_test_for_size(n_test, passed)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) then
    write(*,*) 'PASS: All sizes completed successfully'
  else
    write(*,*) 'FAIL: One or more sizes had derivative errors'
  end if

contains

  subroutine run_test_for_size(n, passed)
    implicit none
    integer, intent(in) :: n
    logical, intent(out) :: passed

    character :: transa, transb
    integer :: msize, nsize, ksize, lda_val, ldb_val, ldc_val
    real(8) :: alpha, beta
    real(8), dimension(n,n) :: a, b, c
    real(8) :: alphab, betab
    real(8), dimension(n,n) :: ab, bb, cb
    real(8) :: alpha_orig, beta_orig
    real(8), dimension(n,n) :: a_orig, b_orig, c_orig, cb_orig
    integer :: i, j

    transa = 'N'
    transb = 'N'
    msize = n
    nsize = n
    ksize = n
    lda_val = n
    ldb_val = n
    ldc_val = n

    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(a)
    a = a * 2.0d0 - 1.0d0
    call random_number(b)
    b = b * 2.0d0 - 1.0d0
    call random_number(beta)
    beta = beta * 2.0d0 - 1.0d0
    call random_number(c)
    c = c * 2.0d0 - 1.0d0

    alpha_orig = alpha
    a_orig = a
    b_orig = b
    beta_orig = beta
    c_orig = c

    call random_number(cb)
    cb = cb * 2.0d0 - 1.0d0
    cb_orig = cb

    alphab = 0.0d0
    ab = 0.0d0
    bb = 0.0d0
    betab = 0.0d0

    write(*,*) 'Testing DGEMM (n =', n, ')'

    call set_ISIZE2OFA(n)
    call set_ISIZE2OFB(n)

    call dgemm_b(transa, transb, msize, nsize, ksize, alpha, alphab, a, ab, lda_val, b, bb, ldb_val, beta, betab, c, cb, ldc_val)

    call set_ISIZE2OFA(-1)
    call set_ISIZE2OFB(-1)

    call check_vjp_numerically(n, transa, transb, msize, nsize, ksize, lda_val, ldb_val, ldc_val, alpha_orig, a_orig, b_orig, beta_orig, c_orig, cb_orig, alphab, ab, bb, betab, cb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, transa, transb, msize, nsize, ksize, lda_val, ldb_val, ldc_val, alpha_orig, a_orig, b_orig, beta_orig, c_orig, cb_orig, alphab, ab, bb, betab, cb, passed)
    implicit none
    integer, intent(in) :: n
    character, intent(in) :: transa, transb
    integer, intent(in) :: msize, nsize, ksize, lda_val, ldb_val, ldc_val
    real(8), intent(in) :: alpha_orig, beta_orig
    real(8), intent(in) :: a_orig(n,n), b_orig(n,n), c_orig(n,n), cb_orig(n,n)
    real(8), intent(in) :: alphab, betab
    real(8), intent(in) :: ab(n,n), bb(n,n), cb(n,n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    real(8) :: alpha_dir, beta_dir
    real(8), dimension(n,n) :: a_dir, b_dir, c_dir
    real(8), dimension(n,n) :: c_plus, c_minus, c_central_diff
    real(8) :: alpha, beta
    real(8), dimension(n,n) :: a, b, c
    real(8), dimension(n*n) :: temp_products
    integer :: n_products, i, j
    logical :: has_large_errors

    max_error = 0.0d0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    call random_number(alpha_dir)
    alpha_dir = alpha_dir * 2.0d0 - 1.0d0
    call random_number(a_dir)
    a_dir = a_dir * 2.0d0 - 1.0d0
    call random_number(b_dir)
    b_dir = b_dir * 2.0d0 - 1.0d0
    call random_number(beta_dir)
    beta_dir = beta_dir * 2.0d0 - 1.0d0
    call random_number(c_dir)
    c_dir = c_dir * 2.0d0 - 1.0d0

    alpha = alpha_orig + h * alpha_dir
    a = a_orig + h * a_dir
    b = b_orig + h * b_dir
    beta = beta_orig + h * beta_dir
    c = c_orig + h * c_dir
    call dgemm(transa, transb, msize, nsize, ksize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
    c_plus = c

    alpha = alpha_orig - h * alpha_dir
    a = a_orig - h * a_dir
    b = b_orig - h * b_dir
    beta = beta_orig - h * beta_dir
    c = c_orig - h * c_dir
    call dgemm(transa, transb, msize, nsize, ksize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
    c_minus = c

    c_central_diff = (c_plus - c_minus) / (2.0d0 * h)

    vjp_fd = 0.0d0
    n_products = 0
    do j = 1, n
      do i = 1, n
        n_products = n_products + 1
        temp_products(n_products) = cb_orig(i,j) * c_central_diff(i,j)
      end do
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do

    vjp_ad = 0.0d0
    vjp_ad = vjp_ad + alpha_dir * alphab
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
    n_products = 0
    do j = 1, n
      do i = 1, n
        n_products = n_products + 1
        temp_products(n_products) = b_dir(i,j) * bb(i,j)
      end do
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    vjp_ad = vjp_ad + beta_dir * betab
    n_products = 0
    do j = 1, n
      do i = 1, n
        n_products = n_products + 1
        temp_products(n_products) = c_dir(i,j) * cb(i,j)
      end do
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do

    abs_error = abs(vjp_fd - vjp_ad)
    abs_reference = abs(vjp_ad)
    error_bound = 1.0e-5 + 1.0e-5 * abs_reference
    if (abs_error > error_bound) has_large_errors = .true.
    if (abs_reference > 1.0e-10) then
      relative_error = abs_error / abs_reference
    else
      relative_error = abs_error
    end if
    max_error = relative_error
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
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
    do i = 1, n-1
      min_idx = i
      do j = i+1, n
        if (abs(arr(j)) < abs(arr(min_idx))) min_idx = j
      end do
      if (min_idx /= i) then
        temp = arr(i)
        arr(i) = arr(min_idx)
        arr(min_idx) = temp
      end if
    end do
  end subroutine sort_array

end program test_dgemm_reverse