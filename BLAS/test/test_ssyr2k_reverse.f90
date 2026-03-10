! Test program for SSYR2K reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_ssyr2k_reverse
  implicit none

  external :: ssyr2k
  external :: ssyr2k_b

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing SSYR2K (multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 1
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

    character :: uplo
    character :: trans
    integer :: nsize
    integer :: ksize
    real(4) :: alpha
    real(4), dimension(n,n) :: a
    integer :: lda_val
    real(4), dimension(n,n) :: b
    integer :: ldb_val
    real(4) :: beta
    real(4), dimension(n,n) :: c
    integer :: ldc_val
    real(4) :: alphab
    real(4), dimension(n,n) :: ab
    real(4), dimension(n,n) :: bb
    real(4) :: betab
    real(4), dimension(n,n) :: cb
    real(4) :: alpha_orig
    real(4), dimension(n,n) :: a_orig
    real(4), dimension(n,n) :: b_orig
    real(4) :: beta_orig
    real(4), dimension(n,n) :: c_orig
    real(4), dimension(n,n) :: cb_orig
    integer :: i, j

    nsize = n
    ksize = n
    lda_val = n
    ldb_val = n
    ldc_val = n
    uplo = 'U'
    trans = 'N'

    call random_number(alpha)
    alpha = alpha * 2.0 - 1.0
    call random_number(a)
    a = a * 2.0 - 1.0
    call random_number(b)
    b = b * 2.0 - 1.0
    call random_number(beta)
    beta = beta * 2.0 - 1.0
    call random_number(c)
    c = c * 2.0 - 1.0

    alpha_orig = alpha
    a_orig = a
    b_orig = b
    beta_orig = beta
    c_orig = c

    call random_number(cb)
    cb = cb * 2.0 - 1.0
    cb_orig = cb

    alphab = 0.0
    ab = 0.0
    bb = 0.0
    betab = 0.0

    write(*,*) 'Testing SSYR2K (n =', n, ')'

    call set_ISIZE2OFA(n)
    call set_ISIZE2OFB(n)

    call ssyr2k_b(uplo, trans, nsize, ksize, alpha, alphab, a, ab, lda_val, b, bb, ldb_val, beta, betab, c, cb, ldc_val)

    call set_ISIZE2OFA(-1)
    call set_ISIZE2OFB(-1)

    call check_vjp_numerically(n, uplo, trans, nsize, ksize, lda_val, ldb_val, ldc_val, alpha_orig, a_orig, b_orig, beta_orig, c_orig, cb_orig, alphab, ab, bb, betab, cb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, uplo, trans, nsize, ksize, lda_val, ldb_val, ldc_val, alpha_orig, a_orig, b_orig, beta_orig, c_orig, cb_orig, alphab, ab, bb, betab, cb, passed)
    implicit none
    integer, intent(in) :: n
    character, intent(in) :: uplo
    character, intent(in) :: trans
    integer, intent(in) :: nsize
    integer, intent(in) :: ksize
    integer, intent(in) :: lda_val
    integer, intent(in) :: ldb_val
    integer, intent(in) :: ldc_val
    real(4), intent(in) :: alpha_orig
    real(4), intent(in) :: a_orig(n,n)
    real(4), intent(in) :: b_orig(n,n)
    real(4), intent(in) :: beta_orig
    real(4), intent(in) :: c_orig(n,n)
    real(4), intent(in) :: cb_orig(n,n)
    real(4), intent(in) :: alphab
    real(4), intent(in) :: ab(n,n)
    real(4), intent(in) :: bb(n,n)
    real(4), intent(in) :: betab
    real(4), intent(in) :: cb(n,n)
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j, n_products
    real(4), dimension(n) :: temp_products

    real(4) :: alpha_dir
    real(4), dimension(n,n) :: a_dir
    real(4), dimension(n,n) :: b_dir
    real(4) :: beta_dir
    real(4), dimension(n,n) :: c_dir

    real(4), dimension(n,n) :: c_plus, c_minus, c_central_diff

    real(4) :: alpha
    real(4), dimension(n,n) :: a
    real(4), dimension(n,n) :: b
    real(4) :: beta
    real(4), dimension(n,n) :: c

    max_error = 0.0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    call random_number(alpha_dir)
    alpha_dir = alpha_dir * 2.0 - 1.0
    call random_number(a_dir)
    a_dir = a_dir * 2.0 - 1.0
    call random_number(b_dir)
    b_dir = b_dir * 2.0 - 1.0
    call random_number(beta_dir)
    beta_dir = beta_dir * 2.0 - 1.0
    call random_number(c_dir)
    c_dir = c_dir * 2.0 - 1.0

    alpha = alpha_orig + h * alpha_dir
    a = a_orig + h * a_dir
    b = b_orig + h * b_dir
    beta = beta_orig + h * beta_dir
    c = c_orig + h * c_dir
    call ssyr2k(uplo, trans, nsize, ksize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
    c_plus = c

    alpha = alpha_orig - h * alpha_dir
    a = a_orig - h * a_dir
    b = b_orig - h * b_dir
    beta = beta_orig - h * beta_dir
    c = c_orig - h * c_dir
    call ssyr2k(uplo, trans, nsize, ksize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
    c_minus = c

    c_central_diff = (c_plus - c_minus) / (2.0 * h)

    vjp_fd = 0.0
    do j = 1, n
      do i = 1, n
        vjp_fd = vjp_fd + cb_orig(i,j) * c_central_diff(i,j)
      end do
    end do

    vjp_ad = 0.0
    vjp_ad = vjp_ad + alpha_dir * alphab
    do j = 1, n
      do i = 1, n
        vjp_ad = vjp_ad + a_dir(i,j) * ab(i,j)
      end do
    end do
    do j = 1, n
      do i = 1, n
        vjp_ad = vjp_ad + b_dir(i,j) * bb(i,j)
      end do
    end do
    vjp_ad = vjp_ad + beta_dir * betab
    do j = 1, n
      do i = 1, n
        vjp_ad = vjp_ad + c_dir(i,j) * cb(i,j)
      end do
    end do

    abs_error = abs(vjp_fd - vjp_ad)
    abs_reference = abs(vjp_ad)
    error_bound = 1.0e-3 + 1.0e-3 * abs_reference
    if (abs_error > error_bound) has_large_errors = .true.
    if (abs_reference > 1.0e-10) then
      relative_error = abs_error / abs_reference
    else
      relative_error = abs_error
    end if
    max_error = relative_error

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

end program test_ssyr2k_reverse