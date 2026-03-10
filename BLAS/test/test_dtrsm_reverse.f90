! Test program for DTRSM reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_dtrsm_reverse
  implicit none

  external :: dtrsm
  external :: dtrsm_b

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing DTRSM (multi-size: n = 4)'
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

    character :: side
    character :: uplo
    character :: transa
    character :: diag
    integer :: msize
    integer :: nsize
    real(8) :: alpha
    real(8), dimension(n,n) :: a
    integer :: lda_val
    real(8), dimension(n,n) :: b
    integer :: ldb_val
    real(8) :: alphab
    real(8), dimension(n,n) :: ab
    real(8), dimension(n,n) :: bb
    real(8) :: alpha_orig
    real(8), dimension(n,n) :: a_orig
    real(8), dimension(n,n) :: b_orig
    real(8), dimension(n,n) :: bb_orig
    integer :: i, j

    nsize = n
    msize = n
    lda_val = n
    ldb_val = n
    side = 'L'
    uplo = 'U'
    transa = 'N'
    diag = 'N'

    call random_number(alpha)
    alpha = alpha * 2.0 - 1.0
    call random_number(a)
    a = a * 2.0 - 1.0
    call random_number(b)
    b = b * 2.0 - 1.0

    alpha_orig = alpha
    a_orig = a
    b_orig = b

    call random_number(bb)
    bb = bb * 2.0 - 1.0
    bb_orig = bb

    alphab = 0.0
    ab = 0.0

    write(*,*) 'Testing DTRSM (n =', n, ')'

    call set_ISIZE2OFA(n)

    call dtrsm_b(side, uplo, transa, diag, msize, nsize, alpha, alphab, a, ab, lda_val, b, bb, ldb_val)

    call set_ISIZE2OFA(-1)

    call check_vjp_numerically(n, side, uplo, transa, diag, msize, nsize, lda_val, ldb_val, alpha_orig, a_orig, b_orig, bb_orig, alphab, ab, bb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, side, uplo, transa, diag, msize, nsize, lda_val, ldb_val, alpha_orig, a_orig, b_orig, bb_orig, alphab, ab, bb, passed)
    implicit none
    integer, intent(in) :: n
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    integer, intent(in) :: msize
    integer, intent(in) :: nsize
    integer, intent(in) :: lda_val
    integer, intent(in) :: ldb_val
    real(8), intent(in) :: alpha_orig
    real(8), intent(in) :: a_orig(n,n)
    real(8), intent(in) :: b_orig(n,n)
    real(8), intent(in) :: bb_orig(n,n)
    real(8), intent(in) :: alphab
    real(8), intent(in) :: ab(n,n)
    real(8), intent(in) :: bb(n,n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j, n_products
    real(8), dimension(n) :: temp_products

    real(8) :: alpha_dir
    real(8), dimension(n,n) :: a_dir
    real(8), dimension(n,n) :: b_dir

    real(8), dimension(n,n) :: b_plus, b_minus, b_central_diff

    real(8) :: alpha
    real(8), dimension(n,n) :: a
    real(8), dimension(n,n) :: b

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

    alpha = alpha_orig + h * alpha_dir
    a = a_orig + h * a_dir
    b = b_orig + h * b_dir
    call dtrsm(side, uplo, transa, diag, msize, nsize, alpha, a, lda_val, b, ldb_val)
    b_plus = b

    alpha = alpha_orig - h * alpha_dir
    a = a_orig - h * a_dir
    b = b_orig - h * b_dir
    call dtrsm(side, uplo, transa, diag, msize, nsize, alpha, a, lda_val, b, ldb_val)
    b_minus = b

    b_central_diff = (b_plus - b_minus) / (2.0 * h)

    vjp_fd = 0.0
    do j = 1, n
      do i = 1, n
        vjp_fd = vjp_fd + bb_orig(i,j) * b_central_diff(i,j)
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

end program test_dtrsm_reverse