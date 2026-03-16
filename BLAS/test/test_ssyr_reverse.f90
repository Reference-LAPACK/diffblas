! Test program for SSYR reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_ssyr_reverse
  implicit none

  external :: ssyr
  external :: ssyr_b

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(3)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing SSYR (multi-size: n = 4)'
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

    character :: uplo
    integer :: nsize
    real(4) :: alpha
    real(4), dimension(n) :: x
    integer :: incx_val
    real(4), dimension(n,n) :: a
    integer :: lda_val
    real(4) :: alphab
    real(4), dimension(n) :: xb
    real(4), dimension(n,n) :: ab
    real(4) :: alpha_orig
    real(4), dimension(n) :: x_orig
    real(4), dimension(n,n) :: a_orig
    real(4), dimension(n,n) :: ab_orig
    integer :: i, j

    nsize = n
    incx_val = 1
    lda_val = n
    uplo = 'U'

    call random_number(alpha)
    alpha = alpha * 2.0 - 1.0
    call random_number(x)
    x = x * 2.0 - 1.0
    call random_number(a)
    a = a * 2.0 - 1.0

    alpha_orig = alpha
    x_orig = x
    a_orig = a

    call random_number(ab)
    ab = ab * 2.0 - 1.0
    ab_orig = ab

    alphab = 0.0
    xb = 0.0

    write(*,*) 'Testing SSYR (n =', n, ')'

    call set_ISIZE1OFX(n)

    call ssyr_b(uplo, nsize, alpha, alphab, x, xb, incx_val, a, ab, lda_val)

    call set_ISIZE1OFX(-1)

    call check_vjp_numerically(n, uplo, nsize, incx_val, lda_val, alpha_orig, x_orig, a_orig, ab_orig, alphab, xb, ab, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, uplo, nsize, incx_val, lda_val, alpha_orig, x_orig, a_orig, ab_orig, alphab, xb, ab, passed)
    implicit none
    integer, intent(in) :: n
    character, intent(in) :: uplo
    integer, intent(in) :: nsize
    integer, intent(in) :: incx_val
    integer, intent(in) :: lda_val
    real(4), intent(in) :: alpha_orig
    real(4), intent(in) :: x_orig(n)
    real(4), intent(in) :: a_orig(n,n)
    real(4), intent(in) :: ab_orig(n,n)
    real(4), intent(in) :: alphab
    real(4), intent(in) :: xb(n)
    real(4), intent(in) :: ab(n,n)
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j, n_products
    real(4), dimension(n) :: temp_products

    real(4) :: alpha_dir
    real(4), dimension(n) :: x_dir
    real(4), dimension(n,n) :: a_dir

    real(4), dimension(n,n) :: a_plus, a_minus, a_central_diff

    real(4) :: alpha
    real(4), dimension(n) :: x
    real(4), dimension(n,n) :: a

    max_error = 0.0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    call random_number(alpha_dir)
    alpha_dir = alpha_dir * 2.0 - 1.0
    call random_number(x_dir)
    x_dir = x_dir * 2.0 - 1.0
    call random_number(a_dir)
    a_dir = a_dir * 2.0 - 1.0

    alpha = alpha_orig + h * alpha_dir
    x = x_orig + h * x_dir
    a = a_orig + h * a_dir
    call ssyr(uplo, nsize, alpha, x, incx_val, a, lda_val)
    a_plus = a

    alpha = alpha_orig - h * alpha_dir
    x = x_orig - h * x_dir
    a = a_orig - h * a_dir
    call ssyr(uplo, nsize, alpha, x, incx_val, a, lda_val)
    a_minus = a

    a_central_diff = (a_plus - a_minus) / (2.0 * h)

    vjp_fd = 0.0
    do j = 1, n
      do i = 1, n
        vjp_fd = vjp_fd + ab_orig(i,j) * a_central_diff(i,j)
      end do
    end do

    vjp_ad = 0.0
    vjp_ad = vjp_ad + alpha_dir * alphab
    n_products = n
    do i = 1, n
      temp_products(i) = x_dir(i) * xb(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    do j = 1, n
      do i = 1, n
        vjp_ad = vjp_ad + a_dir(i,j) * ab(i,j)
      end do
    end do

    abs_error = abs(vjp_fd - vjp_ad)
    abs_reference = abs(vjp_ad)
    error_bound = 2.0e-3 + 2.0e-3 * abs_reference
    if (abs_error > error_bound) has_large_errors = .true.
    if (abs_reference > 1.0e-10) then
      relative_error = abs_error / abs_reference
    else
      relative_error = abs_error
    end if
    max_error = relative_error
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=2.0e-3, atol=2.0e-3'
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

end program test_ssyr_reverse