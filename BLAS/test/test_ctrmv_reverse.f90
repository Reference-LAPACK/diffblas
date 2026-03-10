! Test program for CTRMV reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_ctrmv_reverse
  implicit none

  external :: ctrmv
  external :: ctrmv_b

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing CTRMV (multi-size: n = 4)'
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
    character :: diag
    integer :: nsize
    complex(4), dimension(n,n) :: a
    integer :: lda_val
    complex(4), dimension(n) :: x
    integer :: incx_val
    complex(4), dimension(n,n) :: ab
    complex(4), dimension(n) :: xb
    complex(4), dimension(n,n) :: a_orig
    complex(4), dimension(n) :: x_orig
    complex(4), dimension(n) :: xb_orig
    real(4) :: temp_re, temp_im
    integer :: i, j

    nsize = n
    lda_val = n
    incx_val = 1
    uplo = 'U'
    trans = 'N'
    diag = 'N'

    do j = 1, n
      do i = 1, n
        call random_number(temp_re)
        call random_number(temp_im)
        a(i,j) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
      end do
    end do
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      x(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do

    a_orig = a
    x_orig = x

    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      xb(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do
    xb_orig = xb

    ab = 0.0

    write(*,*) 'Testing CTRMV (n =', n, ')'

    call set_ISIZE2OFA(n)

    call ctrmv_b(uplo, trans, diag, nsize, a, ab, lda_val, x, xb, incx_val)

    call set_ISIZE2OFA(-1)

    call check_vjp_numerically(n, uplo, trans, diag, nsize, lda_val, incx_val, a_orig, x_orig, xb_orig, ab, xb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, uplo, trans, diag, nsize, lda_val, incx_val, a_orig, x_orig, xb_orig, ab, xb, passed)
    implicit none
    integer, intent(in) :: n
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: nsize
    integer, intent(in) :: lda_val
    integer, intent(in) :: incx_val
    complex(4), intent(in) :: a_orig(n,n)
    complex(4), intent(in) :: x_orig(n)
    complex(4), intent(in) :: xb_orig(n)
    complex(4), intent(in) :: ab(n,n)
    complex(4), intent(in) :: xb(n)
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j, n_products
    real(4), dimension(n) :: temp_products
    real(4) :: temp_re, temp_im

    complex(4), dimension(n,n) :: a_dir
    complex(4), dimension(n) :: x_dir

    complex(4), dimension(n) :: x_plus, x_minus, x_central_diff

    complex(4), dimension(n,n) :: a
    complex(4), dimension(n) :: x

    max_error = 0.0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    do j = 1, n
      do i = 1, n
        call random_number(temp_re)
        call random_number(temp_im)
        a_dir(i,j) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
      end do
    end do
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      x_dir(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do

    a = a_orig + cmplx(h, 0.0) * a_dir
    x = x_orig + cmplx(h, 0.0) * x_dir
    call ctrmv(uplo, trans, diag, nsize, a, lda_val, x, incx_val)
    x_plus = x

    a = a_orig - cmplx(h, 0.0) * a_dir
    x = x_orig - cmplx(h, 0.0) * x_dir
    call ctrmv(uplo, trans, diag, nsize, a, lda_val, x, incx_val)
    x_minus = x

    x_central_diff = (x_plus - x_minus) / (2.0 * h)

    vjp_fd = 0.0
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(xb_orig(i)) * x_central_diff(i))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do

    vjp_ad = 0.0
    do j = 1, n
      do i = 1, n
        vjp_ad = vjp_ad + real(conjg(a_dir(i,j)) * ab(i,j))
      end do
    end do
    n_products = n
    do i = 1, n
      temp_products(i) = real(conjg(x_dir(i)) * xb(i))
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
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

end program test_ctrmv_reverse