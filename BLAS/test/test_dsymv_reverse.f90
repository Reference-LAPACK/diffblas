! Test program for DSYMV reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_dsymv_reverse
  implicit none

  external :: dsymv
  external :: dsymv_b

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing DSYMV (multi-size: n = 4)'
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
    integer :: nsize
    real(8) :: alpha
    real(8), dimension(n,n) :: a
    integer :: lda_val
    real(8), dimension(n) :: x
    integer :: incx_val
    real(8) :: beta
    real(8), dimension(n) :: y
    integer :: incy_val
    real(8) :: alphab
    real(8), dimension(n,n) :: ab
    real(8), dimension(n) :: xb
    real(8) :: betab
    real(8), dimension(n) :: yb
    real(8) :: alpha_orig
    real(8), dimension(n,n) :: a_orig
    real(8), dimension(n) :: x_orig
    real(8) :: beta_orig
    real(8), dimension(n) :: y_orig
    real(8), dimension(n) :: yb_orig
    integer :: i, j

    nsize = n
    lda_val = n
    incx_val = 1
    incy_val = 1
    uplo = 'U'

    call random_number(alpha)
    alpha = alpha * 2.0 - 1.0
    call random_number(a)
    a = a * 2.0 - 1.0
    ! Keep perturbations consistent with symmetric a
    do j = 1, n
    do i = j+1, n
    a(i,j) = a(j,i)
    end do
    end do
    call random_number(x)
    x = x * 2.0 - 1.0
    call random_number(beta)
    beta = beta * 2.0 - 1.0
    call random_number(y)
    y = y * 2.0 - 1.0

    alpha_orig = alpha
    a_orig = a
    x_orig = x
    beta_orig = beta
    y_orig = y

    call random_number(yb)
    yb = yb * 2.0 - 1.0
    yb_orig = yb

    alphab = 0.0
    ab = 0.0
    xb = 0.0
    betab = 0.0

    write(*,*) 'Testing DSYMV (n =', n, ')'

    call set_ISIZE1OFX(n)
    call set_ISIZE2OFA(n)

    call dsymv_b(uplo, nsize, alpha, alphab, a, ab, lda_val, x, xb, incx_val, beta, betab, y, yb, incy_val)

    call set_ISIZE1OFX(-1)
    call set_ISIZE2OFA(-1)

    call check_vjp_numerically(n, uplo, nsize, lda_val, incx_val, incy_val, alpha_orig, a_orig, x_orig, beta_orig, y_orig, yb_orig, alphab, ab, xb, betab, yb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, uplo, nsize, lda_val, incx_val, incy_val, alpha_orig, a_orig, x_orig, beta_orig, y_orig, yb_orig, alphab, ab, xb, betab, yb, passed)
    implicit none
    integer, intent(in) :: n
    character, intent(in) :: uplo
    integer, intent(in) :: nsize
    integer, intent(in) :: lda_val
    integer, intent(in) :: incx_val
    integer, intent(in) :: incy_val
    real(8), intent(in) :: alpha_orig
    real(8), intent(in) :: a_orig(n,n)
    real(8), intent(in) :: x_orig(n)
    real(8), intent(in) :: beta_orig
    real(8), intent(in) :: y_orig(n)
    real(8), intent(in) :: yb_orig(n)
    real(8), intent(in) :: alphab
    real(8), intent(in) :: ab(n,n)
    real(8), intent(in) :: xb(n)
    real(8), intent(in) :: betab
    real(8), intent(in) :: yb(n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j, n_products
    real(8), dimension(n) :: temp_products

    real(8) :: alpha_dir
    real(8), dimension(n,n) :: a_dir
    real(8), dimension(n) :: x_dir
    real(8) :: beta_dir
    real(8), dimension(n) :: y_dir

    real(8), dimension(n) :: y_plus, y_minus, y_central_diff

    real(8) :: alpha
    real(8), dimension(n,n) :: a
    real(8), dimension(n) :: x
    real(8) :: beta
    real(8), dimension(n) :: y

    max_error = 0.0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    call random_number(alpha_dir)
    alpha_dir = alpha_dir * 2.0 - 1.0
    call random_number(a_dir)
    a_dir = a_dir * 2.0 - 1.0
    ! Keep perturbations consistent with symmetric a_dir
    do j = 1, n
    do i = j+1, n
    a_dir(i,j) = a_dir(j,i)
    end do
    end do
    call random_number(x_dir)
    x_dir = x_dir * 2.0 - 1.0
    call random_number(beta_dir)
    beta_dir = beta_dir * 2.0 - 1.0
    call random_number(y_dir)
    y_dir = y_dir * 2.0 - 1.0

    alpha = alpha_orig + h * alpha_dir
    a = a_orig + h * a_dir
    x = x_orig + h * x_dir
    beta = beta_orig + h * beta_dir
    y = y_orig + h * y_dir
    call dsymv(uplo, nsize, alpha, a, lda_val, x, incx_val, beta, y, incy_val)
    y_plus = y

    alpha = alpha_orig - h * alpha_dir
    a = a_orig - h * a_dir
    x = x_orig - h * x_dir
    beta = beta_orig - h * beta_dir
    y = y_orig - h * y_dir
    call dsymv(uplo, nsize, alpha, a, lda_val, x, incx_val, beta, y, incy_val)
    y_minus = y

    y_central_diff = (y_plus - y_minus) / (2.0 * h)

    vjp_fd = 0.0
    n_products = n
    do i = 1, n
      temp_products(i) = yb_orig(i) * y_central_diff(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_fd = vjp_fd + temp_products(i)
    end do

    vjp_ad = 0.0
    vjp_ad = vjp_ad + alpha_dir * alphab
    ! Symmetric A: VJP = sum over upper triangle a_dir*(ab(i,j)+ab(j,i))
    do j = 1, n
      do i = 1, j
        if (i .eq. j) then
          vjp_ad = vjp_ad + a_dir(i,j) * ab(i,j)
        else
          vjp_ad = vjp_ad + a_dir(i,j) * (ab(i,j) + ab(j,i))
        end if
      end do
    end do
    n_products = n
    do i = 1, n
      temp_products(i) = x_dir(i) * xb(i)
    end do
    call sort_array(temp_products, n_products)
    do i = 1, n_products
      vjp_ad = vjp_ad + temp_products(i)
    end do
    vjp_ad = vjp_ad + beta_dir * betab
    n_products = n
    do i = 1, n
      temp_products(i) = y_dir(i) * yb(i)
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

end program test_dsymv_reverse