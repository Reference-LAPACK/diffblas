! Test program for DSYMV vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirs=n
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_dsymv_vector_reverse
  implicit none

  external :: dsymv
  external :: dsymv_bv

  integer :: nbdirs
  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing DSYMV (Vector Reverse, multi-size: n =', test_sizes(1), ')'
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    nbdirs = test_sizes(i)
    call run_test_for_size(n_test, passed, nbdirs)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) then
    write(*,*) 'PASS: All sizes completed successfully'
  else
    write(*,*) 'FAIL: One or more sizes had derivative errors'
  end if

contains

  subroutine run_test_for_size(n, passed, nbdirs)
    implicit none
    integer, intent(in) :: n
    logical, intent(out) :: passed
    integer, intent(in) :: nbdirs

    character :: uplo
    integer :: nsize, lda_val, incx_val, incy_val
    real(8) :: alpha, beta
    real(8), dimension(n,n) :: a
    real(8), dimension(n) :: x, y
    real(8), dimension(nbdirs) :: alphab, betab
    real(8), dimension(nbdirs,n,n) :: ab
    real(8), dimension(nbdirs,n) :: xb, yb
    real(8) :: alpha_orig, beta_orig
    real(8), dimension(n,n) :: a_orig
    real(8), dimension(n) :: x_orig, y_orig
    real(8), dimension(nbdirs,n) :: yb_orig
    integer :: k, ii, jj
    real(4) :: temp_real, temp_imag

    uplo = 'L'
    nsize = n
    lda_val = n
    incx_val = 1
    incy_val = 1

    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(a)
    a = a * 2.0d0 - 1.0d0
    do jj = 1, n
      do ii = jj + 1, n
        a(ii,jj) = a(jj,ii)
      end do
    end do
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    call random_number(beta)
    beta = beta * 2.0d0 - 1.0d0
    call random_number(y)
    y = y * 2.0d0 - 1.0d0

    alpha_orig = alpha
    a_orig = a
    x_orig = x
    beta_orig = beta
    y_orig = y

    do k = 1, nbdirs
      call random_number(yb(k,:))
      yb(k,:) = yb(k,:) * 2.0d0 - 1.0d0
    end do
    yb_orig = yb

    alphab = 0.0d0
    ab = 0.0d0
    xb = 0.0d0
    betab = 0.0d0

    write(*,*) 'Testing DSYMV (Vector Reverse, n =', n, ')'

    call set_ISIZE1OFX(n)
    call set_ISIZE2OFA(n)

    call dsymv_bv(uplo, nsize, alpha, alphab, a, ab, lda_val, x, xb, incx_val, beta, betab, y, yb, incy_val, nbdirs)

    call set_ISIZE1OFX(-1)
    call set_ISIZE2OFA(-1)

    call check_vjp_numerically(n, nbdirs, uplo, nsize, lda_val, incx_val, incy_val, alpha_orig, a_orig, x_orig, beta_orig, y_orig, yb_orig, alphab, ab, xb, betab, yb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nbdirs, uplo, nsize, lda_val, incx_val, incy_val, alpha_orig, a_orig, x_orig, beta_orig, y_orig, yb_orig, alphab, ab, xb, betab, yb, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    character, intent(in) :: uplo
    integer, intent(in) :: nsize, lda_val, incx_val, incy_val
    real(8), intent(in) :: alpha_orig, beta_orig
    real(8), intent(in) :: a_orig(n,n)
    real(8), intent(in) :: x_orig(n), y_orig(n)
    real(8), intent(in) :: yb_orig(nbdirs,n)
    real(8), intent(in) :: alphab(nbdirs), betab(nbdirs)
    real(8), intent(in) :: ab(nbdirs,n,n), xb(nbdirs,n), yb(nbdirs,n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    real(8) :: alpha_dir, beta_dir
    real(8), dimension(n,n) :: a_dir
    real(8), dimension(n) :: x_dir, y_dir
    real(8) :: alpha, beta
    real(8), dimension(n,n) :: a
    real(8), dimension(n) :: x, y, y_plus, y_minus, y_central_diff
    real(8), dimension(n) :: temp_real_fd
    integer :: n_products, i, k, ii, jj
    real(4) :: temp_real, temp_imag
    logical :: has_large_errors

    max_error = 0.0d0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    do k = 1, nbdirs
      call random_number(alpha_dir)
      alpha_dir = alpha_dir * 2.0d0 - 1.0d0
      call random_number(a_dir)
      a_dir = a_dir * 2.0d0 - 1.0d0
      do jj = 1, n
        do ii = jj + 1, n
          a_dir(ii,jj) = a_dir(jj,ii)
        end do
      end do
      call random_number(x_dir)
      x_dir = x_dir * 2.0d0 - 1.0d0
      call random_number(beta_dir)
      beta_dir = beta_dir * 2.0d0 - 1.0d0
      call random_number(y_dir)
      y_dir = y_dir * 2.0d0 - 1.0d0
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
      y_central_diff = (y_plus - y_minus) / (2.0e0 * h)
      vjp_fd = 0.0e0
      n_products = n
      do i = 1, n
        temp_real_fd(i) = yb_orig(k,i) * y_central_diff(i)
      end do
      call sort_array(temp_real_fd, n_products)
      do i = 1, n_products
        vjp_fd = vjp_fd + temp_real_fd(i)
      end do
      vjp_ad = 0.0d0
      vjp_ad = vjp_ad + beta_dir * betab(k)
      vjp_ad = vjp_ad + alpha_dir * alphab(k)
      ! Symmetric A: VJP = sum over upper triangle a_dir*(ab(i,j)+ab(j,i))
      do jj = 1, n
        do ii = 1, jj
          if (ii .eq. jj) then
            vjp_ad = vjp_ad + a_dir(ii,jj) * ab(k,ii,jj)
          else
            vjp_ad = vjp_ad + a_dir(ii,jj) * (ab(k,ii,jj) + ab(k,jj,ii))
          end if
        end do
      end do
      do ii = 1, n
        vjp_ad = vjp_ad + x_dir(ii) * xb(k,ii)
        vjp_ad = vjp_ad + y_dir(ii) * yb(k,ii)
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
      if (relative_error > max_error) max_error = relative_error
    end do

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

end program test_dsymv_vector_reverse