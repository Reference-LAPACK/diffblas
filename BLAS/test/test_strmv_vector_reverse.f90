! Test program for STRMV vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=n
! Multi-size test with outlined run_test_for_size(n) - TRMV/TRSV

program test_strmv_vector_reverse
  implicit none

  external :: strmv
  external :: strmv_bv

  integer :: nbdirs
  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing STRMV (Vector Reverse, multi-size: n =', test_sizes(1), ')'
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    nbdirs = test_sizes(i)
    call run_test_for_size(n_test, passed, nbdirs)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) then
    write(*,*) 'PASS: Vector reverse mode - all sizes completed successfully'
  else
    write(*,*) 'FAIL: Vector reverse mode - one or more sizes had derivative errors'
  end if

contains

  subroutine run_test_for_size(n, passed, nbdirs)
    implicit none
    integer, intent(in) :: n
    logical, intent(out) :: passed
    integer, intent(in) :: nbdirs

    character :: uplo, trans, diag
    integer :: nsize, lda_val, incx_val
    real(4), dimension(n,n) :: a
    real(4), dimension(n) :: x
    real(4), dimension(nbdirs,n,n) :: ab
    real(4), dimension(nbdirs,n) :: xb
    real(4), dimension(n,n) :: a_orig
    real(4), dimension(n) :: x_orig
    real(4), dimension(nbdirs,n) :: xb_orig
    integer :: k, ii, jj
    real(4) :: temp_real, temp_imag

    uplo = 'L'
    trans = 'N'
    diag = 'N'
    nsize = n
    lda_val = n
    incx_val = 1

    ! Lower triangular A (non-unit)
    do jj = 1, n
      do ii = jj, n
        call random_number(temp_real)
        a(ii,jj) = temp_real * 2.0d0 - 1.0d0
      end do
    end do
    do jj = 1, n
      do ii = 1, jj - 1
        a(ii,jj) = 0.0d0
      end do
    end do
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    do k = 1, nbdirs
      call random_number(xb(k,:))
      xb(k,:) = xb(k,:) * 2.0d0 - 1.0d0
    end do

    a_orig = a
    x_orig = x
    xb_orig = xb
    ab = 0.0d0
    xb = xb_orig

    write(*,*) 'Testing STRMV (Vector Reverse, n =', n, ')'

    call set_ISIZE2OFA(n)

    call strmv_bv(uplo, trans, diag, nsize, a, ab, lda_val, x, xb, incx_val, nbdirs)

    call set_ISIZE2OFA(-1)

    call check_vjp_numerically(n, nbdirs, uplo, trans, diag, nsize, lda_val, incx_val, a_orig, x_orig, xb_orig, ab, xb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nbdirs, uplo, trans, diag, nsize, lda_val, incx_val, a_orig, x_orig, xb_orig, ab, xb, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    character, intent(in) :: uplo, trans, diag
    integer, intent(in) :: nsize, lda_val, incx_val
    real(4), intent(in) :: a_orig(n,n)
    real(4), intent(in) :: x_orig(n)
    real(4), intent(in) :: xb_orig(nbdirs,n)
    real(4), intent(in) :: ab(nbdirs,n,n), xb(nbdirs,n)
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    real(4), dimension(n,n) :: a_dir, a
    real(4), dimension(n) :: x_dir, x, x_plus, x_minus, x_central_diff
    real(4), dimension(n) :: temp_real_fd
    integer :: n_products, i, k, ii, jj
    real(4) :: temp_real, temp_imag
    logical :: has_large_errors

    max_error = 0.0d0
    has_large_errors = .false.

    do k = 1, nbdirs
      do jj = 1, n
        do ii = jj, n
          call random_number(temp_real)
          a_dir(ii,jj) = temp_real * 2.0d0 - 1.0d0
        end do
      end do
      do jj = 1, n
        do ii = 1, jj - 1
          a_dir(ii,jj) = 0.0d0
        end do
      end do
      call random_number(x_dir)
      x_dir = x_dir * 2.0d0 - 1.0d0
      a = a_orig + h * a_dir
      x = x_orig + h * x_dir
      call strmv(uplo, trans, diag, nsize, a, lda_val, x, incx_val)
      x_plus = x
      a = a_orig - h * a_dir
      x = x_orig - h * x_dir
      call strmv(uplo, trans, diag, nsize, a, lda_val, x, incx_val)
      x_minus = x
      x_central_diff = (x_plus - x_minus) / (2.0e0 * h)
      vjp_fd = 0.0e0
      n_products = n
      do i = 1, n
        temp_real_fd(i) = xb_orig(k,i) * x_central_diff(i)
      end do
      call sort_array(temp_real_fd, n_products)
      do i = 1, n_products
        vjp_fd = vjp_fd + temp_real_fd(i)
      end do
      vjp_ad = 0.0d0
      ! Triangular A: sum over lower triangle only (same as stored)
      do jj = 1, n
        do ii = jj, n
          vjp_ad = vjp_ad + a_dir(ii,jj) * ab(k,ii,jj)
        end do
      end do
      do ii = 1, n
        vjp_ad = vjp_ad + x_dir(ii) * xb(k,ii)
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
      if (relative_error > max_error) max_error = relative_error
    end do

    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance: rtol=atol=1.0e-3'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors in derivatives'
    else
      write(*,*) 'PASS: Derivatives within tolerance'
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

end program test_strmv_vector_reverse