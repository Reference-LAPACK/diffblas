! Test program for DTRSV reverse (adjoint) mode differentiation
! Hand-written driver following the structure of test_dgemv_reverse.f90.
! Validates the adjoint via the dot-product (VJP) identity:
!   <xb_seed, dx> = <a_dir,ab> + <x_dir,xb_out>
! Uses REAL*8. Multi-size, sweeps DIAG in {'N','U'}.
! (UPLO='U', TRANS='N' held fixed for now.)

program test_dtrsv_reverse
  implicit none

  external :: dtrsv
  external :: dtrsv_b

  integer :: n_test, seed_array(33), test_sizes(3), i, id
  logical :: passed, all_passed
  character :: diag

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing DTRSV (reverse mode)'
  all_passed = .true.
  do id = 1, 2
    if (id == 1) then
      diag = 'N'
    else
      diag = 'U'
    end if
    do i = 1, 3
      n_test = test_sizes(i)
      call run_test_for_size(n_test, diag, passed)
      all_passed = all_passed .and. passed
    end do
  end do
  if (all_passed) then
    write(*,*) 'PASS: All sizes/diags completed successfully'
  else
    write(*,*) 'FAIL: One or more cases had derivative errors'
  end if

contains

  subroutine run_test_for_size(n, diag, passed)
    implicit none
    integer, intent(in) :: n
    character, intent(in) :: diag
    logical, intent(out) :: passed

    character :: uplo, trans
    integer :: nsize, lda_val, incx
    real(8), dimension(n,n) :: a, ab
    real(8), dimension(n) :: x, xb
    real(8), dimension(n,n) :: a_orig
    real(8), dimension(n) :: x_orig, xb_orig
    integer :: i

    uplo = 'U'
    trans = 'N'
    nsize = n
    lda_val = n
    incx = 1

    call random_number(a)
    a = (a * 2.0d0 - 1.0d0) / real(n, 8)
    do i = 1, n
      a(i,i) = 2.0d0 + abs(a(i,i))
    end do
    call random_number(x)
    x = x * 2.0d0 - 1.0d0

    a_orig = a
    x_orig = x

    call random_number(xb)
    xb = xb * 2.0d0 - 1.0d0
    xb_orig = xb

    ab = 0.0d0

    write(*,*) 'Testing DTRSV (n =', n, ', diag = ', diag, ')'

    call dtrsv_b(uplo, trans, diag, nsize, a, ab, lda_val, x, xb, incx)

    call check_vjp_numerically(n, uplo, trans, diag, nsize, lda_val, incx, &
         a_orig, x_orig, xb_orig, ab, xb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, uplo, trans, diag, nsize, lda_val, &
       incx, a_orig, x_orig, xb_orig, ab, xb, passed)
    implicit none
    integer, intent(in) :: n, nsize, lda_val, incx
    character, intent(in) :: uplo, trans, diag
    real(8), intent(in) :: a_orig(n,n), x_orig(n), xb_orig(n)
    real(8), intent(in) :: ab(n,n), xb(n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0d-7
    real(8) :: vjp_ad, vjp_fd, relative_error, max_error
    real(8) :: abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j

    real(8), dimension(n,n) :: a_dir, a
    real(8), dimension(n) :: x_dir, x, x_plus, x_minus, x_cdiff

    has_large_errors = .false.

    call random_number(a_dir)
    a_dir = a_dir * 2.0d0 - 1.0d0
    call random_number(x_dir)
    x_dir = x_dir * 2.0d0 - 1.0d0

    a = a_orig + h * a_dir
    x = x_orig + h * x_dir
    call dtrsv(uplo, trans, diag, nsize, a, lda_val, x, incx)
    x_plus = x

    a = a_orig - h * a_dir
    x = x_orig - h * x_dir
    call dtrsv(uplo, trans, diag, nsize, a, lda_val, x, incx)
    x_minus = x

    x_cdiff = (x_plus - x_minus) / (2.0d0 * h)

    vjp_fd = 0.0d0
    do i = 1, n
      vjp_fd = vjp_fd + xb_orig(i) * x_cdiff(i)
    end do

    vjp_ad = 0.0d0
    do j = 1, n
      do i = 1, n
        vjp_ad = vjp_ad + a_dir(i,j) * ab(i,j)
      end do
    end do
    do i = 1, n
      vjp_ad = vjp_ad + x_dir(i) * xb(i)
    end do

    abs_error = abs(vjp_fd - vjp_ad)
    abs_reference = abs(vjp_ad)
    error_bound = 1.0d-5 + 1.0d-5 * abs_reference
    if (abs_error > error_bound) has_large_errors = .true.
    if (abs_reference > 1.0d-10) then
      relative_error = abs_error / abs_reference
    else
      relative_error = abs_error
    end if
    max_error = relative_error

    write(*,*) 'VJP AD =', vjp_ad, '  VJP FD =', vjp_fd
    write(*,*) 'Maximum relative error:', max_error
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Adjoint is outside tolerance'
    else
      write(*,*) 'PASS: Adjoint is within tolerance (rtol + atol)'
    end if

  end subroutine check_vjp_numerically

end program test_dtrsv_reverse
