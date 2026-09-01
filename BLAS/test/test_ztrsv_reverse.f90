! Test program for ZTRSV reverse (adjoint) mode differentiation
! Hand-written driver following the structure of test_zgemv_reverse.f90.
! Validates the adjoint via the real-part Hermitian dot-product identity:
!   Re<xb_seed, dx> = Re<a_dir,ab> + Re<x_dir,xb_out>
! where <u,v> = sum conjg(u)*v. COMPLEX*16. Multi-size, sweeps DIAG in {'N','U'}.
! (UPLO='U', TRANS='N' held fixed for now.)

program test_ztrsv_reverse
  implicit none

  external :: ztrsv
  external :: ztrsv_b

  integer :: n_test, seed_array(33), test_sizes(3), i, id
  logical :: passed, all_passed
  character :: diag

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing ZTRSV (reverse mode)'
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

  subroutine fill_c(z, k)
    integer, intent(in) :: k
    complex(8), intent(out) :: z(k)
    integer :: t
    real(8) :: r, s
    do t = 1, k
      call random_number(r)
      call random_number(s)
      z(t) = cmplx(2.0d0*r - 1.0d0, 2.0d0*s - 1.0d0, kind=8)
    end do
  end subroutine fill_c

  subroutine run_test_for_size(n, diag, passed)
    integer, intent(in) :: n
    character, intent(in) :: diag
    logical, intent(out) :: passed

    character :: uplo, trans
    integer :: nsize, lda_val, incx
    complex(8), dimension(n,n) :: a, ab
    complex(8), dimension(n) :: x, xb
    complex(8), dimension(n,n) :: a_orig
    complex(8), dimension(n) :: x_orig, xb_orig
    integer :: i

    uplo = 'U'; trans = 'N'
    nsize = n; lda_val = n; incx = 1

    call fill_c(a, n*n)
    a = a / real(n, 8)
    do i = 1, n
      a(i,i) = cmplx(2.0d0 + abs(real(a(i,i))), aimag(a(i,i)), kind=8)
    end do
    call fill_c(x, n)

    a_orig = a; x_orig = x

    call fill_c(xb, n)
    xb_orig = xb
    ab = (0.0d0, 0.0d0)

    write(*,*) 'Testing ZTRSV (n =', n, ', diag = ', diag, ')'

    call ztrsv_b(uplo, trans, diag, nsize, a, ab, lda_val, x, xb, incx)

    call check_vjp_numerically(n, uplo, trans, diag, nsize, lda_val, incx, &
         a_orig, x_orig, xb_orig, ab, xb, passed)
  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, uplo, trans, diag, nsize, lda_val, &
       incx, a_orig, x_orig, xb_orig, ab, xb, passed)
    integer, intent(in) :: n, nsize, lda_val, incx
    character, intent(in) :: uplo, trans, diag
    complex(8), intent(in) :: a_orig(n,n), x_orig(n), xb_orig(n)
    complex(8), intent(in) :: ab(n,n), xb(n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0d-7
    real(8) :: vjp_ad, vjp_fd, relative_error, abs_error, abs_reference
    logical :: has_large_errors
    integer :: i, j
    complex(8), dimension(n,n) :: a_dir, a
    complex(8), dimension(n) :: x_dir, x, x_plus, x_minus, x_cdiff

    has_large_errors = .false.

    call fill_c(a_dir, n*n)
    call fill_c(x_dir, n)

    a = a_orig + h * a_dir
    x = x_orig + h * x_dir
    call ztrsv(uplo, trans, diag, nsize, a, lda_val, x, incx)
    x_plus = x

    a = a_orig - h * a_dir
    x = x_orig - h * x_dir
    call ztrsv(uplo, trans, diag, nsize, a, lda_val, x, incx)
    x_minus = x

    x_cdiff = (x_plus - x_minus) / (2.0d0 * h)

    vjp_fd = 0.0d0
    do i = 1, n
      vjp_fd = vjp_fd + real(conjg(xb_orig(i)) * x_cdiff(i))
    end do

    vjp_ad = 0.0d0
    do j = 1, n
      do i = 1, n
        vjp_ad = vjp_ad + real(conjg(a_dir(i,j)) * ab(i,j))
      end do
    end do
    do i = 1, n
      vjp_ad = vjp_ad + real(conjg(x_dir(i)) * xb(i))
    end do

    abs_error = abs(vjp_fd - vjp_ad)
    abs_reference = abs(vjp_ad)
    if (abs_reference > 1.0d-10) then
      relative_error = abs_error / abs_reference
    else
      relative_error = abs_error
    end if
    if (abs_error > 1.0d-5 + 1.0d-5 * abs_reference) has_large_errors = .true.

    write(*,*) 'VJP AD =', vjp_ad, '  VJP FD =', vjp_fd
    write(*,*) 'Maximum relative error:', relative_error
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Adjoint is outside tolerance'
    else
      write(*,*) 'PASS: Adjoint is within tolerance (rtol + atol)'
    end if
  end subroutine check_vjp_numerically

end program test_ztrsv_reverse
