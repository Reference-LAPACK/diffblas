! Test program for CTRSM reverse (adjoint) mode differentiation
! Hand-written driver following the structure of test_zgemv_reverse.f90.
! Validates the adjoint via the real-part Hermitian dot-product identity:
!   Re<bb_seed, dX> = Re<a_dir,ab> + Re<b_dir,bb_out> + Re(conjg(alpha_dir)*alphab)
! where <u,v> = sum conjg(u)*v. COMPLEX*8. Multi-size, sweeps DIAG in {'N','U'}.
! (SIDE='L', UPLO='U', TRANSA='N' held fixed for now.)

program test_ctrsm_reverse
  implicit none

  external :: ctrsm
  external :: ctrsm_b

  integer :: n_test, seed_array(33), test_sizes(3), i, id
  logical :: passed, all_passed
  character :: diag

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing CTRSM (reverse mode)'
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
    complex(4), intent(out) :: z(k)
    integer :: t
    real(8) :: r, s
    do t = 1, k
      call random_number(r)
      call random_number(s)
      z(t) = cmplx(2.0d0*r - 1.0d0, 2.0d0*s - 1.0d0, kind=4)
    end do
  end subroutine fill_c

  subroutine fill_c1(z)
    complex(4), intent(out) :: z
    real(8) :: r, s
    call random_number(r)
    call random_number(s)
    z = cmplx(2.0d0*r - 1.0d0, 2.0d0*s - 1.0d0, kind=4)
  end subroutine fill_c1

  subroutine run_test_for_size(n, diag, passed)
    integer, intent(in) :: n
    character, intent(in) :: diag
    logical, intent(out) :: passed

    character :: side, uplo, transa
    integer :: msize, nsize, lda_val, ldb_val
    complex(4) :: alpha, alphab
    complex(4), dimension(n,n) :: a, b, ab, bb
    complex(4) :: alpha_orig
    complex(4), dimension(n,n) :: a_orig, b_orig, bb_orig
    integer :: i

    side = 'L'; uplo = 'U'; transa = 'N'
    msize = n; nsize = n; lda_val = n; ldb_val = n

    call fill_c1(alpha)
    call fill_c(a, n*n)
    a = a / real(n, 4)
    do i = 1, n
      a(i,i) = cmplx(2.0d0 + abs(real(a(i,i))), aimag(a(i,i)), kind=4)
    end do
    call fill_c(b, n*n)

    alpha_orig = alpha; a_orig = a; b_orig = b

    call fill_c(bb, n*n)
    bb_orig = bb
    ab = (0.0d0, 0.0d0)
    alphab = (0.0d0, 0.0d0)

    write(*,*) 'Testing CTRSM (n =', n, ', diag = ', diag, ')'

    call ctrsm_b(side, uplo, transa, diag, msize, nsize, alpha, alphab, &
                 a, ab, lda_val, b, bb, ldb_val)

    call check_vjp_numerically(n, side, uplo, transa, diag, msize, nsize, &
         lda_val, ldb_val, alpha_orig, a_orig, b_orig, bb_orig, alphab, ab, &
         bb, passed)
  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, side, uplo, transa, diag, msize, nsize, &
       lda_val, ldb_val, alpha_orig, a_orig, b_orig, bb_orig, alphab, ab, &
       bb, passed)
    integer, intent(in) :: n, msize, nsize, lda_val, ldb_val
    character, intent(in) :: side, uplo, transa, diag
    complex(4), intent(in) :: alpha_orig
    complex(4), intent(in) :: a_orig(n,n), b_orig(n,n), bb_orig(n,n)
    complex(4), intent(in) :: alphab, ab(n,n), bb(n,n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0d-4
    real(8) :: vjp_ad, vjp_fd, relative_error, abs_error, abs_reference
    logical :: has_large_errors
    integer :: i, j
    complex(4) :: alpha_dir, alpha
    complex(4), dimension(n,n) :: a_dir, b_dir, a, b, x_plus, x_minus, x_cdiff

    has_large_errors = .false.

    call fill_c1(alpha_dir)
    call fill_c(a_dir, n*n)
    call fill_c(b_dir, n*n)

    alpha = alpha_orig + h * alpha_dir
    a = a_orig + h * a_dir
    b = b_orig + h * b_dir
    call ctrsm(side, uplo, transa, diag, msize, nsize, alpha, a, lda_val, b, ldb_val)
    x_plus = b

    alpha = alpha_orig - h * alpha_dir
    a = a_orig - h * a_dir
    b = b_orig - h * b_dir
    call ctrsm(side, uplo, transa, diag, msize, nsize, alpha, a, lda_val, b, ldb_val)
    x_minus = b

    x_cdiff = (x_plus - x_minus) / (2.0d0 * h)

    vjp_fd = 0.0d0
    do j = 1, n
      do i = 1, n
        vjp_fd = vjp_fd + real(conjg(bb_orig(i,j)) * x_cdiff(i,j))
      end do
    end do

    vjp_ad = real(conjg(alpha_dir) * alphab)
    do j = 1, n
      do i = 1, n
        vjp_ad = vjp_ad + real(conjg(a_dir(i,j)) * ab(i,j)) &
                        + real(conjg(b_dir(i,j)) * bb(i,j))
      end do
    end do

    abs_error = abs(vjp_fd - vjp_ad)
    abs_reference = abs(vjp_ad)
    if (abs_reference > 1.0d-6) then
      relative_error = abs_error / abs_reference
    else
      relative_error = abs_error
    end if
    if (abs_error > 1.0d-2 + 1.0d-2 * abs_reference) has_large_errors = .true.

    write(*,*) 'VJP AD =', vjp_ad, '  VJP FD =', vjp_fd
    write(*,*) 'Maximum relative error:', relative_error
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Adjoint is outside tolerance'
    else
      write(*,*) 'PASS: Adjoint is within tolerance (rtol + atol)'
    end if
  end subroutine check_vjp_numerically

end program test_ctrsm_reverse
