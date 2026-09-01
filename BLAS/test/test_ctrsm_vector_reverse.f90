! Test program for CTRSM vector reverse (adjoint) mode differentiation
! Hand-written driver following the structure of test_zgemv_vector_reverse.f90.
! Per-direction real-part Hermitian VJP identity:
!   Re<bb_seed(k), dX> = Re<a_dir,ab(k)> + Re<b_dir,bb_out(k)>
!                        + Re(conjg(alpha_dir)*alphab(k))
! COMPLEX*8, nbdirs directions (runtime). Sweeps DIAG in {'N','U'}.
! (SIDE='L', UPLO='U', TRANSA='N' held fixed for now.)

program test_ctrsm_vector_reverse
  implicit none

  external :: ctrsm
  external :: ctrsm_bv

  integer :: n_test, seed_array(33), test_sizes(3), i, id
  logical :: passed, all_passed
  character :: diag

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing CTRSM (vector reverse mode)'
  all_passed = .true.
  do id = 1, 2
    if (id == 1) then
      diag = 'N'
    else
      diag = 'U'
    end if
    do i = 1, 3
      n_test = test_sizes(i)
      call run_test_for_size(n_test, n_test, diag, passed)
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

  subroutine run_test_for_size(n, nbdirs, diag, passed)
    integer, intent(in) :: n, nbdirs
    character, intent(in) :: diag
    logical, intent(out) :: passed

    character :: side, uplo, transa
    integer :: msize, nsize, lda_val, ldb_val
    complex(4) :: alpha
    complex(4), dimension(n,n) :: a, b, a0, b0
    complex(4) :: alphab(nbdirs)
    complex(4) :: ab(nbdirs,n,n), bb(nbdirs,n,n), bb_seed(nbdirs,n,n)
    complex(4) :: alpha_dir, a_dir(n,n), b_dir(n,n), altmp
    complex(4), dimension(n,n) :: xp, xm, cdiff, atmp
    real(8), parameter :: h = 1.0d-4
    real(8) :: max_error, vjp_fd, vjp_ad, abs_error, abs_reference, relerr
    logical :: has_large_errors
    integer :: nd, i, j

    side = 'L'; uplo = 'U'; transa = 'N'
    msize = n; nsize = n; lda_val = n; ldb_val = n

    call fill_c1(alpha)
    call fill_c(a, n*n)
    a = a / real(n, 4)
    do i = 1, n
      a(i,i) = cmplx(2.0d0 + abs(real(a(i,i))), aimag(a(i,i)), kind=4)
    end do
    call fill_c(b, n*n)
    a0 = a; b0 = b

    call fill_c(bb, nbdirs*n*n)
    bb_seed = bb
    ab = (0.0d0, 0.0d0)
    alphab = (0.0d0, 0.0d0)

    write(*,*) 'Testing CTRSM (n =', n, ', nbdirs =', nbdirs, ', diag = ', diag, ')'

    call ctrsm_bv(side, uplo, transa, diag, msize, nsize, alpha, alphab, &
                  a, ab, lda_val, b, bb, ldb_val, nbdirs)

    max_error = 0.0d0
    has_large_errors = .false.
    do nd = 1, nbdirs
      call fill_c1(alpha_dir)
      call fill_c(a_dir, n*n)
      call fill_c(b_dir, n*n)

      altmp = alpha + h*alpha_dir
      atmp = a0 + h*a_dir
      xp = b0 + h * b_dir
      call ctrsm(side, uplo, transa, diag, msize, nsize, altmp, atmp, lda_val, xp, ldb_val)
      altmp = alpha - h*alpha_dir
      atmp = a0 - h*a_dir
      xm = b0 - h * b_dir
      call ctrsm(side, uplo, transa, diag, msize, nsize, altmp, atmp, lda_val, xm, ldb_val)
      cdiff = (xp - xm) / (2.0d0 * h)

      vjp_fd = 0.0d0
      do j = 1, n
        do i = 1, n
          vjp_fd = vjp_fd + real(conjg(bb_seed(nd,i,j)) * cdiff(i,j))
        end do
      end do

      vjp_ad = real(conjg(alpha_dir) * alphab(nd))
      do j = 1, n
        do i = 1, n
          vjp_ad = vjp_ad + real(conjg(a_dir(i,j)) * ab(nd,i,j)) &
                          + real(conjg(b_dir(i,j)) * bb(nd,i,j))
        end do
      end do

      abs_error = abs(vjp_fd - vjp_ad)
      abs_reference = abs(vjp_ad)
      if (abs_error > 1.0d-2 + 1.0d-2 * abs_reference) has_large_errors = .true.
      relerr = abs_error / max(abs_reference, 1.0d-6)
      max_error = max(max_error, relerr)
    end do

    write(*,*) 'Maximum relative error:', max_error
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Adjoint is outside tolerance'
    else
      write(*,*) 'PASS: Adjoint is within tolerance (rtol + atol)'
    end if
  end subroutine run_test_for_size

end program test_ctrsm_vector_reverse
