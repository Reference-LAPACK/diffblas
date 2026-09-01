! Test program for CTRSV vector reverse (adjoint) mode differentiation
! Hand-written driver following the structure of test_zgemv_vector_reverse.f90.
! Per-direction real-part Hermitian VJP identity:
!   Re<xb_seed(k), dx> = Re<a_dir,ab(k)> + Re<x_dir,xb_out(k)>
! COMPLEX*8, nbdirs directions (runtime). Sweeps DIAG in {'N','U'}.
! (UPLO='U', TRANS='N' held fixed for now.)

program test_ctrsv_vector_reverse
  implicit none

  external :: ctrsv
  external :: ctrsv_bv

  integer :: n_test, seed_array(33), test_sizes(3), i, id
  logical :: passed, all_passed
  character :: diag

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing CTRSV (vector reverse mode)'
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

  subroutine run_test_for_size(n, nbdirs, diag, passed)
    integer, intent(in) :: n, nbdirs
    character, intent(in) :: diag
    logical, intent(out) :: passed

    character :: uplo, trans
    integer :: nsize, lda_val, incx
    complex(4), dimension(n,n) :: a, a0
    complex(4), dimension(n) :: x, x0
    complex(4) :: ab(nbdirs,n,n), xb(nbdirs,n), xb_seed(nbdirs,n)
    complex(4) :: a_dir(n,n), x_dir(n)
    complex(4), dimension(n,n) :: atmp
    complex(4), dimension(n) :: xp, xm, cdiff
    real(8), parameter :: h = 1.0d-4
    real(8) :: max_error, vjp_fd, vjp_ad, abs_error, abs_reference, relerr
    logical :: has_large_errors
    integer :: nd, i, j

    uplo = 'U'; trans = 'N'
    nsize = n; lda_val = n; incx = 1

    call fill_c(a, n*n)
    a = a / real(n, 4)
    do i = 1, n
      a(i,i) = cmplx(2.0d0 + abs(real(a(i,i))), aimag(a(i,i)), kind=4)
    end do
    call fill_c(x, n)
    a0 = a; x0 = x

    call fill_c(xb, nbdirs*n)
    xb_seed = xb
    ab = (0.0d0, 0.0d0)

    write(*,*) 'Testing CTRSV (n =', n, ', nbdirs =', nbdirs, ', diag = ', diag, ')'

    call ctrsv_bv(uplo, trans, diag, nsize, a, ab, lda_val, x, xb, incx, nbdirs)

    max_error = 0.0d0
    has_large_errors = .false.
    do nd = 1, nbdirs
      call fill_c(a_dir, n*n)
      call fill_c(x_dir, n)

      atmp = a0 + h*a_dir
      xp = x0 + h * x_dir
      call ctrsv(uplo, trans, diag, nsize, atmp, lda_val, xp, incx)
      atmp = a0 - h*a_dir
      xm = x0 - h * x_dir
      call ctrsv(uplo, trans, diag, nsize, atmp, lda_val, xm, incx)
      cdiff = (xp - xm) / (2.0d0 * h)

      vjp_fd = 0.0d0
      do i = 1, n
        vjp_fd = vjp_fd + real(conjg(xb_seed(nd,i)) * cdiff(i))
      end do

      vjp_ad = 0.0d0
      do j = 1, n
        do i = 1, n
          vjp_ad = vjp_ad + real(conjg(a_dir(i,j)) * ab(nd,i,j))
        end do
      end do
      do i = 1, n
        vjp_ad = vjp_ad + real(conjg(x_dir(i)) * xb(nd,i))
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

end program test_ctrsv_vector_reverse
