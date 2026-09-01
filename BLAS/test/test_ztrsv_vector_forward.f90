! Test program for ZTRSV vector forward (tangent) mode differentiation
! Hand-written driver following the structure of test_zgemv_vector_forward.f90.
! COMPLEX*16, nbdirs directions (runtime). Sweeps DIAG in {'N','U'}.
! (UPLO='U', TRANS='N' held fixed for now.)

program test_ztrsv_vector_forward
  implicit none

  external :: ztrsv
  external :: ztrsv_dv

  integer :: n_test, seed_array(33), test_sizes(3), i, id
  logical :: passed, all_passed
  character :: diag

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing ZTRSV (vector forward mode)'
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
    complex(8), intent(out) :: z(k)
    integer :: t
    real(8) :: r, s
    do t = 1, k
      call random_number(r)
      call random_number(s)
      z(t) = cmplx(2.0d0*r - 1.0d0, 2.0d0*s - 1.0d0, kind=8)
    end do
  end subroutine fill_c

  subroutine run_test_for_size(n, nbdirs, diag, passed)
    integer, intent(in) :: n, nbdirs
    character, intent(in) :: diag
    logical, intent(out) :: passed

    character :: uplo, trans
    integer :: nsize, lda_val, incx
    complex(8), dimension(n,n) :: a, a0
    complex(8), dimension(n) :: x, x0
    complex(8) :: a_dv(nbdirs,n,n), x_dv(nbdirs,n)
    complex(8) :: a_dir(nbdirs,n,n), x_dir(nbdirs,n)
    complex(8), dimension(n) :: xp, xm, cdiff
    complex(8), dimension(n,n) :: atmp
    real(8) :: max_error, abs_error, abs_reference
    logical :: has_large_errors
    integer :: nd, i

    uplo = 'U'; trans = 'N'
    nsize = n; lda_val = n; incx = 1

    call fill_c(a, n*n)
    a = a / real(n, 8)
    do i = 1, n
      a(i,i) = cmplx(2.0d0 + abs(real(a(i,i))), aimag(a(i,i)), kind=8)
    end do
    call fill_c(x, n)
    a0 = a; x0 = x

    call fill_c(a_dv, nbdirs*n*n)
    call fill_c(x_dv, nbdirs*n)
    if (diag == 'U' .or. diag == 'u') then
      do nd = 1, nbdirs
        do i = 1, n
          a_dv(nd,i,i) = (0.0d0, 0.0d0)
        end do
      end do
    end if
    a_dir = a_dv; x_dir = x_dv

    write(*,*) 'Testing ZTRSV (n =', n, ', nbdirs =', nbdirs, ', diag = ', diag, ')'

    call ztrsv_dv(uplo, trans, diag, nsize, a, a_dv, lda_val, x, x_dv, incx, nbdirs)

    max_error = 0.0d0
    has_large_errors = .false.
    do nd = 1, nbdirs
      atmp = a0 + 1.0d-6*a_dir(nd,:,:)
      xp = x0 + 1.0d-6 * x_dir(nd,:)
      call ztrsv(uplo, trans, diag, nsize, atmp, lda_val, xp, incx)
      atmp = a0 - 1.0d-6*a_dir(nd,:,:)
      xm = x0 - 1.0d-6 * x_dir(nd,:)
      call ztrsv(uplo, trans, diag, nsize, atmp, lda_val, xm, incx)
      cdiff = (xp - xm) / (2.0d0 * 1.0d-6)
      do i = 1, n
        abs_error = abs(cdiff(i) - x_dv(nd,i))
        abs_reference = abs(x_dv(nd,i))
        if (abs_error > 1.0d-5 + 1.0d-5 * abs_reference) has_large_errors = .true.
        max_error = max(max_error, abs_error / max(abs_reference, 1.0d-10))
      end do
    end do

    write(*,*) 'Maximum relative error:', max_error
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
  end subroutine run_test_for_size

end program test_ztrsv_vector_forward
