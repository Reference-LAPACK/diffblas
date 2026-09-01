! Test program for CTRSM vector forward (tangent) mode differentiation
! Hand-written driver following the structure of test_zgemv_vector_forward.f90.
! COMPLEX*8, nbdirs directions (runtime, = matrix size). Sweeps DIAG in {'N','U'}.
! (SIDE='L', UPLO='U', TRANSA='N' held fixed for now.)

program test_ctrsm_vector_forward
  implicit none

  external :: ctrsm
  external :: ctrsm_dv

  integer :: n_test, seed_array(33), test_sizes(3), i, id
  logical :: passed, all_passed
  character :: diag

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing CTRSM (vector forward mode)'
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
    complex(4) :: alpha_dv(nbdirs)
    complex(4) :: a_dv(nbdirs,n,n), b_dv(nbdirs,n,n)
    complex(4) :: a_dir(nbdirs,n,n), b_dir(nbdirs,n,n), alpha_dir(nbdirs)
    complex(4), dimension(n,n) :: xp, xm, cdiff, atmp
    complex(4) :: altmp
    real(8) :: max_error, abs_error, abs_reference
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

    call fill_c(alpha_dv, nbdirs)
    call fill_c(a_dv, nbdirs*n*n)
    call fill_c(b_dv, nbdirs*n*n)
    if (diag == 'U' .or. diag == 'u') then
      do nd = 1, nbdirs
        do i = 1, n
          a_dv(nd,i,i) = (0.0d0, 0.0d0)
        end do
      end do
    end if
    ! keep the input directions (b_dv is overwritten with the output)
    alpha_dir = alpha_dv; a_dir = a_dv; b_dir = b_dv

    write(*,*) 'Testing CTRSM (n =', n, ', nbdirs =', nbdirs, ', diag = ', diag, ')'

    call ctrsm_dv(side, uplo, transa, diag, msize, nsize, alpha, alpha_dv, &
                  a, a_dv, lda_val, b, b_dv, ldb_val, nbdirs)

    max_error = 0.0d0
    has_large_errors = .false.
    do nd = 1, nbdirs
      altmp = alpha + 1.0d-3*alpha_dir(nd)
      atmp = a0 + 1.0d-3*a_dir(nd,:,:)
      xp = b0 + 1.0d-3 * b_dir(nd,:,:)
      call ctrsm(side, uplo, transa, diag, msize, nsize, altmp, atmp, lda_val, xp, ldb_val)
      altmp = alpha - 1.0d-3*alpha_dir(nd)
      atmp = a0 - 1.0d-3*a_dir(nd,:,:)
      xm = b0 - 1.0d-3 * b_dir(nd,:,:)
      call ctrsm(side, uplo, transa, diag, msize, nsize, altmp, atmp, lda_val, xm, ldb_val)
      cdiff = (xp - xm) / (2.0d0 * 1.0d-3)
      do j = 1, n
        do i = 1, n
          abs_error = abs(cdiff(i,j) - b_dv(nd,i,j))
          abs_reference = abs(b_dv(nd,i,j))
          if (abs_error > 1.0d-2 + 1.0d-2 * abs_reference) has_large_errors = .true.
          max_error = max(max_error, abs_error / max(abs_reference, 1.0d-6))
        end do
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

end program test_ctrsm_vector_forward
