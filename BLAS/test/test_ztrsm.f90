! Test program for ZTRSM forward (tangent) mode differentiation
! Hand-written driver following the structure of test_zgemv.f90.
! COMPLEX*16. Multi-size, sweeps DIAG in {'N','U'}.
! (SIDE='L', UPLO='U', TRANSA='N' held fixed for now.)
! For TRANSA='N' the map is holomorphic, so the complex tangent is
! compared against the full complex central difference.

program test_ztrsm
  implicit none

  external :: ztrsm
  external :: ztrsm_d

  integer :: n_test, seed_array(33), test_sizes(3), i, id
  logical :: passed, all_passed
  character :: diag

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing ZTRSM (forward mode)'
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

  subroutine fill_c1(z)
    complex(8), intent(out) :: z
    real(8) :: r, s
    call random_number(r)
    call random_number(s)
    z = cmplx(2.0d0*r - 1.0d0, 2.0d0*s - 1.0d0, kind=8)
  end subroutine fill_c1

  subroutine run_test_for_size(n, diag, passed)
    integer, intent(in) :: n
    character, intent(in) :: diag
    logical, intent(out) :: passed

    character :: side, uplo, transa
    integer :: msize, nsize, lda_val, ldb_val
    complex(8) :: alpha, alpha_d
    complex(8), dimension(n,n) :: a, b, a_d, b_d
    complex(8) :: alpha_orig, alpha_d_orig
    complex(8), dimension(n,n) :: a_orig, a_d_orig, b_orig, b_d_orig
    integer :: i

    side = 'L'; uplo = 'U'; transa = 'N'
    msize = n; nsize = n; lda_val = n; ldb_val = n

    call fill_c1(alpha)
    call fill_c(a, n*n)
    a = a / real(n, 8)
    do i = 1, n
      a(i,i) = cmplx(2.0d0 + abs(real(a(i,i))), aimag(a(i,i)), kind=8)
    end do
    call fill_c(b, n*n)

    call fill_c1(alpha_d)
    call fill_c(a_d, n*n)
    call fill_c(b_d, n*n)
    if (diag == 'U' .or. diag == 'u') then
      do i = 1, n
        a_d(i,i) = (0.0d0, 0.0d0)
      end do
    end if

    alpha_orig = alpha; a_orig = a; b_orig = b
    alpha_d_orig = alpha_d; a_d_orig = a_d; b_d_orig = b_d

    write(*,*) 'Testing ZTRSM (n =', n, ', diag = ', diag, ')'

    call ztrsm_d(side, uplo, transa, diag, msize, nsize, alpha, alpha_d, &
                 a, a_d, lda_val, b, b_d, ldb_val)

    call check_derivatives_numerically(n, side, uplo, transa, diag, msize, &
         nsize, lda_val, ldb_val, alpha_orig, a_orig, b_orig, alpha_d_orig, &
         a_d_orig, b_d_orig, b_d, passed)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, side, uplo, transa, diag, &
       msize, nsize, lda_val, ldb_val, alpha_orig, a_orig, b_orig, &
       alpha_d_orig, a_d_orig, b_d_orig, b_d, passed)
    integer, intent(in) :: n, msize, nsize, lda_val, ldb_val
    character, intent(in) :: side, uplo, transa, diag
    complex(8), intent(in) :: alpha_orig, alpha_d_orig
    complex(8), intent(in) :: a_orig(n,n), a_d_orig(n,n)
    complex(8), intent(in) :: b_orig(n,n), b_d_orig(n,n), b_d(n,n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0d-6
    real(8) :: max_error, abs_error, abs_reference, error_bound
    complex(8) :: central_diff
    logical :: has_large_errors
    complex(8), dimension(n,n) :: b_forward, b_backward, a, b
    complex(8) :: alpha
    integer :: i, j

    max_error = 0.0d0
    has_large_errors = .false.

    alpha = alpha_orig + h * alpha_d_orig
    a = a_orig + h * a_d_orig
    b = b_orig + h * b_d_orig
    call ztrsm(side, uplo, transa, diag, msize, nsize, alpha, a, lda_val, b, ldb_val)
    b_forward = b

    alpha = alpha_orig - h * alpha_d_orig
    a = a_orig - h * a_d_orig
    b = b_orig - h * b_d_orig
    call ztrsm(side, uplo, transa, diag, msize, nsize, alpha, a, lda_val, b, ldb_val)
    b_backward = b

    do j = 1, n
      do i = 1, n
        central_diff = (b_forward(i,j) - b_backward(i,j)) / (2.0d0 * h)
        abs_error = abs(central_diff - b_d(i,j))
        abs_reference = abs(b_d(i,j))
        error_bound = 1.0d-5 + 1.0d-5 * abs_reference
        if (abs_error > error_bound) then
          has_large_errors = .true.
          write(*,*) 'Large error in B(', i, ',', j, '):', central_diff, b_d(i,j)
        end if
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
  end subroutine check_derivatives_numerically

end program test_ztrsm
