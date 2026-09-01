! Test program for ZTRSV forward (tangent) mode differentiation
! Hand-written driver following the structure of test_zgemv.f90.
! COMPLEX*16. Multi-size, sweeps DIAG in {'N','U'}.
! (UPLO='U', TRANS='N' held fixed for now.)  TRANS='N' => holomorphic,
! so the complex tangent is compared to the full complex central difference.

program test_ztrsv
  implicit none

  external :: ztrsv
  external :: ztrsv_d

  integer :: n_test, seed_array(33), test_sizes(3), i, id
  logical :: passed, all_passed
  character :: diag

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing ZTRSV (forward mode)'
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
    complex(8), dimension(n,n) :: a, a_d
    complex(8), dimension(n) :: x, x_d
    complex(8), dimension(n,n) :: a_orig, a_d_orig
    complex(8), dimension(n) :: x_orig, x_d_orig
    integer :: i

    uplo = 'U'; trans = 'N'
    nsize = n; lda_val = n; incx = 1

    call fill_c(a, n*n)
    a = a / real(n, 8)
    do i = 1, n
      a(i,i) = cmplx(2.0d0 + abs(real(a(i,i))), aimag(a(i,i)), kind=8)
    end do
    call fill_c(x, n)

    call fill_c(a_d, n*n)
    call fill_c(x_d, n)
    if (diag == 'U' .or. diag == 'u') then
      do i = 1, n
        a_d(i,i) = (0.0d0, 0.0d0)
      end do
    end if

    a_orig = a; x_orig = x; a_d_orig = a_d; x_d_orig = x_d

    write(*,*) 'Testing ZTRSV (n =', n, ', diag = ', diag, ')'

    call ztrsv_d(uplo, trans, diag, nsize, a, a_d, lda_val, x, x_d, incx)

    call check_derivatives_numerically(n, uplo, trans, diag, nsize, lda_val, &
         incx, a_orig, x_orig, a_d_orig, x_d_orig, x_d, passed)
  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, uplo, trans, diag, nsize, &
       lda_val, incx, a_orig, x_orig, a_d_orig, x_d_orig, x_d, passed)
    integer, intent(in) :: n, nsize, lda_val, incx
    character, intent(in) :: uplo, trans, diag
    complex(8), intent(in) :: a_orig(n,n), a_d_orig(n,n)
    complex(8), intent(in) :: x_orig(n), x_d_orig(n), x_d(n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0d-6
    real(8) :: max_error, abs_error, abs_reference, error_bound
    complex(8) :: central_diff
    logical :: has_large_errors
    complex(8), dimension(n,n) :: a
    complex(8), dimension(n) :: x, x_forward, x_backward
    integer :: i

    max_error = 0.0d0
    has_large_errors = .false.

    a = a_orig + h * a_d_orig
    x = x_orig + h * x_d_orig
    call ztrsv(uplo, trans, diag, nsize, a, lda_val, x, incx)
    x_forward = x

    a = a_orig - h * a_d_orig
    x = x_orig - h * x_d_orig
    call ztrsv(uplo, trans, diag, nsize, a, lda_val, x, incx)
    x_backward = x

    do i = 1, n
      central_diff = (x_forward(i) - x_backward(i)) / (2.0d0 * h)
      abs_error = abs(central_diff - x_d(i))
      abs_reference = abs(x_d(i))
      error_bound = 1.0d-5 + 1.0d-5 * abs_reference
      if (abs_error > error_bound) then
        has_large_errors = .true.
        write(*,*) 'Large error in X(', i, '):', central_diff, x_d(i)
      end if
      max_error = max(max_error, abs_error / max(abs_reference, 1.0d-10))
    end do

    write(*,*) 'Maximum relative error:', max_error
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
  end subroutine check_derivatives_numerically

end program test_ztrsv
