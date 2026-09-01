! Test program for DTRSV forward (tangent) mode differentiation
! Hand-written driver following the structure of test_dgemv.f90.
! Uses REAL*8. Multi-size, sweeps DIAG in {'N','U'}.
! (UPLO='U', TRANS='N' held fixed for now.)

program test_dtrsv
  implicit none

  external :: dtrsv
  external :: dtrsv_d

  integer :: n_test, seed_array(33), test_sizes(3), i, id
  logical :: passed, all_passed
  character :: diag

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing DTRSV (forward mode)'
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
    real(8), dimension(n,n) :: a, a_d
    real(8), dimension(n) :: x, x_d
    real(8), dimension(n,n) :: a_orig, a_d_orig
    real(8), dimension(n) :: x_orig, x_d_orig
    integer :: i, j

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

    call random_number(a_d)
    a_d = a_d * 2.0d0 - 1.0d0
    call random_number(x_d)
    x_d = x_d * 2.0d0 - 1.0d0
    if (diag == 'U' .or. diag == 'u') then
      do i = 1, n
        a_d(i,i) = 0.0d0
      end do
    end if

    a_orig = a
    x_orig = x
    a_d_orig = a_d
    x_d_orig = x_d

    write(*,*) 'Testing DTRSV (n =', n, ', diag = ', diag, ')'

    call dtrsv_d(uplo, trans, diag, nsize, a, a_d, lda_val, x, x_d, incx)

    call check_derivatives_numerically(n, uplo, trans, diag, nsize, lda_val, &
         incx, a_orig, x_orig, a_d_orig, x_d_orig, x_d, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, uplo, trans, diag, nsize, &
       lda_val, incx, a_orig, x_orig, a_d_orig, x_d_orig, x_d, passed)
    implicit none
    integer, intent(in) :: n, nsize, lda_val, incx
    character, intent(in) :: uplo, trans, diag
    real(8), intent(in) :: a_orig(n,n), a_d_orig(n,n)
    real(8), intent(in) :: x_orig(n), x_d_orig(n), x_d(n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0d-6
    real(8) :: relative_error, max_error, abs_error, abs_reference, error_bound
    real(8) :: central_diff, ad_result
    logical :: has_large_errors
    real(8), dimension(n,n) :: a
    real(8), dimension(n) :: x, x_forward, x_backward
    integer :: i

    max_error = 0.0d0
    has_large_errors = .false.

    a = a_orig + h * a_d_orig
    x = x_orig + h * x_d_orig
    call dtrsv(uplo, trans, diag, nsize, a, lda_val, x, incx)
    x_forward = x

    a = a_orig - h * a_d_orig
    x = x_orig - h * x_d_orig
    call dtrsv(uplo, trans, diag, nsize, a, lda_val, x, incx)
    x_backward = x

    do i = 1, n
      central_diff = (x_forward(i) - x_backward(i)) / (2.0d0 * h)
      ad_result = x_d(i)
      abs_error = abs(central_diff - ad_result)
      abs_reference = abs(ad_result)
      error_bound = 1.0d-5 + 1.0d-5 * abs_reference
      if (abs_error > error_bound) then
        has_large_errors = .true.
        write(*,*) 'Large error in X(', i, '):', central_diff, ad_result
      end if
      relative_error = abs_error / max(abs_reference, 1.0d-10)
      max_error = max(max_error, relative_error)
    end do

    write(*,*) 'Maximum relative error:', max_error
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_dtrsv
