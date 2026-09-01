! Test program for DTRSM reverse (adjoint) mode differentiation
! Hand-written driver following the structure of test_dgemv_reverse.f90.
! Validates the adjoint via the dot-product (VJP) identity:
!   <bb_seed, dX> = <a_dir,ab> + <b_dir,bb_out> + alpha_dir*alphab
! Uses REAL*8. Multi-size, sweeps DIAG in {'N','U'}.
! (SIDE='L', UPLO='U', TRANSA='N' held fixed for now.)

program test_dtrsm_reverse
  implicit none

  external :: dtrsm
  external :: dtrsm_b

  integer :: n_test, seed_array(33), test_sizes(3), i, id
  logical :: passed, all_passed
  character :: diag

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing DTRSM (reverse mode)'
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

    character :: side, uplo, transa
    integer :: msize, nsize, lda_val, ldb_val
    real(8) :: alpha, alphab
    real(8), dimension(n,n) :: a, b, ab, bb
    real(8) :: alpha_orig
    real(8), dimension(n,n) :: a_orig, b_orig, bb_orig
    integer :: i, j

    side = 'L'
    uplo = 'U'
    transa = 'N'
    msize = n
    nsize = n
    lda_val = n
    ldb_val = n

    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(a)
    a = (a * 2.0d0 - 1.0d0) / real(n, 8)
    do i = 1, n
      a(i,i) = 2.0d0 + abs(a(i,i))
    end do
    call random_number(b)
    b = b * 2.0d0 - 1.0d0

    alpha_orig = alpha
    a_orig = a
    b_orig = b

    ! Output-space seed dF/dX
    call random_number(bb)
    bb = bb * 2.0d0 - 1.0d0
    bb_orig = bb

    alphab = 0.0d0
    ab = 0.0d0

    write(*,*) 'Testing DTRSM (n =', n, ', diag = ', diag, ')'

    call dtrsm_b(side, uplo, transa, diag, msize, nsize, alpha, alphab, &
                 a, ab, lda_val, b, bb, ldb_val)

    call check_vjp_numerically(n, side, uplo, transa, diag, msize, nsize, &
         lda_val, ldb_val, alpha_orig, a_orig, b_orig, bb_orig, alphab, ab, &
         bb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, side, uplo, transa, diag, msize, nsize, &
       lda_val, ldb_val, alpha_orig, a_orig, b_orig, bb_orig, alphab, ab, &
       bb, passed)
    implicit none
    integer, intent(in) :: n, msize, nsize, lda_val, ldb_val
    character, intent(in) :: side, uplo, transa, diag
    real(8), intent(in) :: alpha_orig
    real(8), intent(in) :: a_orig(n,n), b_orig(n,n), bb_orig(n,n)
    real(8), intent(in) :: alphab, ab(n,n), bb(n,n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0d-7
    real(8) :: vjp_ad, vjp_fd, relative_error, max_error
    real(8) :: abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j

    real(8) :: alpha_dir, alpha
    real(8), dimension(n,n) :: a_dir, b_dir, a, b
    real(8), dimension(n,n) :: x_plus, x_minus, x_cdiff

    has_large_errors = .false.

    ! Random perturbation directions for the inputs
    call random_number(alpha_dir)
    alpha_dir = alpha_dir * 2.0d0 - 1.0d0
    call random_number(a_dir)
    a_dir = a_dir * 2.0d0 - 1.0d0
    call random_number(b_dir)
    b_dir = b_dir * 2.0d0 - 1.0d0

    ! Central difference of the solution X w.r.t. the input direction
    alpha = alpha_orig + h * alpha_dir
    a = a_orig + h * a_dir
    b = b_orig + h * b_dir
    call dtrsm(side, uplo, transa, diag, msize, nsize, alpha, a, lda_val, b, ldb_val)
    x_plus = b

    alpha = alpha_orig - h * alpha_dir
    a = a_orig - h * a_dir
    b = b_orig - h * b_dir
    call dtrsm(side, uplo, transa, diag, msize, nsize, alpha, a, lda_val, b, ldb_val)
    x_minus = b

    x_cdiff = (x_plus - x_minus) / (2.0d0 * h)

    ! <bb_seed, dX>
    vjp_fd = 0.0d0
    do j = 1, n
      do i = 1, n
        vjp_fd = vjp_fd + bb_orig(i,j) * x_cdiff(i,j)
      end do
    end do

    ! <a_dir,ab> + <b_dir,bb_out> + alpha_dir*alphab
    vjp_ad = alpha_dir * alphab
    do j = 1, n
      do i = 1, n
        vjp_ad = vjp_ad + a_dir(i,j) * ab(i,j) + b_dir(i,j) * bb(i,j)
      end do
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

end program test_dtrsm_reverse
