! Test program for DTRSM forward (tangent) mode differentiation
! Hand-written driver following the structure of test_dgemv.f90.
! Uses REAL*8 precision. Multi-size, and sweeps DIAG in {'N','U'}.
! (SIDE='L', UPLO='U', TRANSA='N' held fixed for now.)

program test_dtrsm
  implicit none

  external :: dtrsm
  external :: dtrsm_d

  integer :: n_test, seed_array(33), test_sizes(3), i, id
  logical :: passed, all_passed
  character :: diag

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing DTRSM (forward mode)'
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
    real(8) :: alpha
    real(8), dimension(n,n) :: a, b
    ! Derivative variables
    real(8) :: alpha_d
    real(8), dimension(n,n) :: a_d, b_d
    ! Restoration / stored directions
    real(8) :: alpha_orig, alpha_d_orig
    real(8), dimension(n,n) :: a_orig, a_d_orig, b_orig, b_d_orig
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

    ! Well-conditioned triangular A: small off-diagonals, dominant diagonal
    call random_number(a)
    a = (a * 2.0d0 - 1.0d0) / real(n, 8)
    do i = 1, n
      a(i,i) = 2.0d0 + abs(a(i,i))
    end do

    call random_number(b)
    b = b * 2.0d0 - 1.0d0

    ! Input derivative directions
    call random_number(alpha_d)
    alpha_d = alpha_d * 2.0d0 - 1.0d0
    call random_number(a_d)
    a_d = a_d * 2.0d0 - 1.0d0
    call random_number(b_d)
    b_d = b_d * 2.0d0 - 1.0d0

    ! For unit-triangular A the diagonal is constant, so its derivative is 0
    if (diag == 'U' .or. diag == 'u') then
      do i = 1, n
        a_d(i,i) = 0.0d0
      end do
    end if

    alpha_orig = alpha
    a_orig = a
    b_orig = b
    alpha_d_orig = alpha_d
    a_d_orig = a_d
    b_d_orig = b_d

    write(*,*) 'Testing DTRSM (n =', n, ', diag = ', diag, ')'

    call dtrsm_d(side, uplo, transa, diag, msize, nsize, alpha, alpha_d, &
                 a, a_d, lda_val, b, b_d, ldb_val)

    call check_derivatives_numerically(n, side, uplo, transa, diag, msize, &
         nsize, lda_val, ldb_val, alpha_orig, a_orig, b_orig, alpha_d_orig, &
         a_d_orig, b_d_orig, b_d, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, side, uplo, transa, diag, &
       msize, nsize, lda_val, ldb_val, alpha_orig, a_orig, b_orig, &
       alpha_d_orig, a_d_orig, b_d_orig, b_d, passed)
    implicit none
    integer, intent(in) :: n, msize, nsize, lda_val, ldb_val
    character, intent(in) :: side, uplo, transa, diag
    real(8), intent(in) :: alpha_orig, alpha_d_orig
    real(8), intent(in) :: a_orig(n,n), a_d_orig(n,n)
    real(8), intent(in) :: b_orig(n,n), b_d_orig(n,n)
    real(8), intent(in) :: b_d(n,n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0d-6
    real(8) :: relative_error, max_error, abs_error, abs_reference, error_bound
    real(8) :: central_diff, ad_result
    logical :: has_large_errors
    real(8), dimension(n,n) :: b_forward, b_backward, a, b
    real(8) :: alpha
    integer :: i, j

    max_error = 0.0d0
    has_large_errors = .false.

    ! Forward perturbation
    alpha = alpha_orig + h * alpha_d_orig
    a = a_orig + h * a_d_orig
    b = b_orig + h * b_d_orig
    call dtrsm(side, uplo, transa, diag, msize, nsize, alpha, a, lda_val, b, ldb_val)
    b_forward = b

    ! Backward perturbation
    alpha = alpha_orig - h * alpha_d_orig
    a = a_orig - h * a_d_orig
    b = b_orig - h * b_d_orig
    call dtrsm(side, uplo, transa, diag, msize, nsize, alpha, a, lda_val, b, ldb_val)
    b_backward = b

    do j = 1, n
      do i = 1, n
        central_diff = (b_forward(i,j) - b_backward(i,j)) / (2.0d0 * h)
        ad_result = b_d(i,j)
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 1.0d-5 + 1.0d-5 * abs_reference
        if (abs_error > error_bound) then
          has_large_errors = .true.
          write(*,*) 'Large error in B(', i, ',', j, '):', central_diff, ad_result
        end if
        relative_error = abs_error / max(abs_reference, 1.0d-10)
        max_error = max(max_error, relative_error)
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

end program test_dtrsm
