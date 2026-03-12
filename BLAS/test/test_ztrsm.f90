! Test program for ZTRSM differentiation (BLAS3 outlined)
! Generated automatically by run_tapenade_blas.py
! Multi-size run_test_for_size(n) - BLAS3

program test_ztrsm
  implicit none
  external :: ztrsm
  external :: ztrsm_d
  integer :: n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing ZTRSM (multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    call run_test_for_size(n_test, passed)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) write(*,*) 'PASS: All sizes OK'
  if (.not. all_passed) write(*,*) 'FAIL: Derivative errors'
contains
  subroutine run_test_for_size(n, passed)
    implicit none
    integer, intent(in) :: n
    logical, intent(out) :: passed
    integer :: msize, nsize, ksize, lda_val, ldb_val, ldc_val
    character :: side, uplo, transa
    character :: diag
    complex(8) :: alpha, alpha_d, beta, beta_d
    complex(8), dimension(n,n) :: a, a_d, b, b_d
    complex(8), dimension(n,n) :: b_orig, b_plus, b_minus
    real(8), parameter :: h = 1.0e-7
    real(8) :: max_err, abs_err, ref_c
    integer :: ii, jj
    real(4) :: tr, ti
    msize = n
    nsize = n
    ksize = n
    lda_val = n
    ldb_val = n
    ldc_val = n
    side = 'L'
    uplo = 'U'
    transa = 'N'
    diag = 'N'
    call random_number(tr)
    call random_number(ti)
    alpha = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(alpha))
    alpha_d = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(alpha_d))
    call random_number(tr)
    call random_number(ti)
    beta = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(beta))
    beta_d = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(beta_d))
    do jj = 1, n
      do ii = 1, n
        call random_number(tr)
        call random_number(ti)
        a(ii,jj) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(a))
        a_d(ii,jj) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(a_d))
      end do
    end do
    do jj = 1, n
      do ii = 1, n
        call random_number(tr)
        call random_number(ti)
        b(ii,jj) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(b))
        b_d(ii,jj) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(b_d))
      end do
    end do
    ! Set direction for derivative w.r.t. alpha only; FD check below
    alpha_d = cmplx(1.0, 0.0, kind=kind(alpha_d))
    a_d = 0.0d0
    b_d = 0.0d0
    b_orig = b
    call ztrsm_d(side, uplo, transa, diag, msize, nsize, alpha, alpha_d, a, a_d, lda_val, b, b_d, ldb_val)
    ! Finite-difference check: (output(alpha+h) - output(alpha-h))/(2h) vs derivative
    b_plus = b_orig
    call ztrsm(side, uplo, transa, diag, msize, nsize, alpha + h, a, lda_val, b_plus, ldb_val)
    b_minus = b_orig
    call ztrsm(side, uplo, transa, diag, msize, nsize, alpha - h, a, lda_val, b_minus, ldb_val)
    max_err = 0.0d0
    do jj = 1, n
      do ii = 1, n
        abs_err = abs((b_plus(ii,jj) - b_minus(ii,jj)) / (2.0d0 * h) - b_d(ii,jj))
        if (abs_err > max_err) max_err = abs_err
      end do
    end do
    ref_c = maxval(abs(b_d)) + 1.0d0
    passed = (max_err <= 1.0e-5 * ref_c)
    if (.not. passed) write(*,*) 'FAIL: BLAS3 scalar forward FD max_err =', max_err
    if (passed) write(*,*) 'PASS: BLAS3 scalar forward FD check'
  end subroutine run_test_for_size
end program test_ztrsm