! Test program for ZSYMM differentiation (BLAS3 outlined)
! Generated automatically by run_tapenade_blas.py
! Multi-size run_test_for_size(n) - BLAS3

program test_zsymm
  implicit none
  external :: zsymm
  external :: zsymm_d
  integer :: n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing ZSYMM (multi-size: n = 4)'
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
    complex(8) :: alpha, alpha_d, beta, beta_d
    complex(8), dimension(n,n) :: a, a_d, b, b_d, c, c_d
    complex(8), dimension(n,n) :: c_orig, c_plus, c_minus
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
    do jj = 1, n
      do ii = 1, n
        call random_number(tr)
        call random_number(ti)
        c(ii,jj) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(c))
        c_d(ii,jj) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(c_d))
      end do
    end do
    do jj = 1, n
      do ii = jj+1, n
        a(ii,jj) = a(jj,ii)
        a_d(ii,jj) = a_d(jj,ii)
      end do
    end do
    ! Set direction for derivative w.r.t. alpha only; FD check below
    alpha_d = cmplx(1.0, 0.0, kind=kind(alpha_d))
    a_d = 0.0d0
    b_d = 0.0d0
    beta_d = 0.0d0
    c_d = 0.0d0
    c_orig = c
    call zsymm_d(side, uplo, msize, nsize, alpha, alpha_d, a, a_d, lda_val, b, b_d, ldb_val, beta, beta_d, c, c_d, ldc_val)
    ! Finite-difference check: (output(alpha+h) - output(alpha-h))/(2h) vs derivative
    c_plus = c_orig
    call zsymm(side, uplo, msize, nsize, alpha + h, a, lda_val, b, ldb_val, beta, c_plus, ldc_val)
    c_minus = c_orig
    call zsymm(side, uplo, msize, nsize, alpha - h, a, lda_val, b, ldb_val, beta, c_minus, ldc_val)
    max_err = 0.0d0
    do jj = 1, n
      do ii = 1, n
        abs_err = abs((c_plus(ii,jj) - c_minus(ii,jj)) / (2.0d0 * h) - c_d(ii,jj))
        if (abs_err > max_err) max_err = abs_err
      end do
    end do
    ref_c = maxval(abs(c_d)) + 1.0d0
    passed = (max_err <= 1.0e-5 * ref_c)
    if (.not. passed) write(*,*) 'FAIL: BLAS3 scalar forward FD max_err =', max_err
    if (passed) write(*,*) 'PASS: BLAS3 scalar forward FD check'
  end subroutine run_test_for_size
end program test_zsymm