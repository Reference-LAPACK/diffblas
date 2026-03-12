! Test program for SSYR2K differentiation (BLAS3 outlined)
! Generated automatically by run_tapenade_blas.py
! Multi-size run_test_for_size(n) - BLAS3

program test_ssyr2k
  implicit none
  external :: ssyr2k
  external :: ssyr2k_d
  integer :: n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing SSYR2K (multi-size: n = 4)'
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
    real(4) :: alpha, alpha_d, beta, beta_d
    real(4), dimension(n,n) :: a, a_d, b, b_d, c, c_d
    real(4), dimension(n,n) :: c_orig, c_plus, c_minus
    real(4), parameter :: h = 1.0e-3
    real(4) :: max_err, abs_err, ref_c
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
    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(alpha_d)
    alpha_d = alpha_d * 2.0d0 - 1.0d0
    call random_number(beta)
    beta = beta * 2.0d0 - 1.0d0
    call random_number(beta_d)
    call random_number(a)
    a = a * 2.0d0 - 1.0d0
    call random_number(a_d)
    a_d = a_d * 2.0d0 - 1.0d0
    call random_number(b)
    b = b * 2.0d0 - 1.0d0
    call random_number(b_d)
    b_d = b_d * 2.0d0 - 1.0d0
    call random_number(c)
    c = c * 2.0d0 - 1.0d0
    call random_number(c_d)
    c_d = c_d * 2.0d0 - 1.0d0
    ! Set direction for derivative w.r.t. alpha only; FD check below
    alpha_d = 1.0d0
    a_d = 0.0d0
    b_d = 0.0d0
    beta_d = 0.0d0
    c_d = 0.0d0
    c_orig = c
    call ssyr2k_d(uplo, transa, nsize, ksize, alpha, alpha_d, a, a_d, lda_val, b, b_d, ldb_val, beta, beta_d, c, c_d, ldc_val)
    ! Finite-difference check: (output(alpha+h) - output(alpha-h))/(2h) vs derivative
    c_plus = c_orig
    call ssyr2k(uplo, transa, nsize, ksize, alpha + h, a, lda_val, b, ldb_val, beta, c_plus, ldc_val)
    c_minus = c_orig
    call ssyr2k(uplo, transa, nsize, ksize, alpha - h, a, lda_val, b, ldb_val, beta, c_minus, ldc_val)
    max_err = 0.0d0
    do jj = 1, n
      do ii = 1, n
        abs_err = abs((c_plus(ii,jj) - c_minus(ii,jj)) / (2.0d0 * h) - c_d(ii,jj))
        if (abs_err > max_err) max_err = abs_err
      end do
    end do
    ref_c = maxval(abs(c_d)) + 1.0d0
    passed = (max_err <= 1.0e-3 * ref_c)
    if (.not. passed) write(*,*) 'FAIL: BLAS3 scalar forward FD max_err =', max_err
    if (passed) write(*,*) 'PASS: BLAS3 scalar forward FD check'
  end subroutine run_test_for_size
end program test_ssyr2k