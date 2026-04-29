! Test program for DTRMM reverse (BLAS3 outlined)
program test_dtrmm_reverse
  implicit none
  external :: dtrmm
  external :: dtrmm_b
  integer :: n_test, test_sizes(3), i
  integer :: seed_array(33)
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing DTRMM (multi-size: n =', test_sizes(1), ')'
  all_passed = .true.
  do i = 1, 3
    call run_test_for_size(test_sizes(i), passed)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) write(*,*) 'PASS: All sizes completed successfully'
  if (.not. all_passed) write(*,*) 'FAIL: One or more sizes had derivative errors'
contains
  subroutine run_test_for_size(n, passed)
    integer, intent(in) :: n
    logical, intent(out) :: passed
    integer :: msize, nsize, ksize, lda_val, ldb_val, ldc_val
    character :: side, uplo, transa
    character :: diag
    real(8) :: alpha, alphab, beta, betab
    real(8), dimension(n,n) :: a, ab, b, bb
    real(8), dimension(n,n) :: bb_seed, b_orig, b_plus, b_minus
    real(8) :: alpha_dir
    real(8), dimension(n,n) :: a_dir, b_dir, a_fd
    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_fd, vjp_ad, abs_error, ref_c, relative_error, abs_reference
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
    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(beta)
    beta = beta * 2.0d0 - 1.0d0
    call random_number(a)
    a = a * 2.0d0 - 1.0d0
    call random_number(b)
    b = b * 2.0d0 - 1.0d0
    ! Save primal inputs for VJP base point (before _b overwrites INOUT)
    b_orig = b
    ! Seed direction on output (C or B) for VJP; then zero input adjoints
    call random_number(bb)
    bb = bb * 2.0d0 - 1.0d0
    bb_seed = bb
    write(*,*) 'Testing DTRMM (n =', n, ')'
    alphab = 0.0d0
    betab = 0.0d0
    ab = 0.0d0
    call set_ISIZE2OFA(n)
    call dtrmm_b(side, uplo, transa, diag, msize, nsize, alpha, alphab, a, ab, lda_val, b, bb, ldb_val)
    call set_ISIZE2OFA(-1)
    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    ! VJP finite-difference check: perturb inputs by (alphab, ab, bb, betab), compare d(cb_seed*C)/d_dir
    call random_number(tr)
    alpha_dir = tr * 2.0d0 - 1.0d0
    call random_number(b_dir)
    b_dir = b_dir * 2.0d0 - 1.0d0
    call random_number(a_dir)
    a_dir = a_dir * 2.0d0 - 1.0d0
    do jj = 1, n
      do ii = 1, n
        if (ii > jj) a_dir(ii,jj) = 0.0d0
      end do
    end do
    a_fd = a + h * a_dir
    b_plus = b_orig + h * b_dir
    call dtrmm(side, uplo, transa, diag, msize, nsize, alpha + h*alpha_dir, a_fd, lda_val, b_plus, ldb_val)
    a_fd = a - h * a_dir
    b_minus = b_orig - h * b_dir
    call dtrmm(side, uplo, transa, diag, msize, nsize, alpha - h*alpha_dir, a_fd, lda_val, b_minus, ldb_val)
    vjp_fd = 0.0d0
    vjp_fd = sum(bb_seed * (b_plus - b_minus)) / (2.0d0 * h)
    vjp_ad = 0.0d0
    vjp_ad = vjp_ad + alpha_dir * alphab
    vjp_ad = vjp_ad + sum(a_dir * ab)
    vjp_ad = vjp_ad + sum(b_dir * bb)
    abs_error = abs(vjp_fd - vjp_ad)
    abs_reference = abs(vjp_ad)
    if (abs_reference > 1.0e-10) then
      relative_error = abs_error / abs_reference
    else
      relative_error = abs_error
    end if
    ref_c = abs(vjp_ad) + 1.0d0
    passed = (abs_error <= 1.0e-5 * ref_c)
    write(*,*) 'Maximum relative error:', relative_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    if (.not. passed) write(*,*) 'FAIL: Derivatives are outside tolerance'
    if (passed) write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
  end subroutine run_test_for_size
end program test_dtrmm_reverse