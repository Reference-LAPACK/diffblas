! Test program for DSYR2K reverse (BLAS3 outlined)
program test_dsyr2k_reverse
  implicit none
  external :: dsyr2k
  external :: dsyr2k_b
  integer :: n_test, test_sizes(1), i
  integer :: seed_array(33)
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing DSYR2K (multi-size: n =', test_sizes(1), ')'
  all_passed = .true.
  do i = 1, 1
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
    real(8) :: alpha, alphab, beta, betab
    real(8), dimension(n,n) :: a, ab, b, bb, c, cb
    real(8), dimension(n,n) :: cb_seed, c_plus, c_minus
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
    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(beta)
    beta = beta * 2.0d0 - 1.0d0
    call random_number(a)
    a = a * 2.0d0 - 1.0d0
    call random_number(b)
    b = b * 2.0d0 - 1.0d0
    call random_number(c)
    c = c * 2.0d0 - 1.0d0
    ! Save primal inputs for VJP base point (before _b overwrites INOUT)
    ! Seed direction on output (C or B) for VJP; then zero input adjoints
    call random_number(cb)
    cb = cb * 2.0d0 - 1.0d0
    cb_seed = cb
    write(*,*) 'Testing DSYR2K (n =', n, ')'
    alphab = 0.0d0
    betab = 0.0d0
    ab = 0.0d0
    bb = 0.0d0
    call set_ISIZE2OFA(n)
    call set_ISIZE2OFB(n)
    call dsyr2k_b(uplo, transa, nsize, ksize, alpha, alphab, a, ab, lda_val, b, bb, ldb_val, beta, betab, c, cb, ldc_val)
    call set_ISIZE2OFA(-1)
    call set_ISIZE2OFB(-1)
    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    ! VJP finite-difference check: perturb inputs by (alphab, ab, bb, betab), compare d(cb_seed*C)/d_dir
    c_plus = c
    call dsyr2k(uplo, transa, nsize, ksize, alpha + h*alphab, a + h*ab, lda_val, b + h*bb, ldb_val, beta + h*betab, c_plus, ldc_val)
    c_minus = c
    call dsyr2k(uplo, transa, nsize, ksize, alpha - h*alphab, a - h*ab, lda_val, b - h*bb, ldb_val, beta - h*betab, c_minus, ldc_val)
    vjp_fd = 0.0d0
    do jj = 1, n
      do ii = 1, n
        vjp_fd = vjp_fd + cb_seed(ii,jj) * (c_plus(ii,jj) - c_minus(ii,jj)) / (2.0d0 * h)
      end do
    end do
    vjp_ad = 0.0d0
    vjp_ad = alphab*alphab + betab*betab + sum(ab*ab)
    vjp_ad = vjp_ad + sum(bb*bb)
    abs_error = abs(vjp_fd - vjp_ad)
    abs_reference = abs(vjp_ad)
    if (abs_reference > 1.0e-10) then
      relative_error = abs_error / abs_reference
    else
      relative_error = abs_error
    end if
    ref_c = abs(vjp_ad) + 1.0d0
    passed = (abs_error <= 1.0e-5 * ref_c)
    write(*,*) ''
    write(*,*) 'Maximum relative error:', relative_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    if (.not. passed) write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    if (passed) write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
  end subroutine run_test_for_size
end program test_dsyr2k_reverse