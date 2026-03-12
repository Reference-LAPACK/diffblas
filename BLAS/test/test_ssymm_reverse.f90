! Test program for SSYMM reverse (BLAS3 outlined)
program test_ssymm_reverse
  implicit none
  external :: ssymm
  external :: ssymm_b
  integer :: n_test, test_sizes(1), i
  integer :: seed_array(33)
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing SSYMM (multi-size: n =', test_sizes(1), ')'
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
    real(4) :: alpha, alphab, beta, betab
    real(4), dimension(n,n) :: a, ab, b, bb, c, cb
    real(4), dimension(n,n) :: cb_seed, c_plus, c_minus
    real(4), dimension(n,n) :: c_orig
    real(4) :: alpha_dir, beta_dir
    real(4), dimension(n,n) :: a_dir, b_dir, c_dir, a_fd, b_fd
    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_fd, vjp_ad, abs_error, ref_c, relative_error, abs_reference
    real(4) :: vjp_ad_alpha, vjp_ad_beta, vjp_ad_a, vjp_ad_b, vjp_ad_c
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
    c_orig = c
    ! Seed direction on output (C or B) for VJP; then zero input adjoints
    call random_number(cb)
    cb = cb * 2.0d0 - 1.0d0
    cb_seed = cb
    write(*,*) 'Testing SSYMM (n =', n, ')'
    alphab = 0.0d0
    betab = 0.0d0
    ab = 0.0d0
    bb = 0.0d0
    call set_ISIZE2OFA(n)
    call set_ISIZE2OFB(n)
    call ssymm_b(side, uplo, msize, nsize, alpha, alphab, a, ab, lda_val, b, bb, ldb_val, beta, betab, c, cb, ldc_val)
    call set_ISIZE2OFA(-1)
    call set_ISIZE2OFB(-1)
    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    ! VJP finite-difference check: perturb inputs by (alphab, ab, bb, betab), compare d(cb_seed*C)/d_dir
    call random_number(tr)
    alpha_dir = tr * 2.0d0 - 1.0d0
    call random_number(tr)
    beta_dir = tr * 2.0d0 - 1.0d0
    call random_number(a_dir)
    a_dir = a_dir * 2.0d0 - 1.0d0
    do jj = 1, n
      do ii = 1, n
        if (ii > jj) a_dir(ii,jj) = 0.0d0
      end do
    end do
    call random_number(b_dir)
    b_dir = b_dir * 2.0d0 - 1.0d0
    call random_number(c_dir)
    c_dir = c_dir * 2.0d0 - 1.0d0
    a_fd = a + h * a_dir
    b_fd = b + h * b_dir
    c_plus = c_orig + h * c_dir
    call ssymm(side, uplo, msize, nsize, alpha + h*alpha_dir, a_fd, lda_val, b_fd, ldb_val, beta + h*beta_dir, c_plus, ldc_val)
    a_fd = a - h * a_dir
    b_fd = b - h * b_dir
    c_minus = c_orig - h * c_dir
    call ssymm(side, uplo, msize, nsize, alpha - h*alpha_dir, a_fd, lda_val, b_fd, ldb_val, beta - h*beta_dir, c_minus, ldc_val)
    vjp_fd = 0.0d0
    do jj = 1, n
      do ii = 1, n
        vjp_fd = vjp_fd + cb_seed(ii,jj) * (c_plus(ii,jj) - c_minus(ii,jj)) / (2.0d0 * h)
      end do
    end do
    vjp_ad = 0.0d0
    vjp_ad_alpha = 0.0d0
    vjp_ad_beta  = 0.0d0
    vjp_ad_a     = 0.0d0
    vjp_ad_b     = 0.0d0
    vjp_ad_c     = 0.0d0
    vjp_ad_alpha = alpha_dir * alphab
    vjp_ad_beta  = beta_dir * betab
    vjp_ad = vjp_ad + vjp_ad_alpha + vjp_ad_beta
    do jj = 1, n
      do ii = 1, n
        if (ii <= jj) then
          vjp_ad_a = vjp_ad_a + a_dir(ii,jj) * ab(ii,jj)
        end if
      end do
    end do
    vjp_ad_b = sum(b_dir * bb)
    vjp_ad_c = sum(c_dir * cb)
    vjp_ad = vjp_ad + vjp_ad_a + vjp_ad_b + vjp_ad_c
    write(*,*) 'VJP components: fd=', vjp_fd, ' ad=', vjp_ad
    write(*,*) '  ad_alpha=', vjp_ad_alpha, ' ad_beta=', vjp_ad_beta
    write(*,*) '  ad_A=', vjp_ad_a, ' ad_B=', vjp_ad_b, ' ad_C=', vjp_ad_c
    abs_error = abs(vjp_fd - vjp_ad)
    abs_reference = abs(vjp_ad)
    if (abs_reference > 1.0e-10) then
      relative_error = abs_error / abs_reference
    else
      relative_error = abs_error
    end if
    ref_c = abs(vjp_ad) + 1.0d0
    passed = (abs_error <= 1.0e-3 * ref_c)
    write(*,*) ''
    write(*,*) 'Maximum relative error:', relative_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    if (.not. passed) write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    if (passed) write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
  end subroutine run_test_for_size
end program test_ssymm_reverse