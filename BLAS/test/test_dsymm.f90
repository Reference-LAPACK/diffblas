! Test program for DSYMM differentiation (BLAS3 outlined)
! Generated automatically by run_tapenade_blas.py
! Multi-size run_test_for_size(n) - BLAS3

program test_dsymm
  implicit none
  external :: dsymm
  external :: dsymm_d
  integer :: n_test, seed_array(33), test_sizes(3), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing DSYMM (multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 3
    n_test = test_sizes(i)
    call run_test_for_size(n_test, passed)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) write(*,*) 'PASS: All sizes completed successfully'
  if (.not. all_passed) write(*,*) 'FAIL: One or more sizes had derivative errors'
contains
  subroutine run_test_for_size(n, passed)
    implicit none
    integer, intent(in) :: n
    logical, intent(out) :: passed
    integer :: msize, nsize, ksize, lda_val, ldb_val, ldc_val
    character :: side, uplo, transa
    real(8) :: alpha, alpha_d, beta, beta_d
    real(8), dimension(n,n) :: a, a_d, b, b_d, c, c_d
    real(8), dimension(n,n) :: c_orig, c_plus, c_minus
    real(8), parameter :: h = 1.0e-7
    real(8) :: max_err, abs_err, ref_c, relative_error
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
    do jj = 1, n
      do ii = jj+1, n
        a(ii,jj) = a(jj,ii)
        a_d(ii,jj) = a_d(jj,ii)
      end do
    end do
    ! Set direction for derivative w.r.t. alpha only; FD check below
    alpha_d = 1.0d0
    a_d = 0.0d0
    b_d = 0.0d0
    beta_d = 0.0d0
    c_d = 0.0d0
    c_orig = c
    call dsymm_d(side, uplo, msize, nsize, alpha, alpha_d, a, a_d, lda_val, b, b_d, ldb_val, beta, beta_d, c, c_d, ldc_val)
    write(*,*) 'Testing DSYMM (n =', n, ')'
    write(*,*) 'Function calls completed successfully'
    ! Finite-difference check: (output(alpha+h) - output(alpha-h))/(2h) vs derivative
    c_plus = c_orig
    call dsymm(side, uplo, msize, nsize, alpha + h, a, lda_val, b, ldb_val, beta, c_plus, ldc_val)
    c_minus = c_orig
    call dsymm(side, uplo, msize, nsize, alpha - h, a, lda_val, b, ldb_val, beta, c_minus, ldc_val)
    max_err = 0.0d0
    do jj = 1, n
      do ii = 1, n
        abs_err = abs((c_plus(ii,jj) - c_minus(ii,jj)) / (2.0d0 * h) - c_d(ii,jj))
        if (abs_err > max_err) max_err = abs_err
      end do
    end do
    ref_c = maxval(abs(c_d)) + 1.0d0
    relative_error = 0.0d0
    if (ref_c > 1.0d-10) relative_error = max_err / ref_c
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Maximum relative error:', relative_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = (max_err <= 1.0e-5 * ref_c)
    if (.not. passed) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if
  end subroutine run_test_for_size
end program test_dsymm