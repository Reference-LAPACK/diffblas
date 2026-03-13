! Test program for STRSM vector forward (BLAS3 outlined)
program test_strsm_vector_forward
  implicit none
  external :: strsm
  external :: strsm_dv
  integer :: nbdirs, n_test, test_sizes(1), i
  integer :: seed_array(33)
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing STRSM (Vector Forward, multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    nbdirs = n_test
    call run_test_for_size(n_test, passed, nbdirs)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) write(*,*) 'PASS: All sizes completed successfully'
  if (.not. all_passed) write(*,*) 'FAIL: One or more sizes had derivative errors'
contains
  subroutine run_test_for_size(n, passed, nbdirs)
    integer, intent(in) :: n, nbdirs
    logical, intent(out) :: passed
    integer :: msize, nsize, ksize, lda_val, ldb_val, ldc_val
    character :: side, uplo, transa
    character :: diag
    real(4) :: alpha, beta
    real(4), dimension(n,n) :: a, b, c
    real(4), dimension(nbdirs) :: alpha_dv, beta_dv
    real(4), dimension(nbdirs,n,n) :: a_dv, b_dv, c_dv
    real(4), dimension(nbdirs,n,n) :: b_dv_seed
    real(4), dimension(n,n) :: b_orig, b_plus, b_minus
    real(4), dimension(n,n) :: a_t, b_t
    real(4), parameter :: h = 1.0e-3
    real(4) :: max_err, abs_err, ref_c, max_err_over_dirs, worst_ref_c, relative_error
    integer :: ii, jj, idir, k
    real(4) :: tr, ti
    msize = n
    nsize = n
    ksize = n
    lda_val = n
    ldb_val = n
    ldc_val = n
    side = 'L'
    uplo = 'L'
    transa = 'N'
    diag = 'N'
    write(*,*) 'Testing STRSM (Vector Forward, n =', n, ')'
    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(beta)
    beta = beta * 2.0d0 - 1.0d0
    call random_number(a)
    a = a * 2.0d0 - 1.0d0
    call random_number(b)
    b = b * 2.0d0 - 1.0d0
    call random_number(alpha_dv)
    alpha_dv = alpha_dv * 2.0d0 - 1.0d0
    call random_number(beta_dv)
    beta_dv = beta_dv * 2.0d0 - 1.0d0
    call random_number(a_dv)
    a_dv = a_dv * 2.0d0 - 1.0d0
    call random_number(b_dv)
    b_dv = b_dv * 2.0d0 - 1.0d0
    b_orig = b
    b_dv_seed = b_dv
    call strsm_dv(side, uplo, transa, diag, msize, nsize, alpha, alpha_dv, a, a_dv, lda_val, b, b_dv, ldb_val, nbdirs)
    write(*,*) 'Function calls completed successfully'
    ! Finite-difference check per direction k: (output(+h) - output(-h))/(2h) vs c_dv(k,:,:)
    passed = .true.
    max_err_over_dirs = 0.0d0
    worst_ref_c = 1.0d0
    do k = 1, nbdirs
      a_t = a + h * a_dv(k,:,:)
      b_plus = b_orig + h * b_dv_seed(k,:,:)
      call strsm(side, uplo, transa, diag, msize, nsize, alpha + h*alpha_dv(k), a_t, lda_val, b_plus, ldb_val)
      a_t = a - h * a_dv(k,:,:)
      b_minus = b_orig - h * b_dv_seed(k,:,:)
      call strsm(side, uplo, transa, diag, msize, nsize, alpha - h*alpha_dv(k), a_t, lda_val, b_minus, ldb_val)
      max_err = 0.0d0
      do jj = 1, n
        do ii = 1, n
          abs_err = abs((b_plus(ii,jj) - b_minus(ii,jj)) / (2.0d0 * h) - b_dv(k,ii,jj))
          if (abs_err > max_err) max_err = abs_err
        end do
      end do
      ref_c = maxval(abs(b_dv(k,:,:))) + 1.0d0
      if (max_err > 1.0e-3 * ref_c) then
        passed = .false.
        write(*,*) '  direction k=', k, ' max_err=', max_err, ' ref_c=', ref_c, ' tol=', (1.0e-3)*ref_c
      end if
      if (max_err > max_err_over_dirs) then
        max_err_over_dirs = max_err
        worst_ref_c = ref_c
      end if
    end do
    relative_error = 0.0d0
    if (worst_ref_c > 1.0d-10) relative_error = max_err_over_dirs / worst_ref_c
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    write(*,*) 'Maximum relative error:', relative_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    if (.not. passed) write(*,*) 'FAIL: Derivatives are outside tolerance'
    if (passed) write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
  end subroutine run_test_for_size
end program test_strsm_vector_forward