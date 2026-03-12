! Test program for DSYRK vector forward (BLAS3 outlined)
program test_dsyrk_vector_forward
  implicit none
  external :: dsyrk
  external :: dsyrk_dv
  integer :: nbdirs, n_test, test_sizes(1), i
  integer :: seed_array(33)
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    nbdirs = n_test
    call run_test_for_size(n_test, passed, nbdirs)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) write(*,*) 'PASS: BLAS3 vector forward'
  if (.not. all_passed) write(*,*) 'FAIL: BLAS3 vector forward'
contains
  subroutine run_test_for_size(n, passed, nbdirs)
    integer, intent(in) :: n, nbdirs
    logical, intent(out) :: passed
    integer :: msize, nsize, ksize, lda_val, ldb_val, ldc_val
    character :: side, uplo, transa
    real(8) :: alpha, beta
    real(8), dimension(n,n) :: a, b, c
    real(8), dimension(nbdirs) :: alpha_dv, beta_dv
    real(8), dimension(nbdirs,n,n) :: a_dv, b_dv, c_dv
    real(8), dimension(nbdirs,n,n) :: c_dv_seed
    real(8), dimension(n,n) :: c_orig, c_plus, c_minus
    real(8), dimension(n,n) :: a_t, b_t
    real(8), parameter :: h = 1.0e-7
    real(8) :: max_err, abs_err, ref_c
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
    write(*,*) 'Testing DSYRK (Vector Forward, n =', n, ')'
    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(beta)
    beta = beta * 2.0d0 - 1.0d0
    call random_number(a)
    a = a * 2.0d0 - 1.0d0
    call random_number(c)
    c = c * 2.0d0 - 1.0d0
    call random_number(alpha_dv)
    alpha_dv = alpha_dv * 2.0d0 - 1.0d0
    call random_number(beta_dv)
    beta_dv = beta_dv * 2.0d0 - 1.0d0
    call random_number(a_dv)
    a_dv = a_dv * 2.0d0 - 1.0d0
    call random_number(c_dv)
    c_dv = c_dv * 2.0d0 - 1.0d0
    c_orig = c
    c_dv_seed = c_dv
    call dsyrk_dv(uplo, transa, nsize, ksize, alpha, alpha_dv, a, a_dv, lda_val, beta, beta_dv, c, c_dv, ldc_val, nbdirs)
    ! Finite-difference check per direction k: (output(+h) - output(-h))/(2h) vs c_dv(k,:,:)
    passed = .true.
    do k = 1, nbdirs
      a_t = a + h * a_dv(k,:,:)
      c_plus = c_orig + h * c_dv_seed(k,:,:)
      call dsyrk(uplo, transa, nsize, ksize, alpha + h*alpha_dv(k), a_t, lda_val, beta + h*beta_dv(k), c_plus, ldc_val)
      a_t = a - h * a_dv(k,:,:)
      c_minus = c_orig - h * c_dv_seed(k,:,:)
      call dsyrk(uplo, transa, nsize, ksize, alpha - h*alpha_dv(k), a_t, lda_val, beta - h*beta_dv(k), c_minus, ldc_val)
      max_err = 0.0d0
      do jj = 1, n
        do ii = 1, n
          abs_err = abs((c_plus(ii,jj) - c_minus(ii,jj)) / (2.0d0 * h) - c_dv(k,ii,jj))
          if (abs_err > max_err) max_err = abs_err
        end do
      end do
      ref_c = maxval(abs(c_dv(k,:,:))) + 1.0d0
      if (max_err > 1.0e-5 * ref_c) then
        passed = .false.
        write(*,*) '  direction k=', k, ' max_err=', max_err, ' ref_c=', ref_c, ' tol=', (1.0e-5)*ref_c
      end if
    end do
    if (.not. passed) write(*,*) 'FAIL: BLAS3 vector forward FD check'
    if (passed) write(*,*) 'PASS: BLAS3 vector forward FD check'
  end subroutine run_test_for_size
end program test_dsyrk_vector_forward