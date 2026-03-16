! Test program for DSYMM vector reverse (BLAS3 outlined)
program test_dsymm_vector_reverse
  implicit none
  external :: dsymm
  external :: dsymm_bv
  integer :: nbdirs, n_test, test_sizes(3), i
  integer :: seed_array(33)
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing DSYMM (Vector Reverse, multi-size: n =', test_sizes(1), ')'
  all_passed = .true.
  do i = 1, 3
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
    real(8) :: alpha, beta
    real(8), dimension(n,n) :: a, b, c
    real(8), dimension(nbdirs) :: alphab, betab
    real(8), dimension(nbdirs,n,n) :: ab, bb, cb
    real(8), dimension(nbdirs,n,n) :: cb_seed
    real(8), dimension(n,n) :: c_orig, c_plus, c_minus
    real(8) :: alpha_dir, beta_dir
    real(8), dimension(n,n) :: a_dir, b_dir, c_dir
    real(8), dimension(n,n) :: a_t, b_t
    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_fd, vjp_ad, abs_error, ref_c, relative_error, abs_reference, max_error
    integer :: ii, jj, k
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
    c_orig = c
    alphab = 0.0d0
    betab = 0.0d0
    ab = 0.0d0
    bb = 0.0d0
    call set_ISIZE2OFA(n)
    call set_ISIZE2OFB(n)
    call dsymm_bv(side, uplo, msize, nsize, alpha, alphab, a, ab, lda_val, b, bb, ldb_val, beta, betab, c, cb, ldc_val, nbdirs)
    call set_ISIZE2OFA(-1)
    call set_ISIZE2OFB(-1)
    write(*,*) 'Testing DSYMM (Vector Reverse, n =', n, ')'
    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    ! VJP finite-difference check per direction k
    passed = .true.
    max_error = 0.0d0
    do k = 1, nbdirs
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
      a_t = a + h * a_dir
      b_t = b + h * b_dir
      c_plus = c_orig + h * c_dir
      call dsymm(side, uplo, msize, nsize, alpha + h*alpha_dir, a_t, lda_val, b_t, ldb_val, beta + h*beta_dir, c_plus, ldc_val)
      a_t = a - h * a_dir
      b_t = b - h * b_dir
      c_minus = c_orig - h * c_dir
      call dsymm(side, uplo, msize, nsize, alpha - h*alpha_dir, a_t, lda_val, b_t, ldb_val, beta - h*beta_dir, c_minus, ldc_val)
      vjp_fd = 0.0d0
      vjp_fd = sum(cb_seed(k,:,:) * (c_plus - c_minus)) / (2.0d0 * h)
      vjp_ad = alpha_dir * alphab(k) + beta_dir * betab(k)
      vjp_ad = vjp_ad + sum(a_dir * ab(k,:,:)) + sum(b_dir * bb(k,:,:)) + sum(c_dir * cb(k,:,:))
      abs_error = abs(vjp_fd - vjp_ad)
      abs_reference = abs(vjp_ad)
      if (abs_reference > 1.0e-10) then
        relative_error = abs_error / abs_reference
      else
        relative_error = abs_error
      end if
      if (relative_error > max_error) max_error = relative_error
      ref_c = abs(vjp_ad) + 1.0d0
      if (abs_error > 1.0e-5 * ref_c) passed = .false.
    end do
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    if (.not. passed) write(*,*) 'FAIL: Derivatives are outside tolerance'
    if (passed) write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
  end subroutine run_test_for_size
end program test_dsymm_vector_reverse