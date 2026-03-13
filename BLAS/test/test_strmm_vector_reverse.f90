! Test program for STRMM vector reverse (BLAS3 outlined)
program test_strmm_vector_reverse
  implicit none
  external :: strmm
  external :: strmm_bv
  integer :: nbdirs, n_test, test_sizes(1), i
  integer :: seed_array(33)
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing STRMM (Vector Reverse, multi-size: n =', test_sizes(1), ')'
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
    real(4), dimension(nbdirs) :: alphab, betab
    real(4), dimension(nbdirs,n,n) :: ab, bb, cb
    real(4), dimension(nbdirs,n,n) :: bb_seed
    real(4), dimension(n,n) :: b_orig, b_plus, b_minus
    real(4) :: alpha_dir
    real(4), dimension(n,n) :: a_dir, b_dir, a_fd
    real(4), dimension(n,n) :: a_t, b_t
    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_fd, vjp_ad, abs_error, ref_c, relative_error, abs_reference, max_error
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
    diag = 'N'
    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(beta)
    beta = beta * 2.0d0 - 1.0d0
    call random_number(a)
    a = a * 2.0d0 - 1.0d0
    call random_number(b)
    b = b * 2.0d0 - 1.0d0
    call random_number(bb)
    bb = bb * 2.0d0 - 1.0d0
    b_orig = b
    bb_seed = bb
    alphab = 0.0d0
    betab = 0.0d0
    ab = 0.0d0
    call set_ISIZE2OFA(n)
    call strmm_bv(side, uplo, transa, diag, msize, nsize, alpha, alphab, a, ab, lda_val, b, bb, ldb_val, nbdirs)
    call set_ISIZE2OFA(-1)
    write(*,*) 'Testing STRMM (Vector Reverse, n =', n, ')'
    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    ! VJP finite-difference check per direction k
    passed = .true.
    max_error = 0.0d0
    do k = 1, nbdirs
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
      call strmm(side, uplo, transa, diag, msize, nsize, alpha + h*alpha_dir, a_fd, lda_val, b_plus, ldb_val)
      a_fd = a - h * a_dir
      b_minus = b_orig - h * b_dir
      call strmm(side, uplo, transa, diag, msize, nsize, alpha - h*alpha_dir, a_fd, lda_val, b_minus, ldb_val)
      vjp_fd = 0.0d0
      vjp_fd = sum(bb_seed(k,:,:) * (b_plus - b_minus)) / (2.0d0 * h)
      vjp_ad = alpha_dir * alphab(k) + sum(a_dir * ab(k,:,:)) + sum(b_dir * bb(k,:,:))
      abs_error = abs(vjp_fd - vjp_ad)
      abs_reference = abs(vjp_ad)
      if (abs_reference > 1.0e-10) then
        relative_error = abs_error / abs_reference
      else
        relative_error = abs_error
      end if
      if (relative_error > max_error) max_error = relative_error
      ref_c = abs(vjp_ad) + 1.0d0
      if (abs_error > 1.0e-3 * ref_c) passed = .false.
    end do
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    if (.not. passed) write(*,*) 'FAIL: Derivatives are outside tolerance'
    if (passed) write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
  end subroutine run_test_for_size
end program test_strmm_vector_reverse