! Test program for ZSYMM vector reverse (BLAS3 outlined)
program test_zsymm_vector_reverse
  implicit none
  external :: zsymm
  external :: zsymm_bv
  integer :: nbdirs, n_test, test_sizes(3), i
  integer :: seed_array(33)
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing ZSYMM (Vector Reverse, multi-size: n =', test_sizes(1), ')'
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
    complex(8) :: alpha, beta
    complex(8), dimension(n,n) :: a, b, c
    complex(8), dimension(nbdirs) :: alphab, betab
    complex(8), dimension(nbdirs,n,n) :: ab, bb, cb
    complex(8), dimension(nbdirs,n,n) :: cb_seed
    complex(8), dimension(n,n) :: c_orig, c_plus, c_minus
    complex(8) :: alpha_dir, beta_dir
    complex(8), dimension(n,n) :: a_dir, b_dir, c_dir
    complex(8), dimension(n,n) :: a_t, b_t
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
    call random_number(tr)
    call random_number(ti)
    alpha = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(alpha))
    call random_number(tr)
    call random_number(ti)
    beta = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(beta))
    do jj = 1, n
      do ii = 1, n
        call random_number(tr)
        call random_number(ti)
        a(ii,jj) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(a))
      end do
    end do
    do jj = 1, n
      do ii = 1, n
        call random_number(tr)
        call random_number(ti)
        b(ii,jj) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(b))
      end do
    end do
    do jj = 1, n
      do ii = 1, n
        call random_number(tr)
        call random_number(ti)
        c(ii,jj) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(c))
      end do
    end do
    do k = 1, nbdirs
      do jj = 1, n
        do ii = 1, n
          call random_number(tr)
          call random_number(ti)
          cb(k,ii,jj) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(cb))
        end do
      end do
    end do
    cb_seed = cb
    c_orig = c
    alphab = 0.0d0
    betab = 0.0d0
    ab = 0.0d0
    bb = 0.0d0
    call set_ISIZE2OFA(n)
    call set_ISIZE2OFB(n)
    call zsymm_bv(side, uplo, msize, nsize, alpha, alphab, a, ab, lda_val, b, bb, ldb_val, beta, betab, c, cb, ldc_val, nbdirs)
    call set_ISIZE2OFA(-1)
    call set_ISIZE2OFB(-1)
    write(*,*) 'Testing ZSYMM (Vector Reverse, n =', n, ')'
    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    ! VJP finite-difference check per direction k
    passed = .true.
    max_error = 0.0d0
    do k = 1, nbdirs
      call random_number(tr)
      call random_number(ti)
      alpha_dir = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(alpha_dir))
      call random_number(tr)
      call random_number(ti)
      beta_dir = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(beta_dir))
      do ii = 1, n
        call random_number(tr)
        a_dir(ii,ii) = cmplx(tr*2.0-1.0, 0.0, kind=kind(a_dir))
      end do
      do jj = 1, n
        do ii = jj+1, n
          call random_number(tr)
          call random_number(ti)
          a_dir(jj,ii) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(a_dir))
          a_dir(ii,jj) = conjg(a_dir(jj,ii))
        end do
      end do
      do jj = 1, n
        do ii = 1, n
          call random_number(tr)
          call random_number(ti)
          b_dir(ii,jj) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(b_dir))
          call random_number(tr)
          call random_number(ti)
          c_dir(ii,jj) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(c_dir))
        end do
      end do
      a_t = a + h * a_dir
      b_t = b + h * b_dir
      c_plus = c_orig + h * c_dir
      call zsymm(side, uplo, msize, nsize, alpha + h*alpha_dir, a_t, lda_val, b_t, ldb_val, beta + h*beta_dir, c_plus, ldc_val)
      a_t = a - h * a_dir
      b_t = b - h * b_dir
      c_minus = c_orig - h * c_dir
      call zsymm(side, uplo, msize, nsize, alpha - h*alpha_dir, a_t, lda_val, b_t, ldb_val, beta - h*beta_dir, c_minus, ldc_val)
      vjp_fd = 0.0d0
      do jj = 1, n
        do ii = 1, n
          vjp_fd = vjp_fd + real(conjg(cb_seed(k,ii,jj)) * (c_plus(ii,jj) - c_minus(ii,jj)) / (2.0d0 * h))
        end do
      end do
      vjp_ad = real(conjg(alpha_dir) * alphab(k)) + real(conjg(beta_dir) * betab(k))
      vjp_ad = vjp_ad + sum(real(conjg(a_dir) * ab(k,:,:)))
      vjp_ad = vjp_ad + sum(real(conjg(b_dir) * bb(k,:,:)))
      vjp_ad = vjp_ad + sum(real(conjg(c_dir) * cb(k,:,:)))
      abs_error = abs(vjp_fd - vjp_ad)
      abs_reference = abs(vjp_ad)
      if (abs_reference > 1.0e-10) then
        relative_error = abs_error / abs_reference
      else
        relative_error = abs_error
      end if
      if (relative_error > max_error) max_error = relative_error
      ref_c = abs(vjp_ad) + 1.0d0
      if (abs_error > 1.0e-2 * ref_c) passed = .false.
    end do
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-2, atol=1.0e-2'
    if (.not. passed) write(*,*) 'FAIL: Derivatives are outside tolerance'
    if (passed) write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
  end subroutine run_test_for_size
end program test_zsymm_vector_reverse