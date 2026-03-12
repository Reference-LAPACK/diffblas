! Test program for ZTRSM vector forward (BLAS3 outlined)
program test_ztrsm_vector_forward
  implicit none
  external :: ztrsm
  external :: ztrsm_dv
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
    character :: diag
    complex(8) :: alpha, beta
    complex(8), dimension(n,n) :: a, b, c
    complex(8), dimension(nbdirs) :: alpha_dv, beta_dv
    complex(8), dimension(nbdirs,n,n) :: a_dv, b_dv, c_dv
    complex(8), dimension(nbdirs,n,n) :: b_dv_seed
    complex(8), dimension(n,n) :: b_orig, b_plus, b_minus
    complex(8), dimension(n,n) :: a_t, b_t
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
    diag = 'N'
    write(*,*) 'Testing ZTRSM (Vector Forward, n =', n, ')'
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
    do idir = 1, nbdirs
      call random_number(tr)
      call random_number(ti)
      alpha_dv(idir) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(alpha_dv))
      call random_number(tr)
      call random_number(ti)
      beta_dv(idir) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(beta_dv))
      do jj = 1, n
        do ii = 1, n
          call random_number(tr)
          call random_number(ti)
          a_dv(idir,ii,jj) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(a_dv))
        end do
      end do
      do jj = 1, n
        do ii = 1, n
          call random_number(tr)
          call random_number(ti)
          b_dv(idir,ii,jj) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(b_dv))
        end do
      end do
    end do
    b_orig = b
    b_dv_seed = b_dv
    call ztrsm_dv(side, uplo, transa, diag, msize, nsize, alpha, alpha_dv, a, a_dv, lda_val, b, b_dv, ldb_val, nbdirs)
    ! Finite-difference check per direction k: (output(+h) - output(-h))/(2h) vs c_dv(k,:,:)
    passed = .true.
    do k = 1, nbdirs
      a_t = a + h * a_dv(k,:,:)
      b_plus = b_orig + h * b_dv_seed(k,:,:)
      call ztrsm(side, uplo, transa, diag, msize, nsize, alpha + h*alpha_dv(k), a_t, lda_val, b_plus, ldb_val)
      a_t = a - h * a_dv(k,:,:)
      b_minus = b_orig - h * b_dv_seed(k,:,:)
      call ztrsm(side, uplo, transa, diag, msize, nsize, alpha - h*alpha_dv(k), a_t, lda_val, b_minus, ldb_val)
      max_err = 0.0d0
      do jj = 1, n
        do ii = 1, n
          abs_err = abs((b_plus(ii,jj) - b_minus(ii,jj)) / (2.0d0 * h) - b_dv(k,ii,jj))
          if (abs_err > max_err) max_err = abs_err
        end do
      end do
      ref_c = maxval(abs(b_dv(k,:,:))) + 1.0d0
      if (max_err > 1.0e-5 * ref_c) then
        passed = .false.
        write(*,*) '  direction k=', k, ' max_err=', max_err, ' ref_c=', ref_c, ' tol=', (1.0e-5)*ref_c
      end if
    end do
    if (.not. passed) write(*,*) 'FAIL: BLAS3 vector forward FD check'
    if (passed) write(*,*) 'PASS: BLAS3 vector forward FD check'
  end subroutine run_test_for_size
end program test_ztrsm_vector_forward