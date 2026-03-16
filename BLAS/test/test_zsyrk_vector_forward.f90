! Test program for ZSYRK vector forward (BLAS3 outlined)
program test_zsyrk_vector_forward
  implicit none
  external :: zsyrk
  external :: zsyrk_dv
  integer :: nbdirs, n_test, test_sizes(3), i
  integer :: seed_array(33)
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing ZSYRK (Vector Forward, multi-size: n = 4)'
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
    complex(8), dimension(nbdirs) :: alpha_dv, beta_dv
    complex(8), dimension(nbdirs,n,n) :: a_dv, b_dv, c_dv
    complex(8), dimension(nbdirs,n,n) :: c_dv_seed
    complex(8), dimension(n,n) :: c_orig, c_plus, c_minus
    complex(8), dimension(n,n) :: a_t, b_t
    real(8), parameter :: h = 1.0e-7
    real(8) :: max_err, abs_err, ref_c, max_err_over_dirs, worst_ref_c, relative_error
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
    write(*,*) 'Testing ZSYRK (Vector Forward, n =', n, ')'
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
        c(ii,jj) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(c))
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
          c_dv(idir,ii,jj) = cmplx(tr*2.0-1.0, ti*2.0-1.0, kind=kind(c_dv))
        end do
      end do
    end do
    c_orig = c
    c_dv_seed = c_dv
    call zsyrk_dv(uplo, transa, nsize, ksize, alpha, alpha_dv, a, a_dv, lda_val, beta, beta_dv, c, c_dv, ldc_val, nbdirs)
    write(*,*) 'Function calls completed successfully'
    ! Finite-difference check per direction k: (output(+h) - output(-h))/(2h) vs c_dv(k,:,:)
    passed = .true.
    max_err_over_dirs = 0.0d0
    worst_ref_c = 1.0d0
    do k = 1, nbdirs
      a_t = a + h * a_dv(k,:,:)
      c_plus = c_orig + h * c_dv_seed(k,:,:)
      call zsyrk(uplo, transa, nsize, ksize, alpha + h*alpha_dv(k), a_t, lda_val, beta + h*beta_dv(k), c_plus, ldc_val)
      a_t = a - h * a_dv(k,:,:)
      c_minus = c_orig - h * c_dv_seed(k,:,:)
      call zsyrk(uplo, transa, nsize, ksize, alpha - h*alpha_dv(k), a_t, lda_val, beta - h*beta_dv(k), c_minus, ldc_val)
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
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    if (.not. passed) write(*,*) 'FAIL: Derivatives are outside tolerance'
    if (passed) write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
  end subroutine run_test_for_size
end program test_zsyrk_vector_forward