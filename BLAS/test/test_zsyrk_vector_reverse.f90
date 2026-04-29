! Test program for ZSYRK vector reverse (BLAS3 outlined)
program test_zsyrk_vector_reverse
  implicit none
  external :: zsyrk
  external :: zsyrk_bv
  integer :: nbdirs, n_test, test_sizes(3), i
  integer :: seed_array(33)
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing ZSYRK (Vector Reverse, multi-size: n =', test_sizes(1), ')'
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
    complex(8), dimension(n,n) :: c_plus, c_minus
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
    alphab = 0.0d0
    betab = 0.0d0
    ab = 0.0d0
    call set_ISIZE2OFA(n)
    call zsyrk_bv(uplo, transa, nsize, ksize, alpha, alphab, a, ab, lda_val, beta, betab, c, cb, ldc_val, nbdirs)
    call set_ISIZE2OFA(-1)
    write(*,*) 'Testing ZSYRK (Vector Reverse, n =', n, ')'
    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    ! VJP finite-difference check per direction k
    passed = .true.
    max_error = 0.0d0
    do k = 1, nbdirs
      a_t = a + h * ab(k,:,:)
      c_plus = c
      call zsyrk(uplo, transa, nsize, ksize, alpha + h*alphab(k), a_t, lda_val, beta + h*betab(k), c_plus, ldc_val)
      a_t = a - h * ab(k,:,:)
      c_minus = c
      call zsyrk(uplo, transa, nsize, ksize, alpha - h*alphab(k), a_t, lda_val, beta - h*betab(k), c_minus, ldc_val)
      vjp_fd = 0.0d0
      do jj = 1, n
        do ii = 1, n
          vjp_fd = vjp_fd + real(conjg(cb_seed(k,ii,jj)) * (c_plus(ii,jj) - c_minus(ii,jj)) / (2.0d0 * h))
        end do
      end do
      vjp_ad = real(conjg(alphab(k))*alphab(k)) + real(conjg(betab(k))*betab(k)) + sum(real(conjg(ab(k,:,:))*ab(k,:,:)))
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
end program test_zsyrk_vector_reverse