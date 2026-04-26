! Test program for CGEMM vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=n
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_cgemm_vector_forward
  implicit none

  external :: cgemm
  external :: cgemm_dv

  integer :: nbdirs
  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(3)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing CGEMM (Vector Forward, multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 3
    n_test = test_sizes(i)
    nbdirs = test_sizes(i)
    call run_test_for_size(n_test, passed, nbdirs)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) then
    write(*,*) 'PASS: All sizes completed successfully'
  else
    write(*,*) 'FAIL: One or more sizes had derivative errors'
  end if

contains

  subroutine run_test_for_size(n, passed, nbdirs)
    implicit none
    integer, intent(in) :: n
    logical, intent(out) :: passed
    integer, intent(in) :: nbdirs

    character :: transa, transb
    integer :: msize, nsize, ksize, lda_val, ldb_val, ldc_val
    complex(4) :: alpha, beta
    complex(4), dimension(n,n) :: a, b, c
    complex(4), dimension(nbdirs) :: alpha_dv, beta_dv
    complex(4), dimension(nbdirs,n,n) :: a_dv, b_dv, c_dv
    complex(4) :: alpha_orig, beta_orig
    complex(4), dimension(nbdirs) :: alpha_dv_orig, beta_dv_orig
    complex(4), dimension(n,n) :: a_orig, b_orig, c_orig
    complex(4), dimension(nbdirs,n,n) :: a_dv_orig, b_dv_orig, c_dv_orig
    integer :: idir, ii, jj
    real(4) :: temp_real, temp_imag

    transa = 'N'
    transb = 'N'
    msize = n
    nsize = n
    ksize = n
    lda_val = n
    ldb_val = n
    ldc_val = n

    call random_number(temp_real)
    call random_number(temp_imag)
    alpha = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(alpha))
    do jj = 1, n
      do ii = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        a(ii,jj) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(a))
      end do
    end do
    do jj = 1, n
      do ii = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        b(ii,jj) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(b))
      end do
    end do
    call random_number(temp_real)
    call random_number(temp_imag)
    beta = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(beta))
    do jj = 1, n
      do ii = 1, n
        call random_number(temp_real)
        call random_number(temp_imag)
        c(ii,jj) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(c))
      end do
    end do

    do idir = 1, nbdirs
      call random_number(temp_real)
      call random_number(temp_imag)
      alpha_dv(idir) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(alpha_dv))
    end do
    do idir = 1, nbdirs
      do jj = 1, n
        do ii = 1, n
          call random_number(temp_real)
          call random_number(temp_imag)
          a_dv(idir,ii,jj) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(a_dv))
        end do
      end do
    end do
    do idir = 1, nbdirs
      do jj = 1, n
        do ii = 1, n
          call random_number(temp_real)
          call random_number(temp_imag)
          b_dv(idir,ii,jj) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(b_dv))
        end do
      end do
    end do
    do idir = 1, nbdirs
      call random_number(temp_real)
      call random_number(temp_imag)
      beta_dv(idir) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(beta_dv))
    end do
    do idir = 1, nbdirs
      do jj = 1, n
        do ii = 1, n
          call random_number(temp_real)
          call random_number(temp_imag)
          c_dv(idir,ii,jj) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(c_dv))
        end do
      end do
    end do

    alpha_orig = alpha
    alpha_dv_orig = alpha_dv
    a_orig = a
    a_dv_orig = a_dv
    b_orig = b
    b_dv_orig = b_dv
    beta_orig = beta
    beta_dv_orig = beta_dv
    c_orig = c
    c_dv_orig = c_dv

    write(*,*) 'Testing CGEMM (Vector Forward, n =', n, ')'

    call cgemm_dv(transa, transb, msize, nsize, ksize, alpha, alpha_dv, a, a_dv, lda_val, b, b_dv, ldb_val, beta, beta_dv, c, c_dv, ldc_val, nbdirs)

    write(*,*) 'Function calls completed successfully'

    call check_derivatives_numerically(n, nbdirs, transa, transb, msize, nsize, ksize, lda_val, ldb_val, ldc_val, alpha_orig, alpha_dv_orig, a_orig, a_dv_orig, b_orig, b_dv_orig, beta_orig, beta_dv_orig, c_orig, c_dv_orig, c_dv, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, nbdirs, transa, transb, msize, nsize, ksize, lda_val, ldb_val, ldc_val, alpha_orig, alpha_dv_orig, a_orig, a_dv_orig, b_orig, b_dv_orig, beta_orig, beta_dv_orig, c_orig, c_dv_orig, c_dv, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    character, intent(in) :: transa, transb
    integer, intent(in) :: msize, nsize, ksize, lda_val, ldb_val, ldc_val
    complex(4), intent(in) :: alpha_orig, beta_orig
    complex(4), intent(in) :: alpha_dv_orig(nbdirs), beta_dv_orig(nbdirs)
    complex(4), intent(in) :: a_orig(n,n), a_dv_orig(nbdirs,n,n)
    complex(4), intent(in) :: b_orig(n,n), b_dv_orig(nbdirs,n,n)
    complex(4), intent(in) :: c_orig(n,n), c_dv_orig(nbdirs,n,n)
    complex(4), intent(in) :: c_dv(nbdirs,n,n)
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3
    real(4) :: relative_error, max_error, abs_error, abs_reference, error_bound
    complex(4) :: central_diff, ad_result
    logical :: has_large_errors
    complex(4), dimension(n,n) :: c_forward, c_backward
    integer :: i, j, idir
    complex(4) :: alpha, beta
    complex(4), dimension(n,n) :: a, b, c

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    do idir = 1, nbdirs
      alpha = alpha_orig + h * alpha_dv_orig(idir)
      a = a_orig + h * a_dv_orig(idir,:,:)
      b = b_orig + h * b_dv_orig(idir,:,:)
      beta = beta_orig + h * beta_dv_orig(idir)
      c = c_orig + h * c_dv_orig(idir,:,:)
      call cgemm(transa, transb, msize, nsize, ksize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
      c_forward = c
      alpha = alpha_orig - h * alpha_dv_orig(idir)
      a = a_orig - h * a_dv_orig(idir,:,:)
      b = b_orig - h * b_dv_orig(idir,:,:)
      beta = beta_orig - h * beta_dv_orig(idir)
      c = c_orig - h * c_dv_orig(idir,:,:)
      call cgemm(transa, transb, msize, nsize, ksize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
      c_backward = c
      do j = 1, min(2, n)
        do i = 1, min(2, n)
          central_diff = (c_forward(i,j) - c_backward(i,j)) / (2.0e0 * h)
          ad_result = c_dv(idir,i,j)
          abs_error = abs(central_diff - ad_result)
          abs_reference = abs(ad_result)
          error_bound = 1.0e-3 + 1.0e-3 * abs_reference
          if (abs_error > error_bound) then
            has_large_errors = .true.
    write(*,*) '  Large error in direction', idir, ' output C(', i, ',', j, '):'
    write(*,*) '    Central diff: ', central_diff
    write(*,*) '    AD result:   ', ad_result
          end if
          relative_error = abs_error / max(abs_reference, 1.0e-10)
          max_error = max(max_error, relative_error)
        end do
      end do
    end do

    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_cgemm_vector_forward