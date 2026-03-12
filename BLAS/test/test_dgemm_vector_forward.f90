! Test program for DGEMM vector forward mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirs=n
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_dgemm_vector_forward
  implicit none

  external :: dgemm
  external :: dgemm_dv

  integer :: nbdirs
  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing DGEMM (Vector Forward, multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    nbdirs = test_sizes(i)
    call run_test_for_size(n_test, passed, nbdirs)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) then
    write(*,*) 'PASS: Vector forward mode - all sizes completed successfully'
  else
    write(*,*) 'FAIL: Vector forward mode - one or more sizes had derivative errors'
  end if

contains

  subroutine run_test_for_size(n, passed, nbdirs)
    implicit none
    integer, intent(in) :: n
    logical, intent(out) :: passed
    integer, intent(in) :: nbdirs

    character :: transa, transb
    integer :: msize, nsize, ksize, lda_val, ldb_val, ldc_val
    real(8) :: alpha, beta
    real(8), dimension(n,n) :: a, b, c
    real(8), dimension(nbdirs) :: alpha_dv, beta_dv
    real(8), dimension(nbdirs,n,n) :: a_dv, b_dv, c_dv
    real(8) :: alpha_orig, beta_orig
    real(8), dimension(nbdirs) :: alpha_dv_orig, beta_dv_orig
    real(8), dimension(n,n) :: a_orig, b_orig, c_orig
    real(8), dimension(nbdirs,n,n) :: a_dv_orig, b_dv_orig, c_dv_orig
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

    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(a)
    a = a * 2.0d0 - 1.0d0
    call random_number(b)
    b = b * 2.0d0 - 1.0d0
    call random_number(beta)
    beta = beta * 2.0d0 - 1.0d0
    call random_number(c)
    c = c * 2.0d0 - 1.0d0

    do idir = 1, nbdirs
      call random_number(temp_real)
      alpha_dv(idir) = temp_real * 2.0d0 - 1.0d0
    end do
    do idir = 1, nbdirs
      call random_number(a_dv(idir,:,:))
      a_dv(idir,:,:) = a_dv(idir,:,:) * 2.0d0 - 1.0d0
    end do
    do idir = 1, nbdirs
      call random_number(b_dv(idir,:,:))
      b_dv(idir,:,:) = b_dv(idir,:,:) * 2.0d0 - 1.0d0
    end do
    do idir = 1, nbdirs
      call random_number(temp_real)
      beta_dv(idir) = temp_real * 2.0d0 - 1.0d0
    end do
    do idir = 1, nbdirs
      call random_number(c_dv(idir,:,:))
      c_dv(idir,:,:) = c_dv(idir,:,:) * 2.0d0 - 1.0d0
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

    write(*,*) 'Testing DGEMM (Vector Forward, n =', n, ')'

    call dgemm_dv(transa, transb, msize, nsize, ksize, alpha, alpha_dv, a, a_dv, lda_val, b, b_dv, ldb_val, beta, beta_dv, c, c_dv, ldc_val, nbdirs)

    write(*,*) 'Function calls completed successfully'

    call check_derivatives_numerically(n, nbdirs, transa, transb, msize, nsize, ksize, lda_val, ldb_val, ldc_val, alpha_orig, alpha_dv_orig, a_orig, a_dv_orig, b_orig, b_dv_orig, beta_orig, beta_dv_orig, c_orig, c_dv_orig, c_dv, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, nbdirs, transa, transb, msize, nsize, ksize, lda_val, ldb_val, ldc_val, alpha_orig, alpha_dv_orig, a_orig, a_dv_orig, b_orig, b_dv_orig, beta_orig, beta_dv_orig, c_orig, c_dv_orig, c_dv, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    character, intent(in) :: transa, transb
    integer, intent(in) :: msize, nsize, ksize, lda_val, ldb_val, ldc_val
    real(8), intent(in) :: alpha_orig, beta_orig
    real(8), intent(in) :: alpha_dv_orig(nbdirs), beta_dv_orig(nbdirs)
    real(8), intent(in) :: a_orig(n,n), a_dv_orig(nbdirs,n,n)
    real(8), intent(in) :: b_orig(n,n), b_dv_orig(nbdirs,n,n)
    real(8), intent(in) :: c_orig(n,n), c_dv_orig(nbdirs,n,n)
    real(8), intent(in) :: c_dv(nbdirs,n,n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: relative_error, max_error, abs_error, abs_reference, error_bound
    real(8) :: central_diff, ad_result
    logical :: has_large_errors
    real(8), dimension(n,n) :: c_forward, c_backward
    integer :: i, j, idir
    real(8) :: alpha, beta
    real(8), dimension(n,n) :: a, b, c

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking vector derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    do idir = 1, nbdirs
      alpha = alpha_orig + h * alpha_dv_orig(idir)
      a = a_orig + h * a_dv_orig(idir,:,:)
      b = b_orig + h * b_dv_orig(idir,:,:)
      beta = beta_orig + h * beta_dv_orig(idir)
      c = c_orig + h * c_dv_orig(idir,:,:)
      call dgemm(transa, transb, msize, nsize, ksize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
      c_forward = c
      alpha = alpha_orig - h * alpha_dv_orig(idir)
      a = a_orig - h * a_dv_orig(idir,:,:)
      b = b_orig - h * b_dv_orig(idir,:,:)
      beta = beta_orig - h * beta_dv_orig(idir)
      c = c_orig - h * c_dv_orig(idir,:,:)
      call dgemm(transa, transb, msize, nsize, ksize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
      c_backward = c
      do j = 1, min(2, n)
        do i = 1, min(2, n)
          central_diff = (c_forward(i,j) - c_backward(i,j)) / (2.0e0 * h)
          ad_result = c_dv(idir,i,j)
          abs_error = abs(central_diff - ad_result)
          abs_reference = abs(ad_result)
          error_bound = 1.0e-5 + 1.0e-5 * abs_reference
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

    write(*,*) 'Maximum relative error across all directions:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in vector derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Vector derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_dgemm_vector_forward