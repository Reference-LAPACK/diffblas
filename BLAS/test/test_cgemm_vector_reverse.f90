! Test program for CGEMM vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision with nbdirs=n
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_cgemm_vector_reverse
  implicit none

  external :: cgemm
  external :: cgemm_bv

  integer :: nbdirs
  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(3)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing CGEMM (Vector Reverse, multi-size: n =', test_sizes(1), ')'
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
    complex(4), dimension(nbdirs) :: alphab, betab
    complex(4), dimension(nbdirs,n,n) :: ab, bb, cb
    complex(4) :: alpha_orig, beta_orig
    complex(4), dimension(n,n) :: a_orig, b_orig, c_orig
    complex(4), dimension(nbdirs,n,n) :: cb_orig
    integer :: k, ii, jj
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

    alpha_orig = alpha
    a_orig = a
    b_orig = b
    beta_orig = beta
    c_orig = c

    do k = 1, nbdirs
      do jj = 1, n
        do ii = 1, n
          call random_number(temp_real)
          call random_number(temp_imag)
          cb(k,ii,jj) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(cb))
        end do
      end do
    end do
    cb_orig = cb

    alphab = 0.0d0
    ab = 0.0d0
    bb = 0.0d0
    betab = 0.0d0

    write(*,*) 'Testing CGEMM (Vector Reverse, n =', n, ')'

    call set_ISIZE2OFA(n)
    call set_ISIZE2OFB(n)

    call cgemm_bv(transa, transb, msize, nsize, ksize, alpha, alphab, a, ab, lda_val, b, bb, ldb_val, beta, betab, c, cb, ldc_val, nbdirs)

    call set_ISIZE2OFA(-1)
    call set_ISIZE2OFB(-1)

    call check_vjp_numerically(n, nbdirs, transa, transb, msize, nsize, ksize, lda_val, ldb_val, ldc_val, alpha_orig, a_orig, b_orig, beta_orig, c_orig, cb_orig, alphab, ab, bb, betab, cb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, nbdirs, transa, transb, msize, nsize, ksize, lda_val, ldb_val, ldc_val, alpha_orig, a_orig, b_orig, beta_orig, c_orig, cb_orig, alphab, ab, bb, betab, cb, passed)
    implicit none
    integer, intent(in) :: n, nbdirs
    character, intent(in) :: transa, transb
    integer, intent(in) :: msize, nsize, ksize, lda_val, ldb_val, ldc_val
    complex(4), intent(in) :: alpha_orig, beta_orig
    complex(4), intent(in) :: a_orig(n,n), b_orig(n,n), c_orig(n,n)
    complex(4), intent(in) :: cb_orig(nbdirs,n,n)
    complex(4), intent(in) :: alphab(nbdirs), betab(nbdirs)
    complex(4), intent(in) :: ab(nbdirs,n,n), bb(nbdirs,n,n), cb(nbdirs,n,n)
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3
    real(4) :: relative_error, max_error, abs_error, abs_reference, error_bound
    real(4) :: vjp_ad, vjp_fd
    complex(4) :: alpha_dir, beta_dir
    complex(4), dimension(n,n) :: a_dir, b_dir, c_dir
    complex(4), dimension(n,n) :: c_plus, c_minus, c_central_diff
    complex(4) :: alpha, beta
    complex(4), dimension(n,n) :: a, b, c
    complex(4), dimension(n*n) :: temp_products
    integer :: n_products, i, j, k
    integer :: ii, jj
    real(4) :: temp_real, temp_imag
    logical :: has_large_errors

    max_error = 0.0d0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    do k = 1, nbdirs
      call random_number(temp_real)
      call random_number(temp_imag)
      alpha_dir = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(alpha_dir))
      do jj = 1, n
        do ii = 1, n
          call random_number(temp_real)
          call random_number(temp_imag)
          a_dir(ii,jj) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(a_dir))
        end do
      end do
      do jj = 1, n
        do ii = 1, n
          call random_number(temp_real)
          call random_number(temp_imag)
          b_dir(ii,jj) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(b_dir))
        end do
      end do
      call random_number(temp_real)
      call random_number(temp_imag)
      beta_dir = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(beta_dir))
      do jj = 1, n
        do ii = 1, n
          call random_number(temp_real)
          call random_number(temp_imag)
          c_dir(ii,jj) = cmplx(temp_real*2.0 - 1.0, temp_imag*2.0 - 1.0, kind=kind(c_dir))
        end do
      end do
      alpha = alpha_orig + h * alpha_dir
      a = a_orig + h * a_dir
      b = b_orig + h * b_dir
      beta = beta_orig + h * beta_dir
      c = c_orig + h * c_dir
      call cgemm(transa, transb, msize, nsize, ksize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
      c_plus = c
      alpha = alpha_orig - h * alpha_dir
      a = a_orig - h * a_dir
      b = b_orig - h * b_dir
      beta = beta_orig - h * beta_dir
      c = c_orig - h * c_dir
      call cgemm(transa, transb, msize, nsize, ksize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
      c_minus = c
      c_central_diff = (c_plus - c_minus) / (2.0d0 * h)
      vjp_fd = 0.0d0
      n_products = 0
      do j = 1, n
        do i = 1, n
          n_products = n_products + 1
          temp_products(n_products) = conjg(cb_orig(k,i,j)) * c_central_diff(i,j)
        end do
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_fd = vjp_fd + real(temp_products(i))
      end do
      vjp_ad = 0.0d0
      n_products = 0
      do j = 1, n
        do i = 1, n
          n_products = n_products + 1
          temp_products(n_products) = conjg(b_dir(i,j)) * bb(k,i,j)
        end do
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_ad = vjp_ad + real(temp_products(i))
      end do
      vjp_ad = vjp_ad + real(conjg(beta_dir) * betab(k))
      n_products = 0
      do j = 1, n
        do i = 1, n
          n_products = n_products + 1
          temp_products(n_products) = conjg(a_dir(i,j)) * ab(k,i,j)
        end do
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_ad = vjp_ad + real(temp_products(i))
      end do
      vjp_ad = vjp_ad + real(conjg(alpha_dir) * alphab(k))
      n_products = 0
      do j = 1, n
        do i = 1, n
          n_products = n_products + 1
          temp_products(n_products) = conjg(c_dir(i,j)) * cb(k,i,j)
        end do
      end do
      call sort_array(temp_products, n_products)
      do i = 1, n_products
        vjp_ad = vjp_ad + real(temp_products(i))
      end do
      abs_error = abs(vjp_fd - vjp_ad)
      abs_reference = abs(vjp_ad)
      error_bound = 1.0e-2 + 1.0e-2 * abs_reference
      if (abs_error > error_bound) has_large_errors = .true.
      if (abs_reference > 1.0e-10) then
        relative_error = abs_error / abs_reference
      else
        relative_error = abs_error
      end if
      if (relative_error > max_error) max_error = relative_error
    end do
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-2, atol=1.0e-2'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_vjp_numerically

  subroutine sort_array(arr, n)
    implicit none
    integer, intent(in) :: n
    complex(4), dimension(n), intent(inout) :: arr
    integer :: i, j, min_idx
    complex(4) :: temp
    do i = 1, n-1
      min_idx = i
      do j = i+1, n
        if (abs(arr(j)) < abs(arr(min_idx))) min_idx = j
      end do
      if (min_idx /= i) then
        temp = arr(i)
        arr(i) = arr(min_idx)
        arr(min_idx) = temp
      end if
    end do
  end subroutine sort_array

end program test_cgemm_vector_reverse