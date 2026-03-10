! Test program for ZHEMM reverse mode (adjoint) differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_zhemm_reverse
  implicit none

  external :: zhemm
  external :: zhemm_b

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing ZHEMM (multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    call run_test_for_size(n_test, passed)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) then
    write(*,*) 'PASS: All sizes completed successfully'
  else
    write(*,*) 'FAIL: One or more sizes had derivative errors'
  end if

contains

  subroutine run_test_for_size(n, passed)
    implicit none
    integer, intent(in) :: n
    logical, intent(out) :: passed

    character :: side
    character :: uplo
    integer :: msize
    integer :: nsize
    complex(8) :: alpha
    complex(8), dimension(n,n) :: a
    integer :: lda_val
    complex(8), dimension(n,n) :: b
    integer :: ldb_val
    complex(8) :: beta
    complex(8), dimension(n,n) :: c
    integer :: ldc_val
    complex(8) :: alphab
    complex(8), dimension(n,n) :: ab
    complex(8), dimension(n,n) :: bb
    complex(8) :: betab
    complex(8), dimension(n,n) :: cb
    complex(8) :: alpha_orig
    complex(8), dimension(n,n) :: a_orig
    complex(8), dimension(n,n) :: b_orig
    complex(8) :: beta_orig
    complex(8), dimension(n,n) :: c_orig
    complex(8), dimension(n,n) :: cb_orig
    real(4) :: temp_re, temp_im
    integer :: i, j

    nsize = n
    msize = n
    lda_val = n
    ldb_val = n
    ldc_val = n
    side = 'L'
    uplo = 'U'

    call random_number(temp_re)
    call random_number(temp_im)
    alpha = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    ! Initialize a as Hermitian matrix
    ! Fill diagonal with real numbers
    do i = 1, n
    call random_number(temp_re)
    a(i,i) = cmplx(temp_re * 2.0 - 1.0, 0.0)  ! Real diagonal
    end do
    
    ! Fill upper triangle with complex numbers
    do i = 1, n
    do j = i+1, n
    call random_number(temp_re)
    call random_number(temp_im)
    a(i,j) = cmplx(temp_re, temp_im) * (2.0,2.0) - (1.0,1.0)
    end do
    end do
    
    ! Fill lower triangle with complex conjugates
    do i = 2, n
    do j = 1, i-1
    a(i,j) = conjg(a(j,i))  ! A(i,j) = conj(A(j,i))
    end do
    end do
    do j = 1, n
      do i = 1, n
        call random_number(temp_re)
        call random_number(temp_im)
        b(i,j) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
      end do
    end do
    call random_number(temp_re)
    call random_number(temp_im)
    beta = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    do j = 1, n
      do i = 1, n
        call random_number(temp_re)
        call random_number(temp_im)
        c(i,j) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
      end do
    end do

    alpha_orig = alpha
    a_orig = a
    b_orig = b
    beta_orig = beta
    c_orig = c

    call random_number(temp_re)
    call random_number(temp_im)
    cb = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    cb_orig = cb

    alphab = 0.0
    ab = 0.0
    bb = 0.0
    betab = 0.0

    write(*,*) 'Testing ZHEMM (n =', n, ')'

    call set_ISIZE2OFA(n)
    call set_ISIZE2OFB(n)

    call zhemm_b(side, uplo, msize, nsize, alpha, alphab, a, ab, lda_val, b, bb, ldb_val, beta, betab, c, cb, ldc_val)

    call set_ISIZE2OFA(-1)
    call set_ISIZE2OFB(-1)

    call check_vjp_numerically(n, side, uplo, msize, nsize, lda_val, ldb_val, ldc_val, alpha_orig, a_orig, b_orig, beta_orig, c_orig, cb_orig, alphab, ab, bb, betab, cb, passed)

  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, side, uplo, msize, nsize, lda_val, ldb_val, ldc_val, alpha_orig, a_orig, b_orig, beta_orig, c_orig, cb_orig, alphab, ab, bb, betab, cb, passed)
    implicit none
    integer, intent(in) :: n
    character, intent(in) :: side
    character, intent(in) :: uplo
    integer, intent(in) :: msize
    integer, intent(in) :: nsize
    integer, intent(in) :: lda_val
    integer, intent(in) :: ldb_val
    integer, intent(in) :: ldc_val
    complex(8), intent(in) :: alpha_orig
    complex(8), intent(in) :: a_orig(n,n)
    complex(8), intent(in) :: b_orig(n,n)
    complex(8), intent(in) :: beta_orig
    complex(8), intent(in) :: c_orig(n,n)
    complex(8), intent(in) :: cb_orig(n,n)
    complex(8), intent(in) :: alphab
    complex(8), intent(in) :: ab(n,n)
    complex(8), intent(in) :: bb(n,n)
    complex(8), intent(in) :: betab
    complex(8), intent(in) :: cb(n,n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    logical :: has_large_errors
    integer :: i, j, n_products
    real(8), dimension(n) :: temp_products
    real(4) :: temp_re, temp_im

    complex(8) :: alpha_dir
    complex(8), dimension(n,n) :: a_dir
    complex(8), dimension(n,n) :: b_dir
    complex(8) :: beta_dir
    complex(8), dimension(n,n) :: c_dir

    complex(8), dimension(n,n) :: c_plus, c_minus, c_central_diff

    complex(8) :: alpha
    complex(8), dimension(n,n) :: a
    complex(8), dimension(n,n) :: b
    complex(8) :: beta
    complex(8), dimension(n,n) :: c

    max_error = 0.0
    has_large_errors = .false.

    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    call random_number(temp_re)
    call random_number(temp_im)
    alpha_dir = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    do j = 1, n
      do i = 1, n
        call random_number(temp_re)
        call random_number(temp_im)
        a_dir(i,j) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
      end do
    end do
    ! Keep perturbations consistent with Hermitian a_dir and imag(diag(a_dir)) = 0
    do i = 1, n
    a_dir(i,i) = cmplx(real(a_dir(i,i)), 0.0d0)
    end do
    do j = 1, n
    do i = j+1, n
    a_dir(i,j) = conjg(a_dir(j,i))
    end do
    end do
    do j = 1, n
      do i = 1, n
        call random_number(temp_re)
        call random_number(temp_im)
        b_dir(i,j) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
      end do
    end do
    call random_number(temp_re)
    call random_number(temp_im)
    beta_dir = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    do j = 1, n
      do i = 1, n
        call random_number(temp_re)
        call random_number(temp_im)
        c_dir(i,j) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
      end do
    end do

    alpha = alpha_orig + cmplx(h, 0.0) * alpha_dir
    a = a_orig + cmplx(h, 0.0) * a_dir
    b = b_orig + cmplx(h, 0.0) * b_dir
    beta = beta_orig + cmplx(h, 0.0) * beta_dir
    c = c_orig + cmplx(h, 0.0) * c_dir
    call zhemm(side, uplo, msize, nsize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
    c_plus = c

    alpha = alpha_orig - cmplx(h, 0.0) * alpha_dir
    a = a_orig - cmplx(h, 0.0) * a_dir
    b = b_orig - cmplx(h, 0.0) * b_dir
    beta = beta_orig - cmplx(h, 0.0) * beta_dir
    c = c_orig - cmplx(h, 0.0) * c_dir
    call zhemm(side, uplo, msize, nsize, alpha, a, lda_val, b, ldb_val, beta, c, ldc_val)
    c_minus = c

    c_central_diff = (c_plus - c_minus) / (2.0 * h)

    vjp_fd = 0.0
    do j = 1, n
      do i = 1, n
        vjp_fd = vjp_fd + real(conjg(cb_orig(i,j)) * c_central_diff(i,j))
      end do
    end do

    vjp_ad = 0.0
    vjp_ad = vjp_ad + real(conjg(alpha_dir) * alphab)
    ! Hermitian A: VJP = sum over upper triangle of conjg(a_dir)*ab + a_dir*ab^T
    do j = 1, n
      do i = 1, j
        if (i .eq. j) then
          vjp_ad = vjp_ad + real(conjg(a_dir(i,j)) * ab(i,j))
        else
          vjp_ad = vjp_ad + real(conjg(a_dir(i,j)) * ab(i,j) + a_dir(i,j) * ab(j,i))
        end if
      end do
    end do
    do j = 1, n
      do i = 1, n
        vjp_ad = vjp_ad + real(conjg(b_dir(i,j)) * bb(i,j))
      end do
    end do
    vjp_ad = vjp_ad + real(conjg(beta_dir) * betab)
    do j = 1, n
      do i = 1, n
        vjp_ad = vjp_ad + real(conjg(c_dir(i,j)) * cb(i,j))
      end do
    end do

    abs_error = abs(vjp_fd - vjp_ad)
    abs_reference = abs(vjp_ad)
    error_bound = 1.0e-5 + 1.0e-5 * abs_reference
    if (abs_error > error_bound) has_large_errors = .true.
    if (abs_reference > 1.0e-10) then
      relative_error = abs_error / abs_reference
    else
      relative_error = abs_error
    end if
    max_error = relative_error

    write(*,*) ''
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_vjp_numerically

  subroutine sort_array(arr, n)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(inout) :: arr
    integer :: i, j, min_idx
    real(8) :: temp
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

end program test_zhemm_reverse