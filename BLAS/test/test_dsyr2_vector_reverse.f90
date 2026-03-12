! Test program for DSYR2 vector reverse mode differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision with nbdirs=n
! Multi-size outlined run_test_for_size(n) - SYR/SYR2

program test_dsyr2_vector_reverse
  implicit none
  external :: dsyr2
  external :: dsyr2_bv
  integer :: nbdirs, n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing DSYR2 (Vector Reverse, multi-size: n =', test_sizes(1), ')'
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    nbdirs = test_sizes(i)
    call run_test_for_size(n_test, passed, nbdirs)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) write(*,*) 'PASS: Vector reverse - all sizes completed successfully'
  if (.not. all_passed) write(*,*) 'FAIL: Vector reverse - one or more sizes had derivative errors'
contains
  subroutine run_test_for_size(n, passed, nbdirs)
    integer, intent(in) :: n, nbdirs
    logical, intent(out) :: passed
    character :: uplo
    integer :: nsize, lda_val, incx_val, incy_val
    real(8) :: alpha
    real(8), dimension(n) :: x
    real(8), dimension(n,n) :: a
    real(8), dimension(nbdirs) :: alphab
    real(8), dimension(nbdirs,n) :: xb
    real(8), dimension(nbdirs,n,n) :: ab
    real(8), dimension(n) :: y
    real(8), dimension(nbdirs,n) :: yb
    real(8), dimension(nbdirs,n,n) :: ab_orig
    real(8) :: alpha_orig
    real(8), dimension(n) :: x_orig
    real(8), dimension(n) :: y_orig
    real(8), dimension(n,n) :: a_orig
    integer :: k, ii, jj
    real(4) :: tr, ti
    uplo = 'U'
    nsize = n
    lda_val = n
    incx_val = 1
    incy_val = 1
    call random_number(tr)
    alpha = tr * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    call random_number(y)
    y = y * 2.0d0 - 1.0d0
    call random_number(a)
    a = a * 2.0d0 - 1.0d0
    do jj = 1, n
      do ii = jj+1, n
        a(ii,jj) = a(jj,ii)
      end do
    end do
    do k = 1, nbdirs
      call random_number(ab(k,:,:))
      ab(k,:,:) = ab(k,:,:) * 2.0d0 - 1.0d0
      do jj = 1, n
        do ii = jj+1, n
          ab(k,ii,jj) = ab(k,jj,ii)
        end do
      end do
    end do
    alpha_orig = alpha
    x_orig = x
    y_orig = y
    a_orig = a
    ab_orig = ab
    alphab = 0.0d0
    xb = 0.0d0
    yb = 0.0d0
    write(*,*) 'Testing DSYR2 (Vector Reverse, n =', n, ')'
    call set_ISIZE1OFX(n)
    call set_ISIZE1OFY(n)
    call dsyr2_bv(uplo, nsize, alpha, alphab, x, xb, incx_val, y, yb, incy_val, a, ab, lda_val, nbdirs)
    call set_ISIZE1OFX(-1)
    call set_ISIZE1OFY(-1)
    call check_vjp_syr_syr2(n, nbdirs, uplo, nsize, lda_val, incx_val, incy_val, alpha_orig, x_orig, a_orig, ab_orig, alphab, xb, ab, passed, y_orig, yb)
  end subroutine run_test_for_size
  subroutine check_vjp_syr_syr2(n, nbdirs, uplo, nsize, lda_val, incx_val, incy_val, alpha, x, a, ab_orig, alphab, xb, ab, passed, y, yb)
    integer, intent(in) :: n, nbdirs
    character, intent(in) :: uplo
    integer, intent(in) :: nsize, lda_val, incx_val, incy_val
    real(8), intent(in) :: alpha, x(n)
    real(8), intent(in) :: a(n,n)
    real(8), intent(in) :: ab_orig(nbdirs,n,n)
    real(8), intent(in) :: alphab(nbdirs), xb(nbdirs,n)
    real(8), intent(in) :: ab(nbdirs,n,n)
    logical, intent(out) :: passed
    real(8), intent(in), optional :: y(n), yb(nbdirs,n)
    real(8), parameter :: h = 1.0e-7
    real(8) :: vjp_fd, vjp_ad, re, err_bnd, tr, ti, relative_error, abs_reference, max_error
    real(8) :: alpha_dir
    real(8), dimension(n,n) :: a_dir, a_t, a_plus, a_minus, a_cdiff
    real(8), dimension(n) :: x_dir, x_t
    real(8), dimension(n) :: y_dir, y_t
    integer :: k, i, j
    logical :: has_err
    has_err = .false.
    max_error = 0.0d0
    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h
    do k = 1, nbdirs
      call random_number(tr)
      call random_number(ti)
      alpha_dir = tr * 2.0d0 - 1.0d0
      call random_number(x_dir)
      x_dir = x_dir * 2.0d0 - 1.0d0
      if (present(y)) call random_number(y_dir)
      if (present(y)) y_dir = y_dir * 2.0d0 - 1.0d0
      call random_number(a_dir)
      a_dir = a_dir * 2.0d0 - 1.0d0
      do j = 1, n
        do i = j+1, n
          a_dir(i,j) = a_dir(j,i)
        end do
      end do
      a_t = a + h * a_dir
      x_t = x + h * x_dir
      if (present(y)) y_t = y + h * y_dir
      if (present(y)) then
        call dsyr2(uplo, nsize, alpha + h*alpha_dir, x_t, incx_val, y_t, incy_val, a_t, lda_val)
      else
        call dsyr2(uplo, nsize, alpha + h*alpha_dir, x_t, incx_val, a_t, lda_val)
      end if
      a_plus = a_t
      a_t = a - h * a_dir
      x_t = x - h * x_dir
      if (present(y)) y_t = y - h * y_dir
      if (present(y)) then
        call dsyr2(uplo, nsize, alpha - h*alpha_dir, x_t, incx_val, y_t, incy_val, a_t, lda_val)
      else
        call dsyr2(uplo, nsize, alpha - h*alpha_dir, x_t, incx_val, a_t, lda_val)
      end if
      a_minus = a_t
      a_cdiff = (a_plus - a_minus) / (2.0d0 * h)
      vjp_fd = 0.0d0
      do j = 1, n
        do i = 1, j
          if (i.eq.j) then
            vjp_fd = vjp_fd + ab_orig(k,i,j) * a_cdiff(i,j)
          else
            vjp_fd = vjp_fd + ab_orig(k,i,j)*(a_cdiff(i,j)+a_cdiff(j,i))
          end if
        end do
      end do
      vjp_ad = alpha_dir * alphab(k)
      vjp_ad = vjp_ad + sum(x_dir*xb(k,:))
      do j = 1, n
        do i = 1, j
          if (i.eq.j) then
            vjp_ad = vjp_ad + a_dir(i,j)*ab(k,i,j)
          else
            vjp_ad = vjp_ad + a_dir(i,j)*(ab(k,i,j)+ab(k,j,i))
          end if
        end do
      end do
      if (present(y)) then
        vjp_ad = vjp_ad + sum(y_dir*yb(k,:))
      end if
      re = abs(vjp_fd - vjp_ad)
      abs_reference = abs(vjp_ad)
      if (abs_reference > 1.0e-10) then
        relative_error = re / abs_reference
      else
        relative_error = re
      end if
      if (relative_error > max_error) max_error = relative_error
      err_bnd = 1.0e-5 + 1.0e-5 * abs(vjp_ad)
      if (re > err_bnd) has_err = .true.
    end do
    write(*,*) ''
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_err
    if (.not. passed) write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    if (passed) write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
  end subroutine check_vjp_syr_syr2
end program test_dsyr2_vector_reverse