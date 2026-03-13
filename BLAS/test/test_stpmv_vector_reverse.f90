! Test program for STPMV vector reverse mode differentiation
! Multi-size outlined run_test_for_size(n) - TPMV/TPSV packed triangular

program test_stpmv_vector_reverse
  implicit none
  external :: stpmv
  external :: stpmv_bv
  integer :: nbdirs, n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing STPMV (Vector Reverse, multi-size: n =', test_sizes(1), ')'
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    nbdirs = test_sizes(i)
    call run_test_for_size(n_test, passed, nbdirs)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) write(*,*) 'PASS: All sizes completed successfully'
  if (.not. all_passed) write(*,*) 'FAIL: One or more sizes had derivative errors'
contains
  subroutine run_test_for_size(n, passed, nbdirs)
    integer, intent(in) :: n, nbdirs
    logical, intent(out) :: passed
    character :: uplo, trans, diag
    integer :: nsize, incx_val, npack
    real(4), allocatable :: ap(:), x(:)
    real(4), allocatable :: apb(:,:), xb(:,:)
    real(4), allocatable :: ap_orig(:), x_orig(:), xb_orig(:,:)
    integer :: idir, ii
    real(4) :: tr, ti
    uplo = 'L'
    trans = 'N'
    diag = 'N'
    nsize = n
    incx_val = 1
    npack = (n * (n + 1)) / 2
    allocate(ap(npack), x(n), apb(nbdirs, npack), xb(nbdirs, n))
    allocate(ap_orig(npack), x_orig(n), xb_orig(nbdirs, n))
    call random_number(ap)
    ap = ap * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    do idir = 1, nbdirs
      call random_number(xb(idir,:))
      xb(idir,:) = xb(idir,:) * 2.0d0 - 1.0d0
    end do
    ap_orig = ap
    x_orig = x
    xb_orig = xb
    apb = 0.0d0
    write(*,*) 'Testing STPMV (Vector Reverse, n =', n, ')'
    call set_ISIZE1OFAp(npack)
    ! xb holds seed (direction on output x); _bv overwrites xb with adjoint
    call stpmv_bv(uplo, trans, diag, nsize, ap, apb, x, xb, incx_val, nbdirs)
    call set_ISIZE1OFAp(-1)
    write(*,*) 'Function calls completed successfully'
    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', 1.0e-3

    call check_vjp_numerically(n, npack, nbdirs, uplo, trans, diag, nsize, incx_val, ap_orig, x_orig, xb_orig, apb, xb, passed)
    if (allocated(ap)) deallocate(ap)
    if (allocated(apb)) deallocate(apb)
    if (allocated(x)) deallocate(x)
    if (allocated(xb)) deallocate(xb)
    if (allocated(ap_orig)) deallocate(ap_orig)
    if (allocated(x_orig)) deallocate(x_orig)
    if (allocated(xb_orig)) deallocate(xb_orig)
  end subroutine run_test_for_size

  subroutine check_vjp_numerically(n, npack, nbdirs, uplo, trans, diag, nsize, incx_val, ap_orig, x_orig, xb_orig, apb, xb, passed)
    implicit none
    integer, intent(in) :: n, npack, nbdirs, nsize, incx_val
    character, intent(in) :: uplo, trans, diag
    real(4), intent(in) :: ap_orig(npack), x_orig(n), xb_orig(nbdirs,n)
    real(4), intent(in) :: apb(nbdirs,npack), xb(nbdirs,n)
    logical, intent(out) :: passed
    real(4), parameter :: h = 1.0e-3
    real(4) :: vjp_ad, vjp_fd, relative_error, max_error, abs_error, abs_reference, error_bound
    real(4), allocatable :: ap(:), x(:), ap_dir(:), x_dir(:), x_plus(:), x_minus(:)
    real(4), dimension(n) :: temp_real_fd
    integer :: k, i, ii, n_products
    real(4) :: temp_real, temp_imag
    logical :: has_large_errors
    allocate(ap(npack), x(n), ap_dir(npack), x_dir(n), x_plus(n), x_minus(n))
    max_error = 0.0d0
    has_large_errors = .false.
    do k = 1, nbdirs
      call random_number(ap_dir)
      ap_dir = ap_dir * 2.0d0 - 1.0d0
      call random_number(x_dir)
      x_dir = x_dir * 2.0d0 - 1.0d0
      ap = ap_orig + h * ap_dir
      x = x_orig + h * x_dir
      call stpmv(uplo, trans, diag, nsize, ap, x, incx_val)
      x_plus = x
      ap = ap_orig - h * ap_dir
      x = x_orig - h * x_dir
      call stpmv(uplo, trans, diag, nsize, ap, x, incx_val)
      x_minus = x
      vjp_fd = 0.0e0
      n_products = n
      do i = 1, n
        temp_real_fd(i) = xb_orig(k,i) * (x_plus(i) - x_minus(i)) / (2.0e0 * h)
      end do
      call sort_array(temp_real_fd, n_products)
      do i = 1, n_products
        vjp_fd = vjp_fd + temp_real_fd(i)
      end do
      vjp_ad = 0.0d0
      do ii = 1, npack
        vjp_ad = vjp_ad + ap_dir(ii) * apb(k,ii)
      end do
      do ii = 1, n
        vjp_ad = vjp_ad + x_dir(ii) * xb(k,ii)
      end do
      abs_error = abs(vjp_fd - vjp_ad)
      abs_reference = abs(vjp_ad)
      error_bound = 1.0e-3 + 1.0e-3 * abs_reference
      if (abs_error > error_bound) has_large_errors = .true.
      if (abs_reference > 1.0e-10) then
        relative_error = abs_error / abs_reference
      else
        relative_error = abs_error
      end if
      if (relative_error > max_error) max_error = relative_error
    end do
    deallocate(ap, x, ap_dir, x_dir, x_plus, x_minus)
    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
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
    real(4), dimension(n), intent(inout) :: arr
    integer :: i, j, min_idx
    real(4) :: temp
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
end program test_stpmv_vector_reverse