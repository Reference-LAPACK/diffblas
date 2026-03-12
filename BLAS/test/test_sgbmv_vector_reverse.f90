! Test program for SGBMV vector reverse - BLAS2 band
! Generated automatically by run_tapenade_blas.py
! Multi-size outlined run_test_for_size(n, passed, nbdirs)

program test_sgbmv_vector_reverse
  implicit none
  external :: sgbmv
  external :: sgbmv_bv
  integer :: nbdirs, n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing SGBMV (Vector Reverse band, multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 1
    n_test = test_sizes(i)
    nbdirs = test_sizes(i)
    call run_test_for_size(n_test, passed, nbdirs)
    all_passed = all_passed .and. passed
  end do
  if (all_passed) write(*,*) 'PASS: Vector reverse band - all sizes OK'
  if (.not. all_passed) write(*,*) 'FAIL: Vector reverse band - errors'
contains
  subroutine run_test_for_size(n, passed, nbdirs)
    implicit none
    integer, intent(in) :: n, nbdirs
    logical, intent(out) :: passed
    character :: uplo, trans, diag
    integer :: nsize, ksize, lda_val, incx_val, incy_val
    integer :: msize, kl, ku
    real(4) :: alpha, alphab, beta, betab
    real(4), dimension(:,:), allocatable :: a
    real(4), dimension(:,:,:), allocatable :: ab
    real(4), dimension(:), allocatable :: x, y
    real(4), dimension(:,:), allocatable :: xb, yb
    integer :: band_row, j
    real(4) :: temp_real
    ksize = max(0, n - 1)
    msize = n
    nsize = n
    kl = 1
    ku = 1
    lda_val = kl + ku + 1
    incx_val = 1
    incy_val = 1
    uplo = 'U'
    trans = 'N'
    diag = 'N'
    allocate(a(lda_val, n), ab(nbdirs, lda_val, n), x(n), xb(nbdirs, n), y(n), yb(nbdirs, n))
    ! Initialize a as general band matrix (kl, ku band storage)
    do j = 1, n
    do band_row = max(1, ku+2-j), min(kl+ku+1, ku+msize-j+1)
    call random_number(temp_real)
    a(band_row, j) = temp_real * 2.0 - 1.0
    end do
    end do
    call random_number(alpha)
    alpha = alpha * 2.0d0 - 1.0d0
    call random_number(beta)
    beta = beta * 2.0d0 - 1.0d0
    call random_number(x)
    x = x * 2.0d0 - 1.0d0
    call random_number(y)
    y = y * 2.0d0 - 1.0d0
    alphab = 0.0d0
    betab = 0.0d0
    xb = 0.0d0
    ab = 0.0d0
    yb = 0.0d0
    write(*,*) 'Testing SGBMV (Vector Reverse band, n =', n, ')'
    call set_ISIZE1OFX(n)
    call set_ISIZE2OFA(n)
    call sgbmv_bv(trans, msize, nsize, kl, ku, alpha, alphab, a, ab, lda_val, x, xb, incx_val, beta, betab, y, yb, incy_val, nbdirs)
    call set_ISIZE1OFX(-1)
    call set_ISIZE2OFA(-1)
    passed = .true.
    if (allocated(a)) deallocate(a)
    if (allocated(ab)) deallocate(ab)
    if (allocated(x)) deallocate(x)
    if (allocated(xb)) deallocate(xb)
    if (allocated(y)) deallocate(y)
    if (allocated(yb)) deallocate(yb)
  end subroutine run_test_for_size
end program test_sgbmv_vector_reverse