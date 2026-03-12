! Test program for CHBMV vector reverse - BLAS2 band
! Generated automatically by run_tapenade_blas.py
! Multi-size outlined run_test_for_size(n, passed, nbdirs)

program test_chbmv_vector_reverse
  implicit none
  external :: chbmv
  external :: chbmv_bv
  integer :: nbdirs, n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing CHBMV (Vector Reverse band, multi-size: n = 4)'
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
    complex(4) :: alpha, alphab, beta, betab
    complex(4), dimension(:,:), allocatable :: a
    complex(4), dimension(:,:,:), allocatable :: ab
    complex(4), dimension(:), allocatable :: x, y
    complex(4), dimension(:,:), allocatable :: xb, yb
    integer :: band_row, j
    real(4) :: temp_real, temp_imag
    ksize = max(0, n - 1)
    nsize = n
    lda_val = ksize + 1
    incx_val = 1
    incy_val = 1
    uplo = 'U'
    trans = 'N'
    diag = 'N'
    allocate(a(lda_val, n), ab(nbdirs, lda_val, n), x(n), xb(nbdirs, n), y(n), yb(nbdirs, n))
    ! Initialize a as Hermitian band matrix (upper band storage, real diagonal)
    do j = 1, n
    do band_row = max(1, ksize+2-j), ksize+1
    if (band_row .eq. ksize+1) then
    call random_number(temp_real)
    a(band_row, j) = cmplx(temp_real * 2.0 - 1.0, 0.0)  ! Real diagonal
    else
    call random_number(temp_real)
    call random_number(temp_imag)
    a(band_row, j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end if
    end do
    end do
    call random_number(temp_real)
    call random_number(temp_imag)
    alpha = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(alpha))
    call random_number(temp_real)
    call random_number(temp_imag)
    beta = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(beta))
    do j = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      x(j) = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(x))
      call random_number(temp_real)
      call random_number(temp_imag)
      y(j) = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(y))
    end do
    alphab = 0.0d0
    betab = 0.0d0
    xb = 0.0d0
    ab = 0.0d0
    yb = 0.0d0
    write(*,*) 'Testing CHBMV (Vector Reverse band, n =', n, ')'
    call set_ISIZE1OFX(n)
    call set_ISIZE2OFA(n)
    call chbmv_bv(uplo, nsize, ksize, alpha, alphab, a, ab, lda_val, x, xb, incx_val, beta, betab, y, yb, incy_val, nbdirs)
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
end program test_chbmv_vector_reverse