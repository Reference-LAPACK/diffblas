! Test program for ZTBMV vector reverse - BLAS2 band
! Generated automatically by run_tapenade_blas.py
! Multi-size outlined run_test_for_size(n, passed, nbdirs)

program test_ztbmv_vector_reverse
  implicit none
  external :: ztbmv
  external :: ztbmv_bv
  integer :: nbdirs, n_test, seed_array(33), test_sizes(1), i
  logical :: passed, all_passed
  seed_array = 42
  call random_seed(put=seed_array)
  test_sizes = (/ 4 /)
  write(*,*) 'Testing ZTBMV (Vector Reverse band, multi-size: n = 4)'
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
    complex(8) :: alpha, alphab, beta, betab
    complex(8), dimension(:,:), allocatable :: a
    complex(8), dimension(:,:,:), allocatable :: ab
    complex(8), dimension(:), allocatable :: x, y
    complex(8), dimension(:,:), allocatable :: xb, yb
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
    allocate(a(lda_val, n), ab(nbdirs, lda_val, n), x(n), xb(nbdirs, n))
    ! Initialize a as triangular band matrix (upper band storage)
    do j = 1, n
    do band_row = max(1, ksize+2-j), ksize+1
    call random_number(temp_real)
    call random_number(temp_imag)
    a(band_row, j) = cmplx(temp_real, temp_imag) * (2.0,2.0) - (1.0,1.0)
    end do
    end do
    do j = 1, n
      call random_number(temp_real)
      call random_number(temp_imag)
      x(j) = cmplx(temp_real*2.0-1.0, temp_imag*2.0-1.0, kind=kind(x))
    end do
    alphab = 0.0d0
    betab = 0.0d0
    xb = 0.0d0
    ab = 0.0d0
    write(*,*) 'Testing ZTBMV (Vector Reverse band, n =', n, ')'
    call set_ISIZE2OFA(n)
    call ztbmv_bv(uplo, trans, diag, nsize, ksize, a, ab, lda_val, x, xb, incx_val, nbdirs)
    call set_ISIZE2OFA(-1)
    passed = .true.
    if (allocated(a)) deallocate(a)
    if (allocated(ab)) deallocate(ab)
    if (allocated(x)) deallocate(x)
    if (allocated(xb)) deallocate(xb)
    if (allocated(y)) deallocate(y)
    if (allocated(yb)) deallocate(yb)
  end subroutine run_test_for_size
end program test_ztbmv_vector_reverse