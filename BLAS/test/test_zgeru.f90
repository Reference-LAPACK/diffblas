! Test program for ZGERU differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_zgeru
  implicit none

  external :: zgeru
  external :: zgeru_d

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(3)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing ZGERU (multi-size: n = 4)'
  all_passed = .true.
  do i = 1, 3
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

    integer :: msize
    integer :: nsize
    complex(8) :: alpha
    complex(8), dimension(n) :: x
    integer :: incx
    complex(8), dimension(n) :: y
    integer :: incy
    complex(8), dimension(n,n) :: a
    integer :: lda_val

    ! Derivative variables
    complex(8), dimension(n) :: x_d
    complex(8), dimension(n,n) :: a_d
    complex(8) :: alpha_d
    complex(8), dimension(n) :: y_d

    ! Array restoration and derivative storage
    complex(8), dimension(n) :: x_orig, x_d_orig
    complex(8), dimension(n,n) :: a_orig, a_d_orig
    complex(8) :: alpha_orig, alpha_d_orig
    complex(8), dimension(n) :: y_orig, y_d_orig
    real(8) :: temp_re, temp_im  ! For complex random init
    integer :: i, j

    msize = n
    nsize = n
    incx = 1
    incy = 1
    lda_val = n

    call random_number(temp_re)
    call random_number(temp_im)
    alpha = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      x(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    end do
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      y(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    end do
    call random_number(temp_re)
    call random_number(temp_im)
    a = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)

    ! Initialize input derivatives
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      x_d(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    end do
    call random_number(temp_re)
    call random_number(temp_im)
    a_d = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    call random_number(temp_re)
    call random_number(temp_im)
    alpha_d = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      y_d(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    end do

    ! Store _orig and _d_orig
    x_d_orig = x_d
    a_d_orig = a_d
    alpha_d_orig = alpha_d
    y_d_orig = y_d
    x_orig = x
    a_orig = a
    alpha_orig = alpha
    y_orig = y

    write(*,*) 'Testing ZGERU (n =', n, ')'
    a_orig = a

    ! Call the differentiated function
    call zgeru_d(msize, nsize, alpha, alpha_d, x, x_d, 1, y, y_d, 1, a, a_d, lda_val)
    x_d = x_d_orig
    alpha_d = alpha_d_orig
    y_d = y_d_orig

    write(*,*) 'Function calls completed successfully'

    ! Numerical differentiation check
    call check_derivatives_numerically(n, msize, nsize, lda_val, x_orig, a_orig, alpha_orig, y_orig, x_d_orig, a_d_orig, alpha_d_orig, y_d_orig, a_d, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, msize, nsize, lda_val, x_orig, a_orig, alpha_orig, y_orig, x_d_orig, a_d_orig, alpha_d_orig, y_d_orig, a_d, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: msize
    integer, intent(in) :: nsize
    integer, intent(in) :: lda_val
    complex(8), intent(in) :: x_orig(n), x_d_orig(n)
    complex(8), intent(in) :: a_orig(n,n), a_d_orig(n,n)
    complex(8), intent(in) :: alpha_orig, alpha_d_orig
    complex(8), intent(in) :: y_orig(n), y_d_orig(n)
    complex(8), intent(in) :: a_d(n,n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-6  ! Step size for finite differences
    real(8) :: relative_error, max_error
    real(8) :: abs_error, abs_reference, error_bound
    real(8) :: central_diff, ad_result
    logical :: has_large_errors
    complex(8), dimension(n,n) :: a_forward, a_backward
    integer :: i, j
    complex(8), dimension(n) :: x
    complex(8), dimension(n,n) :: a
    complex(8) :: alpha
    complex(8), dimension(n) :: y

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    ! Forward perturbation: f(x + h)
    x = x_orig + h * x_d_orig
    a = a_orig + h * a_d_orig
    alpha = alpha_orig + h * alpha_d_orig
    y = y_orig + h * y_d_orig
    call zgeru(msize, nsize, alpha, x, 1, y, 1, a, lda_val)
    a_forward = a

    ! Backward perturbation: f(x - h)
    x = x_orig - h * x_d_orig
    a = a_orig - h * a_d_orig
    alpha = alpha_orig - h * alpha_d_orig
    y = y_orig - h * y_d_orig
    call zgeru(msize, nsize, alpha, x, 1, y, 1, a, lda_val)
    a_backward = a

    ! Compute central differences and compare with AD results
    do j = 1, min(2, n)
      do i = 1, min(2, n)
        central_diff = (a_forward(i,j) - a_backward(i,j)) / (2.0e0 * h)
        ad_result = a_d(i,j)
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 1.0e-5 + 1.0e-5 * abs_reference
        if (abs_error > error_bound) then
          has_large_errors = .true.
          relative_error = abs_error / max(abs_reference, 1.0e-10)
          write(*,*) 'Large error in output A(', i, ',', j, '):'
          write(*,*) '  Central diff: ', central_diff
          write(*,*) '  AD result:   ', ad_result
          write(*,*) '  Absolute error:', abs_error
          write(*,*) '  Error bound:', error_bound
          write(*,*) '  Relative error:', relative_error
        end if
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        max_error = max(max_error, relative_error)
      end do
    end do

    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Derivatives are outside tolerance'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_zgeru