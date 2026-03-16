! Test program for CTRMV differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_ctrmv
  implicit none

  external :: ctrmv
  external :: ctrmv_d

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(3)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4, 10, 25 /)
  write(*,*) 'Testing CTRMV (multi-size: n = 4)'
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

    character :: uplo
    character :: trans
    character :: diag
    integer :: nsize
    complex(4), dimension(n,n) :: a
    integer :: lda_val
    complex(4), dimension(n) :: x
    integer :: incx

    ! Derivative variables
    complex(4), dimension(n,n) :: a_d
    complex(4), dimension(n) :: x_d

    ! Array restoration and derivative storage
    complex(4), dimension(n,n) :: a_orig, a_d_orig
    complex(4), dimension(n) :: x_orig, x_d_orig
    real(4) :: temp_re, temp_im  ! For complex random init
    integer :: i, j

    uplo = 'U'
    trans = 'N'
    diag = 'N'
    nsize = n
    lda_val = n
    incx = 1

    call random_number(temp_re)
    call random_number(temp_im)
    a = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      x(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do

    ! Initialize input derivatives
    call random_number(temp_re)
    call random_number(temp_im)
    a_d = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      x_d(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do

    ! Store _orig and _d_orig
    a_d_orig = a_d
    x_d_orig = x_d
    a_orig = a
    x_orig = x

    write(*,*) 'Testing CTRMV (n =', n, ')'
    x_orig = x

    ! Call the differentiated function
    call ctrmv_d(uplo, trans, diag, nsize, a, a_d, lda_val, x, x_d, 1)
    a_d = a_d_orig

    write(*,*) 'Function calls completed successfully'

    ! Numerical differentiation check
    call check_derivatives_numerically(n, trans, uplo, diag, nsize, lda_val, x_orig, a_orig, x_d_orig, a_d_orig, x_d, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, trans, uplo, diag, nsize, lda_val, x_orig, a_orig, x_d_orig, a_d_orig, x_d, passed)
    implicit none
    integer, intent(in) :: n
    character, intent(in) :: trans
    character, intent(in) :: uplo
    character, intent(in) :: diag
    integer, intent(in) :: nsize
    integer, intent(in) :: lda_val
    complex(4), intent(in) :: x_orig(n), x_d_orig(n)
    complex(4), intent(in) :: a_orig(n,n), a_d_orig(n,n)
    complex(4), intent(in) :: x_d(n)
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3  ! Step size for finite differences
    real(4) :: relative_error, max_error
    real(4) :: abs_error, abs_reference, error_bound
    real(4) :: central_diff, ad_result
    logical :: has_large_errors
    complex(4), dimension(n) :: x_forward, x_backward
    integer :: i, j
    complex(4), dimension(n) :: x
    complex(4), dimension(n,n) :: a

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    ! Forward perturbation: f(x + h)
    x = x_orig + h * x_d_orig
    a = a_orig + h * a_d_orig
    call ctrmv(uplo, trans, diag, nsize, a, lda_val, x, 1)
    x_forward = x

    ! Backward perturbation: f(x - h)
    x = x_orig - h * x_d_orig
    a = a_orig - h * a_d_orig
    call ctrmv(uplo, trans, diag, nsize, a, lda_val, x, 1)
    x_backward = x

    ! Compute central differences and compare with AD results
    do i = 1, n
        central_diff = (x_forward(i) - x_backward(i)) / (2.0e0 * h)
        ad_result = x_d(i)
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 1.0e-3 + 1.0e-3 * abs_reference
        if (abs_error > error_bound) then
          has_large_errors = .true.
          relative_error = abs_error / max(abs_reference, 1.0e-10)
          write(*,*) 'Large error in output X(', i, '):'
          write(*,*) '  Central diff: ', central_diff
          write(*,*) '  AD result:   ', ad_result
          write(*,*) '  Absolute error:', abs_error
          write(*,*) '  Error bound:', error_bound
          write(*,*) '  Relative error:', relative_error
        end if
        relative_error = abs_error / max(abs_reference, 1.0e-10)
        max_error = max(max_error, relative_error)
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

end program test_ctrmv