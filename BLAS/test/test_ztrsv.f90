! Test program for ZTRSV differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_ztrsv
  implicit none

  external :: ztrsv
  external :: ztrsv_d

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing ZTRSV (multi-size: n = 4)'
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

    character :: uplo
    character :: trans
    character :: diag
    integer :: nsize
    complex(8), dimension(n,n) :: a
    integer :: lda_val
    complex(8), dimension(n) :: x
    integer :: incx

    ! Derivative variables
    complex(8), dimension(n,n) :: a_d
    complex(8), dimension(n) :: x_d

    ! Array restoration and derivative storage
    complex(8), dimension(n,n) :: a_orig, a_d_orig
    complex(8), dimension(n) :: x_orig, x_d_orig
    real(8) :: temp_re, temp_im  ! For complex random init
    integer :: i, j

    uplo = 'U'
    trans = 'N'
    diag = 'N'
    nsize = n
    lda_val = n
    incx = 1

    call random_number(temp_re)
    call random_number(temp_im)
    a = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      x(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    end do

    ! Initialize input derivatives
    call random_number(temp_re)
    call random_number(temp_im)
    a_d = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      x_d(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    end do

    ! Store _orig and _d_orig
    a_d_orig = a_d
    x_d_orig = x_d
    a_orig = a
    x_orig = x

    write(*,*) 'Testing ZTRSV (n =', n, ')'
    x_orig = x

    ! Call the differentiated function
    call ztrsv_d(uplo, trans, diag, nsize, a, a_d, lda_val, x, x_d, 1)

    write(*,*) 'Function calls completed successfully'

    ! Numerical differentiation check
    call check_derivatives_numerically(n, trans, uplo, diag, nsize, lda_val, a_orig, x_orig, a_d_orig, x_d_orig, x_d, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, trans, uplo, diag, nsize, lda_val, a_orig, x_orig, a_d_orig, x_d_orig, x_d, passed)
    implicit none
    integer, intent(in) :: n
    character, intent(in) :: trans
    character, intent(in) :: uplo
    character, intent(in) :: diag
    integer, intent(in) :: nsize
    integer, intent(in) :: lda_val
    complex(8), intent(in) :: a_orig(n,n), a_d_orig(n,n)
    complex(8), intent(in) :: x_orig(n), x_d_orig(n)
    complex(8), intent(in) :: x_d(n)
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-6  ! Step size for finite differences
    real(8) :: relative_error, max_error
    real(8) :: abs_error, abs_reference, error_bound
    real(8) :: central_diff, ad_result
    logical :: has_large_errors
    complex(8), dimension(n) :: x_forward, x_backward
    integer :: i, j
    complex(8), dimension(n,n) :: a
    complex(8), dimension(n) :: x

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    ! Forward perturbation: f(x + h)
    a = a_orig + h * a_d_orig
    x = x_orig + h * x_d_orig
    call ztrsv(uplo, trans, diag, nsize, a, lda_val, x, 1)
    x_forward = x

    ! Backward perturbation: f(x - h)
    a = a_orig - h * a_d_orig
    x = x_orig - h * x_d_orig
    call ztrsv(uplo, trans, diag, nsize, a, lda_val, x, 1)
    x_backward = x

    ! Compute central differences and compare with AD results
    do i = 1, n
        central_diff = (x_forward(i) - x_backward(i)) / (2.0e0 * h)
        ad_result = x_d(i)
        abs_error = abs(central_diff - ad_result)
        abs_reference = abs(ad_result)
        error_bound = 1.0e-5 + 1.0e-5 * abs_reference
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
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_ztrsv