! Test program for CDOTU differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*4 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_cdotu
  implicit none

  complex(4), external :: cdotu
  complex(4), external :: cdotu_d

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing CDOTU (multi-size: n = 4)'
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

    integer :: nsize
    complex(4), dimension(n) :: cx
    integer :: incx
    complex(4), dimension(n) :: cy
    integer :: incy

    ! Derivative variables
    complex(4), dimension(n) :: cy_d
    complex(4), dimension(n) :: cx_d
    complex(4) :: cdotu_d_result  ! Derivative of function result (avoid name clash with func_d)

    ! Array restoration and derivative storage
    complex(4), dimension(n) :: cy_orig, cy_d_orig
    complex(4), dimension(n) :: cx_orig, cx_d_orig
    complex(4) :: cdotu_orig  ! Function result (no _d_orig - use _d_result)
    real(4) :: temp_re, temp_im  ! For complex random init
    integer :: i, j

    nsize = n
    incx = 1
    incy = 1

    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      cx(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      cy(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do

    ! Initialize input derivatives
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      cy_d(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      cx_d(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=4)
    end do

    ! Store _orig and _d_orig
    cy_d_orig = cy_d
    cx_d_orig = cx_d
    cy_orig = cy
    cx_orig = cx
    cdotu_orig = cdotu(nsize, cx, 1, cy, 1)

    write(*,*) 'Testing CDOTU (n =', n, ')'

    ! Call the differentiated function
    cdotu_d_result = cdotu_d(nsize, cx, cx_d, 1, cy, cy_d, 1, cdotu_orig)

    write(*,*) 'Function calls completed successfully'

    ! Numerical differentiation check
    call check_derivatives_numerically(n, nsize, cy_orig, cx_orig, cdotu_orig, cy_d_orig, cx_d_orig, cdotu_d_result, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, nsize, cy_orig, cx_orig, cdotu_orig, cy_d_orig, cx_d_orig, cdotu_d_result, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nsize
    complex(4), intent(in) :: cy_orig(n), cy_d_orig(n)
    complex(4), intent(in) :: cx_orig(n), cx_d_orig(n)
    complex(4), intent(in) :: cdotu_orig
    complex(4), intent(in) :: cdotu_d_result
    logical, intent(out) :: passed

    real(4), parameter :: h = 1.0e-3  ! Step size for finite differences
    real(4) :: relative_error, max_error
    real(4) :: abs_error, abs_reference, error_bound
    real(4) :: central_diff, ad_result
    logical :: has_large_errors
    complex(4) :: cdotu_forward, cdotu_backward  ! Function result for FD check
    integer :: i, j
    complex(4), dimension(n) :: cy
    complex(4), dimension(n) :: cx

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    ! Forward perturbation: f(x + h)
    cy = cy_orig + h * cy_d_orig
    cx = cx_orig + h * cx_d_orig
    cdotu_forward = cdotu(nsize, cx, 1, cy, 1)

    ! Backward perturbation: f(x - h)
    cy = cy_orig - h * cy_d_orig
    cx = cx_orig - h * cx_d_orig
    cdotu_backward = cdotu(nsize, cx, 1, cy, 1)

    ! Compute central differences and compare with AD results
    central_diff = (cdotu_forward - cdotu_backward) / (2.0e0 * h)
    ad_result = cdotu_d_result
    abs_error = abs(central_diff - ad_result)
    abs_reference = abs(ad_result)
    error_bound = 1.0e-3 + 1.0e-3 * abs_reference
    if (abs_error > error_bound) then
      has_large_errors = .true.
      relative_error = abs_error / max(abs_reference, 1.0e-10)
      write(*,*) 'Large error in function result CDOTU:'
      write(*,*) '  Central diff: ', central_diff
      write(*,*) '  AD result:   ', ad_result
      write(*,*) '  Absolute error:', abs_error
      write(*,*) '  Error bound:', error_bound
      write(*,*) '  Relative error:', relative_error
    end if
    relative_error = abs_error / max(abs_reference, 1.0e-10)
    max_error = max(max_error, relative_error)

    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-3, atol=1.0e-3'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_cdotu