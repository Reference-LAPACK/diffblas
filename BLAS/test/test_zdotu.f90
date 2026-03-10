! Test program for ZDOTU differentiation
! Generated automatically by run_tapenade_blas.py
! Using REAL*8 precision
! Multi-size test with outlined run_test_for_size(n) - arrays declared to size n

program test_zdotu
  implicit none

  complex(8), external :: zdotu
  complex(8), external :: zdotu_d

  integer :: n_test
  integer :: seed_array(33)
  integer :: test_sizes(1)
  integer :: i
  logical :: passed, all_passed

  seed_array = 42
  call random_seed(put=seed_array)

  test_sizes = (/ 4 /)
  write(*,*) 'Testing ZDOTU (multi-size: n = 4)'
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
    complex(8), dimension(n) :: zx
    integer :: incx
    complex(8), dimension(n) :: zy
    integer :: incy

    ! Derivative variables
    complex(8) :: zdotu_d_result  ! Derivative of function result (avoid name clash with func_d)
    complex(8), dimension(n) :: zx_d
    complex(8), dimension(n) :: zy_d

    ! Array restoration and derivative storage
    complex(8) :: zdotu_orig  ! Function result (no _d_orig - use _d_result)
    complex(8), dimension(n) :: zx_orig, zx_d_orig
    complex(8), dimension(n) :: zy_orig, zy_d_orig
    real(8) :: temp_re, temp_im  ! For complex random init
    integer :: i, j

    nsize = n
    incx = 1
    incy = 1

    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      zx(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    end do
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      zy(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    end do

    ! Initialize input derivatives
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      zx_d(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    end do
    do i = 1, n
      call random_number(temp_re)
      call random_number(temp_im)
      zy_d(i) = cmplx(temp_re * 2.0 - 1.0, temp_im * 2.0 - 1.0, kind=8)
    end do

    ! Store _orig and _d_orig
    zx_d_orig = zx_d
    zy_d_orig = zy_d
    zdotu_orig = zdotu(nsize, zx, 1, zy, 1)
    zx_orig = zx
    zy_orig = zy

    write(*,*) 'Testing ZDOTU (n =', n, ')'

    ! Call the differentiated function
    zdotu_d_result = zdotu_d(nsize, zx, zx_d, 1, zy, zy_d, 1, zdotu_orig)

    write(*,*) 'Function calls completed successfully'

    ! Numerical differentiation check
    call check_derivatives_numerically(n, nsize, zx_orig, zy_orig, zdotu_orig, zx_d_orig, zy_d_orig, zdotu_d_result, passed)

  end subroutine run_test_for_size

  subroutine check_derivatives_numerically(n, nsize, zx_orig, zy_orig, zdotu_orig, zx_d_orig, zy_d_orig, zdotu_d_result, passed)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nsize
    complex(8), intent(in) :: zx_orig(n), zx_d_orig(n)
    complex(8), intent(in) :: zy_orig(n), zy_d_orig(n)
    complex(8), intent(in) :: zdotu_orig
    complex(8), intent(in) :: zdotu_d_result
    logical, intent(out) :: passed

    real(8), parameter :: h = 1.0e-6  ! Step size for finite differences
    real(8) :: relative_error, max_error
    real(8) :: abs_error, abs_reference, error_bound
    real(8) :: central_diff, ad_result
    logical :: has_large_errors
    complex(8) :: zdotu_forward, zdotu_backward  ! Function result for FD check
    integer :: i, j
    complex(8), dimension(n) :: zx
    complex(8), dimension(n) :: zy

    max_error = 0.0e0
    has_large_errors = .false.

    write(*,*) 'Checking derivatives against numerical differentiation:'
    write(*,*) 'Step size h =', h

    ! Forward perturbation: f(x + h)
    zx = zx_orig + h * zx_d_orig
    zy = zy_orig + h * zy_d_orig
    zdotu_forward = zdotu(nsize, zx, 1, zy, 1)

    ! Backward perturbation: f(x - h)
    zx = zx_orig - h * zx_d_orig
    zy = zy_orig - h * zy_d_orig
    zdotu_backward = zdotu(nsize, zx, 1, zy, 1)

    ! Compute central differences and compare with AD results
    central_diff = (zdotu_forward - zdotu_backward) / (2.0e0 * h)
    ad_result = zdotu_d_result
    abs_error = abs(central_diff - ad_result)
    abs_reference = abs(ad_result)
    error_bound = 1.0e-5 + 1.0e-5 * abs_reference
    if (abs_error > error_bound) then
      has_large_errors = .true.
      relative_error = abs_error / max(abs_reference, 1.0e-10)
      write(*,*) 'Large error in function result ZDOTU:'
      write(*,*) '  Central diff: ', central_diff
      write(*,*) '  AD result:   ', ad_result
      write(*,*) '  Absolute error:', abs_error
      write(*,*) '  Error bound:', error_bound
      write(*,*) '  Relative error:', relative_error
    end if
    relative_error = abs_error / max(abs_reference, 1.0e-10)
    max_error = max(max_error, relative_error)

    write(*,*) 'Maximum relative error:', max_error
    write(*,*) 'Tolerance thresholds: rtol=1.0e-5, atol=1.0e-5'
    passed = .not. has_large_errors
    if (has_large_errors) then
      write(*,*) 'FAIL: Large errors detected in derivatives (outside tolerance)'
    else
      write(*,*) 'PASS: Derivatives are within tolerance (rtol + atol)'
    end if

  end subroutine check_derivatives_numerically

end program test_zdotu