! DIFFSIZES_access.f90 - Module storage for ISIZE parameters (no COMMON)
! Used when many ISIZE vars would exceed F77 line limit in COMMON.
MODULE diffsizes_access
  IMPLICIT NONE
  INTEGER, SAVE :: ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global, ISIZE1OFDx_global, &
    ISIZE1OFDy_global, ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global, &
    ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global, ISIZE2OFB_global
  ! Initialize to invalid so we can detect "not set"
  DATA ISIZE1OFAp_global /-1/, ISIZE1OFCx_global /-1/, ISIZE1OFCy_global /-1/, ISIZE1OFDx_global /-1/, &
    ISIZE1OFDy_global /-1/, ISIZE1OFSx_global /-1/, ISIZE1OFSy_global /-1/, ISIZE1OFX_global /-1/, &
    ISIZE1OFY_global /-1/, ISIZE1OFZx_global /-1/, ISIZE1OFZy_global /-1/, ISIZE2OFA_global /-1/, &
    ISIZE2OFB_global /-1/
CONTAINS

  SUBROUTINE set_ISIZE1OFAp(val)
    INTEGER, INTENT(IN) :: val
    ISIZE1OFAp_global = val
  END SUBROUTINE
  INTEGER FUNCTION get_ISIZE1OFAp()
    get_ISIZE1OFAp = ISIZE1OFAp_global
  END FUNCTION
  SUBROUTINE check_ISIZE1OFAp_initialized()
    IF (ISIZE1OFAp_global < 0) THEN
      WRITE(*,'(A)') 'Error: ISIZE1OFAp_global not set.'
      STOP 1
    END IF
  END SUBROUTINE

  SUBROUTINE set_ISIZE1OFCx(val)
    INTEGER, INTENT(IN) :: val
    ISIZE1OFCx_global = val
  END SUBROUTINE
  INTEGER FUNCTION get_ISIZE1OFCx()
    get_ISIZE1OFCx = ISIZE1OFCx_global
  END FUNCTION
  SUBROUTINE check_ISIZE1OFCx_initialized()
    IF (ISIZE1OFCx_global < 0) THEN
      WRITE(*,'(A)') 'Error: ISIZE1OFCx_global not set.'
      STOP 1
    END IF
  END SUBROUTINE

  SUBROUTINE set_ISIZE1OFCy(val)
    INTEGER, INTENT(IN) :: val
    ISIZE1OFCy_global = val
  END SUBROUTINE
  INTEGER FUNCTION get_ISIZE1OFCy()
    get_ISIZE1OFCy = ISIZE1OFCy_global
  END FUNCTION
  SUBROUTINE check_ISIZE1OFCy_initialized()
    IF (ISIZE1OFCy_global < 0) THEN
      WRITE(*,'(A)') 'Error: ISIZE1OFCy_global not set.'
      STOP 1
    END IF
  END SUBROUTINE

  SUBROUTINE set_ISIZE1OFDx(val)
    INTEGER, INTENT(IN) :: val
    ISIZE1OFDx_global = val
  END SUBROUTINE
  INTEGER FUNCTION get_ISIZE1OFDx()
    get_ISIZE1OFDx = ISIZE1OFDx_global
  END FUNCTION
  SUBROUTINE check_ISIZE1OFDx_initialized()
    IF (ISIZE1OFDx_global < 0) THEN
      WRITE(*,'(A)') 'Error: ISIZE1OFDx_global not set.'
      STOP 1
    END IF
  END SUBROUTINE

  SUBROUTINE set_ISIZE1OFDy(val)
    INTEGER, INTENT(IN) :: val
    ISIZE1OFDy_global = val
  END SUBROUTINE
  INTEGER FUNCTION get_ISIZE1OFDy()
    get_ISIZE1OFDy = ISIZE1OFDy_global
  END FUNCTION
  SUBROUTINE check_ISIZE1OFDy_initialized()
    IF (ISIZE1OFDy_global < 0) THEN
      WRITE(*,'(A)') 'Error: ISIZE1OFDy_global not set.'
      STOP 1
    END IF
  END SUBROUTINE

  SUBROUTINE set_ISIZE1OFSx(val)
    INTEGER, INTENT(IN) :: val
    ISIZE1OFSx_global = val
  END SUBROUTINE
  INTEGER FUNCTION get_ISIZE1OFSx()
    get_ISIZE1OFSx = ISIZE1OFSx_global
  END FUNCTION
  SUBROUTINE check_ISIZE1OFSx_initialized()
    IF (ISIZE1OFSx_global < 0) THEN
      WRITE(*,'(A)') 'Error: ISIZE1OFSx_global not set.'
      STOP 1
    END IF
  END SUBROUTINE

  SUBROUTINE set_ISIZE1OFSy(val)
    INTEGER, INTENT(IN) :: val
    ISIZE1OFSy_global = val
  END SUBROUTINE
  INTEGER FUNCTION get_ISIZE1OFSy()
    get_ISIZE1OFSy = ISIZE1OFSy_global
  END FUNCTION
  SUBROUTINE check_ISIZE1OFSy_initialized()
    IF (ISIZE1OFSy_global < 0) THEN
      WRITE(*,'(A)') 'Error: ISIZE1OFSy_global not set.'
      STOP 1
    END IF
  END SUBROUTINE

  SUBROUTINE set_ISIZE1OFZx(val)
    INTEGER, INTENT(IN) :: val
    ISIZE1OFZx_global = val
  END SUBROUTINE
  INTEGER FUNCTION get_ISIZE1OFZx()
    get_ISIZE1OFZx = ISIZE1OFZx_global
  END FUNCTION
  SUBROUTINE check_ISIZE1OFZx_initialized()
    IF (ISIZE1OFZx_global < 0) THEN
      WRITE(*,'(A)') 'Error: ISIZE1OFZx_global not set.'
      STOP 1
    END IF
  END SUBROUTINE

  SUBROUTINE set_ISIZE1OFZy(val)
    INTEGER, INTENT(IN) :: val
    ISIZE1OFZy_global = val
  END SUBROUTINE
  INTEGER FUNCTION get_ISIZE1OFZy()
    get_ISIZE1OFZy = ISIZE1OFZy_global
  END FUNCTION
  SUBROUTINE check_ISIZE1OFZy_initialized()
    IF (ISIZE1OFZy_global < 0) THEN
      WRITE(*,'(A)') 'Error: ISIZE1OFZy_global not set.'
      STOP 1
    END IF
  END SUBROUTINE


  SUBROUTINE set_ISIZE1OFX(val)
    INTEGER, INTENT(IN) :: val
    ISIZE1OFX_global = val
  END SUBROUTINE

  INTEGER FUNCTION get_ISIZE1OFX()
    get_ISIZE1OFX = ISIZE1OFX_global
  END FUNCTION

  SUBROUTINE check_ISIZE1OFX_initialized()
    IF (ISIZE1OFX_global < 0) THEN
      WRITE(*,'(A)') 'Error: ISIZE1OFX_global not set. Call set_ISIZE1OFX before differentiated routine.'
      STOP 1
    END IF
  END SUBROUTINE

  SUBROUTINE set_ISIZE1OFY(val)
    INTEGER, INTENT(IN) :: val
    ISIZE1OFY_global = val
  END SUBROUTINE

  INTEGER FUNCTION get_ISIZE1OFY()
    get_ISIZE1OFY = ISIZE1OFY_global
  END FUNCTION

  SUBROUTINE check_ISIZE1OFY_initialized()
    IF (ISIZE1OFY_global < 0) THEN
      WRITE(*,'(A)') 'Error: ISIZE1OFY_global not set. Call set_ISIZE1OFY before differentiated routine.'
      STOP 1
    END IF
  END SUBROUTINE

  SUBROUTINE set_ISIZE2OFA(val)
    INTEGER, INTENT(IN) :: val
    ISIZE2OFA_global = val
  END SUBROUTINE

  INTEGER FUNCTION get_ISIZE2OFA()
    get_ISIZE2OFA = ISIZE2OFA_global
  END FUNCTION

  SUBROUTINE check_ISIZE2OFA_initialized()
    IF (ISIZE2OFA_global < 0) THEN
      WRITE(*,'(A)') 'Error: ISIZE2OFA_global not set. Call set_ISIZE2OFA before differentiated routine.'
      STOP 1
    END IF
  END SUBROUTINE

  SUBROUTINE set_ISIZE2OFB(val)
    INTEGER, INTENT(IN) :: val
    ISIZE2OFB_global = val
  END SUBROUTINE

  INTEGER FUNCTION get_ISIZE2OFB()
    get_ISIZE2OFB = ISIZE2OFB_global
  END FUNCTION

  SUBROUTINE check_ISIZE2OFB_initialized()
    IF (ISIZE2OFB_global < 0) THEN
      WRITE(*,'(A)') 'Error: ISIZE2OFB_global not set. Call set_ISIZE2OFB before differentiated routine.'
      STOP 1
    END IF
  END SUBROUTINE

END MODULE diffsizes_access

