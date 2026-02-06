MODULE DIFFSIZES
  IMPLICIT NONE
  INTEGER, PARAMETER :: nbdirsmax = 4
  ! ISIZE* are module variables (set via set_ISIZE*(), read via get_ISIZE*() or use directly after check)
  INTEGER, SAVE :: isize1ofap = -1, isize1ofcx = -1, isize1ofcy = -1, isize1ofdx = -1, isize1ofdy = -1, isize1ofsx = -1, &
    & isize1ofsy = -1, isize1ofx = -1, isize1ofy = -1, isize2ofa = -1, isize2ofb = -1
CONTAINS
  SUBROUTINE set_ISIZE1OFAp(val)
    INTEGER, INTENT(IN) :: val
    isize1ofap = val
  END SUBROUTINE set_ISIZE1OFAp

  INTEGER FUNCTION get_ISIZE1OFAp()
    get_ISIZE1OFAp = isize1ofap
  END FUNCTION get_ISIZE1OFAp

  SUBROUTINE check_ISIZE1OFAp_initialized()
    IF (isize1ofap < 0) THEN
      WRITE(*,'(A)') 'Error: isize1ofap not set. Call set_ISIZE1OFAp before differentiated routine.'
      STOP 1
    END IF
  END SUBROUTINE check_ISIZE1OFAp_initialized

  SUBROUTINE set_ISIZE1OFCx(val)
    INTEGER, INTENT(IN) :: val
    isize1ofcx = val
  END SUBROUTINE set_ISIZE1OFCx

  INTEGER FUNCTION get_ISIZE1OFCx()
    get_ISIZE1OFCx = isize1ofcx
  END FUNCTION get_ISIZE1OFCx

  SUBROUTINE check_ISIZE1OFCx_initialized()
    IF (isize1ofcx < 0) THEN
      WRITE(*,'(A)') 'Error: isize1ofcx not set. Call set_ISIZE1OFCx before differentiated routine.'
      STOP 1
    END IF
  END SUBROUTINE check_ISIZE1OFCx_initialized

  SUBROUTINE set_ISIZE1OFCy(val)
    INTEGER, INTENT(IN) :: val
    isize1ofcy = val
  END SUBROUTINE set_ISIZE1OFCy

  INTEGER FUNCTION get_ISIZE1OFCy()
    get_ISIZE1OFCy = isize1ofcy
  END FUNCTION get_ISIZE1OFCy

  SUBROUTINE check_ISIZE1OFCy_initialized()
    IF (isize1ofcy < 0) THEN
      WRITE(*,'(A)') 'Error: isize1ofcy not set. Call set_ISIZE1OFCy before differentiated routine.'
      STOP 1
    END IF
  END SUBROUTINE check_ISIZE1OFCy_initialized

  SUBROUTINE set_ISIZE1OFDx(val)
    INTEGER, INTENT(IN) :: val
    isize1ofdx = val
  END SUBROUTINE set_ISIZE1OFDx

  INTEGER FUNCTION get_ISIZE1OFDx()
    get_ISIZE1OFDx = isize1ofdx
  END FUNCTION get_ISIZE1OFDx

  SUBROUTINE check_ISIZE1OFDx_initialized()
    IF (isize1ofdx < 0) THEN
      WRITE(*,'(A)') 'Error: isize1ofdx not set. Call set_ISIZE1OFDx before differentiated routine.'
      STOP 1
    END IF
  END SUBROUTINE check_ISIZE1OFDx_initialized

  SUBROUTINE set_ISIZE1OFDy(val)
    INTEGER, INTENT(IN) :: val
    isize1ofdy = val
  END SUBROUTINE set_ISIZE1OFDy

  INTEGER FUNCTION get_ISIZE1OFDy()
    get_ISIZE1OFDy = isize1ofdy
  END FUNCTION get_ISIZE1OFDy

  SUBROUTINE check_ISIZE1OFDy_initialized()
    IF (isize1ofdy < 0) THEN
      WRITE(*,'(A)') 'Error: isize1ofdy not set. Call set_ISIZE1OFDy before differentiated routine.'
      STOP 1
    END IF
  END SUBROUTINE check_ISIZE1OFDy_initialized

  SUBROUTINE set_ISIZE1OFSx(val)
    INTEGER, INTENT(IN) :: val
    isize1ofsx = val
  END SUBROUTINE set_ISIZE1OFSx

  INTEGER FUNCTION get_ISIZE1OFSx()
    get_ISIZE1OFSx = isize1ofsx
  END FUNCTION get_ISIZE1OFSx

  SUBROUTINE check_ISIZE1OFSx_initialized()
    IF (isize1ofsx < 0) THEN
      WRITE(*,'(A)') 'Error: isize1ofsx not set. Call set_ISIZE1OFSx before differentiated routine.'
      STOP 1
    END IF
  END SUBROUTINE check_ISIZE1OFSx_initialized

  SUBROUTINE set_ISIZE1OFSy(val)
    INTEGER, INTENT(IN) :: val
    isize1ofsy = val
  END SUBROUTINE set_ISIZE1OFSy

  INTEGER FUNCTION get_ISIZE1OFSy()
    get_ISIZE1OFSy = isize1ofsy
  END FUNCTION get_ISIZE1OFSy

  SUBROUTINE check_ISIZE1OFSy_initialized()
    IF (isize1ofsy < 0) THEN
      WRITE(*,'(A)') 'Error: isize1ofsy not set. Call set_ISIZE1OFSy before differentiated routine.'
      STOP 1
    END IF
  END SUBROUTINE check_ISIZE1OFSy_initialized

  SUBROUTINE set_ISIZE1OFX(val)
    INTEGER, INTENT(IN) :: val
    isize1ofx = val
  END SUBROUTINE set_ISIZE1OFX

  INTEGER FUNCTION get_ISIZE1OFX()
    get_ISIZE1OFX = isize1ofx
  END FUNCTION get_ISIZE1OFX

  SUBROUTINE check_ISIZE1OFX_initialized()
    IF (isize1ofx < 0) THEN
      WRITE(*,'(A)') 'Error: isize1ofx not set. Call set_ISIZE1OFX before differentiated routine.'
      STOP 1
    END IF
  END SUBROUTINE check_ISIZE1OFX_initialized

  SUBROUTINE set_ISIZE1OFY(val)
    INTEGER, INTENT(IN) :: val
    isize1ofy = val
  END SUBROUTINE set_ISIZE1OFY

  INTEGER FUNCTION get_ISIZE1OFY()
    get_ISIZE1OFY = isize1ofy
  END FUNCTION get_ISIZE1OFY

  SUBROUTINE check_ISIZE1OFY_initialized()
    IF (isize1ofy < 0) THEN
      WRITE(*,'(A)') 'Error: isize1ofy not set. Call set_ISIZE1OFY before differentiated routine.'
      STOP 1
    END IF
  END SUBROUTINE check_ISIZE1OFY_initialized

  SUBROUTINE set_ISIZE2OFA(val)
    INTEGER, INTENT(IN) :: val
    isize2ofa = val
  END SUBROUTINE set_ISIZE2OFA

  INTEGER FUNCTION get_ISIZE2OFA()
    get_ISIZE2OFA = isize2ofa
  END FUNCTION get_ISIZE2OFA

  SUBROUTINE check_ISIZE2OFA_initialized()
    IF (isize2ofa < 0) THEN
      WRITE(*,'(A)') 'Error: isize2ofa not set. Call set_ISIZE2OFA before differentiated routine.'
      STOP 1
    END IF
  END SUBROUTINE check_ISIZE2OFA_initialized

  SUBROUTINE set_ISIZE2OFB(val)
    INTEGER, INTENT(IN) :: val
    isize2ofb = val
  END SUBROUTINE set_ISIZE2OFB

  INTEGER FUNCTION get_ISIZE2OFB()
    get_ISIZE2OFB = isize2ofb
  END FUNCTION get_ISIZE2OFB

  SUBROUTINE check_ISIZE2OFB_initialized()
    IF (isize2ofb < 0) THEN
      WRITE(*,'(A)') 'Error: isize2ofb not set. Call set_ISIZE2OFB before differentiated routine.'
      STOP 1
    END IF
  END SUBROUTINE check_ISIZE2OFB_initialized

END MODULE DIFFSIZES
