MODULE DIFFSIZES
  IMPLICIT NONE
  INTEGER, PARAMETER :: nbdirsmax = 4
  ! ISIZE* are module variables (set via set_ISIZE*(), read via get_ISIZE*() or use directly after check)
  INTEGER, SAVE :: isize1ofx = -1, isize1ofy = -1, isize2ofa = -1
CONTAINS
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

END MODULE DIFFSIZES
