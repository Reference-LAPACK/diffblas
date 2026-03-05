C DIFFSIZES_access_wrappers.f - External interface for DIFFSIZES_access module
C C and .f callers expect set_isize*_, get_isize*_, etc.; the F90 module exports
C __diffsizes_access_MOD_* names. These wrappers provide the expected external symbols.
C
      SUBROUTINE set_ISIZE1OFX(val)
      USE diffsizes_access, ONLY: ISIZE1OFX_global
      INTEGER val
      ISIZE1OFX_global = val
      RETURN
      END

      INTEGER FUNCTION get_ISIZE1OFX()
      USE diffsizes_access, ONLY: ISIZE1OFX_global
      get_ISIZE1OFX = ISIZE1OFX_global
      RETURN
      END

      SUBROUTINE check_ISIZE1OFX_initialized()
      USE diffsizes_access, ONLY: ISIZE1OFX_global
      IF (ISIZE1OFX_global .LT. 0) THEN
        WRITE(6,*) 'Error: ISIZE not set before differentiated routine'
        STOP 1
      END IF
      RETURN
      END

      SUBROUTINE set_ISIZE1OFY(val)
      USE diffsizes_access, ONLY: ISIZE1OFY_global
      INTEGER val
      ISIZE1OFY_global = val
      RETURN
      END

      INTEGER FUNCTION get_ISIZE1OFY()
      USE diffsizes_access, ONLY: ISIZE1OFY_global
      get_ISIZE1OFY = ISIZE1OFY_global
      RETURN
      END

      SUBROUTINE check_ISIZE1OFY_initialized()
      USE diffsizes_access, ONLY: ISIZE1OFY_global
      IF (ISIZE1OFY_global .LT. 0) THEN
        WRITE(6,*) 'Error: ISIZE not set before differentiated routine'
        STOP 1
      END IF
      RETURN
      END

      SUBROUTINE set_ISIZE2OFA(val)
      USE diffsizes_access, ONLY: ISIZE2OFA_global
      INTEGER val
      ISIZE2OFA_global = val
      RETURN
      END

      INTEGER FUNCTION get_ISIZE2OFA()
      USE diffsizes_access, ONLY: ISIZE2OFA_global
      get_ISIZE2OFA = ISIZE2OFA_global
      RETURN
      END

      SUBROUTINE check_ISIZE2OFA_initialized()
      USE diffsizes_access, ONLY: ISIZE2OFA_global
      IF (ISIZE2OFA_global .LT. 0) THEN
        WRITE(6,*) 'Error: ISIZE not set before differentiated routine'
        STOP 1
      END IF
      RETURN
      END

      SUBROUTINE set_ISIZE2OFB(val)
      USE diffsizes_access, ONLY: ISIZE2OFB_global
      INTEGER val
      ISIZE2OFB_global = val
      RETURN
      END

      INTEGER FUNCTION get_ISIZE2OFB()
      USE diffsizes_access, ONLY: ISIZE2OFB_global
      get_ISIZE2OFB = ISIZE2OFB_global
      RETURN
      END

      SUBROUTINE check_ISIZE2OFB_initialized()
      USE diffsizes_access, ONLY: ISIZE2OFB_global
      IF (ISIZE2OFB_global .LT. 0) THEN
        WRITE(6,*) 'Error: ISIZE not set before differentiated routine'
        STOP 1
      END IF
      RETURN
      END

      SUBROUTINE set_ISIZE1OFAp(val)
      USE diffsizes_access, ONLY: ISIZE1OFAp_global
      INTEGER val
      ISIZE1OFAp_global = val
      RETURN
      END
      INTEGER FUNCTION get_ISIZE1OFAp()
      USE diffsizes_access, ONLY: ISIZE1OFAp_global
      get_ISIZE1OFAp = ISIZE1OFAp_global
      RETURN
      END
      SUBROUTINE check_ISIZE1OFAp_initialized()
      USE diffsizes_access, ONLY: ISIZE1OFAp_global
      IF (ISIZE1OFAp_global .LT. 0) THEN
        WRITE(6,*) 'Error: ISIZE not set before differentiated routine'
        STOP 1
      END IF
      RETURN
      END

      SUBROUTINE set_ISIZE1OFCx(val)
      USE diffsizes_access, ONLY: ISIZE1OFCx_global
      INTEGER val
      ISIZE1OFCx_global = val
      RETURN
      END
      INTEGER FUNCTION get_ISIZE1OFCx()
      USE diffsizes_access, ONLY: ISIZE1OFCx_global
      get_ISIZE1OFCx = ISIZE1OFCx_global
      RETURN
      END
      SUBROUTINE check_ISIZE1OFCx_initialized()
      USE diffsizes_access, ONLY: ISIZE1OFCx_global
      IF (ISIZE1OFCx_global .LT. 0) THEN
        WRITE(6,*) 'Error: ISIZE not set before differentiated routine'
        STOP 1
      END IF
      RETURN
      END

      SUBROUTINE set_ISIZE1OFCy(val)
      USE diffsizes_access, ONLY: ISIZE1OFCy_global
      INTEGER val
      ISIZE1OFCy_global = val
      RETURN
      END
      INTEGER FUNCTION get_ISIZE1OFCy()
      USE diffsizes_access, ONLY: ISIZE1OFCy_global
      get_ISIZE1OFCy = ISIZE1OFCy_global
      RETURN
      END
      SUBROUTINE check_ISIZE1OFCy_initialized()
      USE diffsizes_access, ONLY: ISIZE1OFCy_global
      IF (ISIZE1OFCy_global .LT. 0) THEN
        WRITE(6,*) 'Error: ISIZE not set before differentiated routine'
        STOP 1
      END IF
      RETURN
      END

      SUBROUTINE set_ISIZE1OFDx(val)
      USE diffsizes_access, ONLY: ISIZE1OFDx_global
      INTEGER val
      ISIZE1OFDx_global = val
      RETURN
      END
      INTEGER FUNCTION get_ISIZE1OFDx()
      USE diffsizes_access, ONLY: ISIZE1OFDx_global
      get_ISIZE1OFDx = ISIZE1OFDx_global
      RETURN
      END
      SUBROUTINE check_ISIZE1OFDx_initialized()
      USE diffsizes_access, ONLY: ISIZE1OFDx_global
      IF (ISIZE1OFDx_global .LT. 0) THEN
        WRITE(6,*) 'Error: ISIZE not set before differentiated routine'
        STOP 1
      END IF
      RETURN
      END

      SUBROUTINE set_ISIZE1OFDy(val)
      USE diffsizes_access, ONLY: ISIZE1OFDy_global
      INTEGER val
      ISIZE1OFDy_global = val
      RETURN
      END
      INTEGER FUNCTION get_ISIZE1OFDy()
      USE diffsizes_access, ONLY: ISIZE1OFDy_global
      get_ISIZE1OFDy = ISIZE1OFDy_global
      RETURN
      END
      SUBROUTINE check_ISIZE1OFDy_initialized()
      USE diffsizes_access, ONLY: ISIZE1OFDy_global
      IF (ISIZE1OFDy_global .LT. 0) THEN
        WRITE(6,*) 'Error: ISIZE not set before differentiated routine'
        STOP 1
      END IF
      RETURN
      END

      SUBROUTINE set_ISIZE1OFSx(val)
      USE diffsizes_access, ONLY: ISIZE1OFSx_global
      INTEGER val
      ISIZE1OFSx_global = val
      RETURN
      END
      INTEGER FUNCTION get_ISIZE1OFSx()
      USE diffsizes_access, ONLY: ISIZE1OFSx_global
      get_ISIZE1OFSx = ISIZE1OFSx_global
      RETURN
      END
      SUBROUTINE check_ISIZE1OFSx_initialized()
      USE diffsizes_access, ONLY: ISIZE1OFSx_global
      IF (ISIZE1OFSx_global .LT. 0) THEN
        WRITE(6,*) 'Error: ISIZE not set before differentiated routine'
        STOP 1
      END IF
      RETURN
      END

      SUBROUTINE set_ISIZE1OFSy(val)
      USE diffsizes_access, ONLY: ISIZE1OFSy_global
      INTEGER val
      ISIZE1OFSy_global = val
      RETURN
      END
      INTEGER FUNCTION get_ISIZE1OFSy()
      USE diffsizes_access, ONLY: ISIZE1OFSy_global
      get_ISIZE1OFSy = ISIZE1OFSy_global
      RETURN
      END
      SUBROUTINE check_ISIZE1OFSy_initialized()
      USE diffsizes_access, ONLY: ISIZE1OFSy_global
      IF (ISIZE1OFSy_global .LT. 0) THEN
        WRITE(6,*) 'Error: ISIZE not set before differentiated routine'
        STOP 1
      END IF
      RETURN
      END

      SUBROUTINE set_ISIZE1OFZx(val)
      USE diffsizes_access, ONLY: ISIZE1OFZx_global
      INTEGER val
      ISIZE1OFZx_global = val
      RETURN
      END
      INTEGER FUNCTION get_ISIZE1OFZx()
      USE diffsizes_access, ONLY: ISIZE1OFZx_global
      get_ISIZE1OFZx = ISIZE1OFZx_global
      RETURN
      END
      SUBROUTINE check_ISIZE1OFZx_initialized()
      USE diffsizes_access, ONLY: ISIZE1OFZx_global
      IF (ISIZE1OFZx_global .LT. 0) THEN
        WRITE(6,*) 'Error: ISIZE not set before differentiated routine'
        STOP 1
      END IF
      RETURN
      END

      SUBROUTINE set_ISIZE1OFZy(val)
      USE diffsizes_access, ONLY: ISIZE1OFZy_global
      INTEGER val
      ISIZE1OFZy_global = val
      RETURN
      END
      INTEGER FUNCTION get_ISIZE1OFZy()
      USE diffsizes_access, ONLY: ISIZE1OFZy_global
      get_ISIZE1OFZy = ISIZE1OFZy_global
      RETURN
      END
      SUBROUTINE check_ISIZE1OFZy_initialized()
      USE diffsizes_access, ONLY: ISIZE1OFZy_global
      IF (ISIZE1OFZy_global .LT. 0) THEN
        WRITE(6,*) 'Error: ISIZE not set before differentiated routine'
        STOP 1
      END IF
      RETURN
      END

