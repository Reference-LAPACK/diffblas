C DIFFSIZES_access.f - Global storage and accessors for ISIZE parameters
C used by differentiated BLAS code. Test code sets these before calling
C the differentiated routine; the routine reads them via getters.
C
      BLOCK DATA diffsizes_init
      INTEGER ISIZE1OFX_global, ISIZE2OFA_global, ISIZE2OFB_global
      COMMON /DIFFSZ/ ISIZE1OFX_global,ISIZE2OFA_global,ISIZE2OFB_global
C     Initialize to invalid value so we can detect "not set"
      DATA ISIZE1OFX_global /-1/
      DATA ISIZE2OFA_global /-1/
      DATA ISIZE2OFB_global /-1/
      END BLOCK DATA

      SUBROUTINE set_ISIZE1OFX(val)
      INTEGER val
      INTEGER ISIZE1OFX_global, ISIZE2OFA_global, ISIZE2OFB_global
      COMMON /DIFFSZ/ ISIZE1OFX_global,ISIZE2OFA_global,ISIZE2OFB_global
      ISIZE1OFX_global = val
      RETURN
      END

      SUBROUTINE set_ISIZE2OFA(val)
      INTEGER val
      INTEGER ISIZE1OFX_global, ISIZE2OFA_global, ISIZE2OFB_global
      COMMON /DIFFSZ/ ISIZE1OFX_global,ISIZE2OFA_global,ISIZE2OFB_global
      ISIZE2OFA_global = val
      RETURN
      END

      SUBROUTINE set_ISIZE2OFB(val)
      INTEGER val
      INTEGER ISIZE1OFX_global, ISIZE2OFA_global, ISIZE2OFB_global
      COMMON /DIFFSZ/ ISIZE1OFX_global,ISIZE2OFA_global,ISIZE2OFB_global
      ISIZE2OFB_global = val
      RETURN
      END

      INTEGER FUNCTION get_ISIZE1OFX()
      INTEGER ISIZE1OFX_global, ISIZE2OFA_global, ISIZE2OFB_global
      COMMON /DIFFSZ/ ISIZE1OFX_global,ISIZE2OFA_global,ISIZE2OFB_global
      get_ISIZE1OFX = ISIZE1OFX_global
      RETURN
      END

      INTEGER FUNCTION get_ISIZE2OFA()
      INTEGER ISIZE1OFX_global, ISIZE2OFA_global, ISIZE2OFB_global
      COMMON /DIFFSZ/ ISIZE1OFX_global,ISIZE2OFA_global,ISIZE2OFB_global
      get_ISIZE2OFA = ISIZE2OFA_global
      RETURN
      END

      INTEGER FUNCTION get_ISIZE2OFB()
      INTEGER ISIZE1OFX_global, ISIZE2OFA_global, ISIZE2OFB_global
      COMMON /DIFFSZ/ ISIZE1OFX_global,ISIZE2OFA_global,ISIZE2OFB_global
      get_ISIZE2OFB = ISIZE2OFB_global
      RETURN
      END

C     Check that ISIZE1OFX_global has been set; stop with message if not.
      SUBROUTINE check_ISIZE1OFX_initialized()
      INTEGER ISIZE1OFX_global, ISIZE2OFA_global, ISIZE2OFB_global
      COMMON /DIFFSZ/ ISIZE1OFX_global,ISIZE2OFA_global,ISIZE2OFB_global
      IF (ISIZE1OFX_global .LT. 0) THEN
        WRITE(*,'(A)') 'Error: ISIZE1OFX_global not set. Call set_ISIZE'
     & // '1OFX before differentiated routine.'
        STOP 1
      END IF
      RETURN
      END

C     Check that ISIZE2OFA_global has been set; stop with message if not.
      SUBROUTINE check_ISIZE2OFA_initialized()
      INTEGER ISIZE1OFX_global, ISIZE2OFA_global, ISIZE2OFB_global
      COMMON /DIFFSZ/ ISIZE1OFX_global,ISIZE2OFA_global,ISIZE2OFB_global
      IF (ISIZE2OFA_global .LT. 0) THEN
        WRITE(*,'(A)') 'Error: ISIZE2OFA_global not set. Call set_ISIZE'
     & // '2OFA before differentiated routine.'
        STOP 1
      END IF
      RETURN
      END

C     Check that ISIZE2OFB_global has been set; stop with message if not.
      SUBROUTINE check_ISIZE2OFB_initialized()
      INTEGER ISIZE1OFX_global, ISIZE2OFA_global, ISIZE2OFB_global
      COMMON /DIFFSZ/ ISIZE1OFX_global,ISIZE2OFA_global,ISIZE2OFB_global
      IF (ISIZE2OFB_global .LT. 0) THEN
        WRITE(*,'(A)') 'Error: ISIZE2OFB_global not set. Call set_ISIZE'
     & // '2OFB before differentiated routine.'
        STOP 1
      END IF
      RETURN
      END

