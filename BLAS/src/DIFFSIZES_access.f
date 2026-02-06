C DIFFSIZES_access.f - Global storage and accessors for ISIZE parameters
C used by differentiated BLAS code. Test code sets these before calling
C the differentiated routine; the routine reads them via getters.
C
      BLOCK DATA diffsizes_init
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFCy_global, ISIZE1OFDx_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFDy_global, ISIZE1OFSx_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFSy_global, ISIZE1OFX_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFY_global, ISIZE1OFZx_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFZy_global, ISIZE2OFA_global
      COMMON /DIFFSIZES_COMMON/ ISIZE2OFB_global
C     Initialize to invalid value so we can detect "not set"
      DATA ISIZE1OFAp_global /-1/
      DATA ISIZE1OFCx_global /-1/
      DATA ISIZE1OFCy_global /-1/
      DATA ISIZE1OFDx_global /-1/
      DATA ISIZE1OFDy_global /-1/
      DATA ISIZE1OFSx_global /-1/
      DATA ISIZE1OFSy_global /-1/
      DATA ISIZE1OFX_global /-1/
      DATA ISIZE1OFY_global /-1/
      DATA ISIZE1OFZx_global /-1/
      DATA ISIZE1OFZy_global /-1/
      DATA ISIZE2OFA_global /-1/
      DATA ISIZE2OFB_global /-1/
      END BLOCK DATA

      SUBROUTINE set_ISIZE1OFAp(val)
      INTEGER val
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      ISIZE1OFAp_global = val
      RETURN
      END

      SUBROUTINE set_ISIZE1OFCx(val)
      INTEGER val
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      ISIZE1OFCx_global = val
      RETURN
      END

      SUBROUTINE set_ISIZE1OFCy(val)
      INTEGER val
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      ISIZE1OFCy_global = val
      RETURN
      END

      SUBROUTINE set_ISIZE1OFDx(val)
      INTEGER val
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      ISIZE1OFDx_global = val
      RETURN
      END

      SUBROUTINE set_ISIZE1OFDy(val)
      INTEGER val
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      ISIZE1OFDy_global = val
      RETURN
      END

      SUBROUTINE set_ISIZE1OFSx(val)
      INTEGER val
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      ISIZE1OFSx_global = val
      RETURN
      END

      SUBROUTINE set_ISIZE1OFSy(val)
      INTEGER val
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      ISIZE1OFSy_global = val
      RETURN
      END

      SUBROUTINE set_ISIZE1OFX(val)
      INTEGER val
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      ISIZE1OFX_global = val
      RETURN
      END

      SUBROUTINE set_ISIZE1OFY(val)
      INTEGER val
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      ISIZE1OFY_global = val
      RETURN
      END

      SUBROUTINE set_ISIZE1OFZx(val)
      INTEGER val
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      ISIZE1OFZx_global = val
      RETURN
      END

      SUBROUTINE set_ISIZE1OFZy(val)
      INTEGER val
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      ISIZE1OFZy_global = val
      RETURN
      END

      SUBROUTINE set_ISIZE2OFA(val)
      INTEGER val
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      ISIZE2OFA_global = val
      RETURN
      END

      SUBROUTINE set_ISIZE2OFB(val)
      INTEGER val
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      ISIZE2OFB_global = val
      RETURN
      END

      INTEGER FUNCTION get_ISIZE1OFAp()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      get_ISIZE1OFAp = ISIZE1OFAp_global
      RETURN
      END

      INTEGER FUNCTION get_ISIZE1OFCx()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      get_ISIZE1OFCx = ISIZE1OFCx_global
      RETURN
      END

      INTEGER FUNCTION get_ISIZE1OFCy()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      get_ISIZE1OFCy = ISIZE1OFCy_global
      RETURN
      END

      INTEGER FUNCTION get_ISIZE1OFDx()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      get_ISIZE1OFDx = ISIZE1OFDx_global
      RETURN
      END

      INTEGER FUNCTION get_ISIZE1OFDy()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      get_ISIZE1OFDy = ISIZE1OFDy_global
      RETURN
      END

      INTEGER FUNCTION get_ISIZE1OFSx()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      get_ISIZE1OFSx = ISIZE1OFSx_global
      RETURN
      END

      INTEGER FUNCTION get_ISIZE1OFSy()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      get_ISIZE1OFSy = ISIZE1OFSy_global
      RETURN
      END

      INTEGER FUNCTION get_ISIZE1OFX()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      get_ISIZE1OFX = ISIZE1OFX_global
      RETURN
      END

      INTEGER FUNCTION get_ISIZE1OFY()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      get_ISIZE1OFY = ISIZE1OFY_global
      RETURN
      END

      INTEGER FUNCTION get_ISIZE1OFZx()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      get_ISIZE1OFZx = ISIZE1OFZx_global
      RETURN
      END

      INTEGER FUNCTION get_ISIZE1OFZy()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      get_ISIZE1OFZy = ISIZE1OFZy_global
      RETURN
      END

      INTEGER FUNCTION get_ISIZE2OFA()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      get_ISIZE2OFA = ISIZE2OFA_global
      RETURN
      END

      INTEGER FUNCTION get_ISIZE2OFB()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      get_ISIZE2OFB = ISIZE2OFB_global
      RETURN
      END

C     Check that ISIZE1OFAp_global has been set; stop with message if not.
      SUBROUTINE check_ISIZE1OFAp_initialized()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      IF (ISIZE1OFAp_global .LT. 0) THEN
        WRITE(*,'(A)') 'Error: ISIZE1OFAp_global not set. Call set_ISIZ'
     & // 'E1OFAp before differentiated routine.'
        STOP 1
      END IF
      RETURN
      END

C     Check that ISIZE1OFCx_global has been set; stop with message if not.
      SUBROUTINE check_ISIZE1OFCx_initialized()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      IF (ISIZE1OFCx_global .LT. 0) THEN
        WRITE(*,'(A)') 'Error: ISIZE1OFCx_global not set. Call set_ISIZ'
     & // 'E1OFCx before differentiated routine.'
        STOP 1
      END IF
      RETURN
      END

C     Check that ISIZE1OFCy_global has been set; stop with message if not.
      SUBROUTINE check_ISIZE1OFCy_initialized()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      IF (ISIZE1OFCy_global .LT. 0) THEN
        WRITE(*,'(A)') 'Error: ISIZE1OFCy_global not set. Call set_ISIZ'
     & // 'E1OFCy before differentiated routine.'
        STOP 1
      END IF
      RETURN
      END

C     Check that ISIZE1OFDx_global has been set; stop with message if not.
      SUBROUTINE check_ISIZE1OFDx_initialized()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      IF (ISIZE1OFDx_global .LT. 0) THEN
        WRITE(*,'(A)') 'Error: ISIZE1OFDx_global not set. Call set_ISIZ'
     & // 'E1OFDx before differentiated routine.'
        STOP 1
      END IF
      RETURN
      END

C     Check that ISIZE1OFDy_global has been set; stop with message if not.
      SUBROUTINE check_ISIZE1OFDy_initialized()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      IF (ISIZE1OFDy_global .LT. 0) THEN
        WRITE(*,'(A)') 'Error: ISIZE1OFDy_global not set. Call set_ISIZ'
     & // 'E1OFDy before differentiated routine.'
        STOP 1
      END IF
      RETURN
      END

C     Check that ISIZE1OFSx_global has been set; stop with message if not.
      SUBROUTINE check_ISIZE1OFSx_initialized()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      IF (ISIZE1OFSx_global .LT. 0) THEN
        WRITE(*,'(A)') 'Error: ISIZE1OFSx_global not set. Call set_ISIZ'
     & // 'E1OFSx before differentiated routine.'
        STOP 1
      END IF
      RETURN
      END

C     Check that ISIZE1OFSy_global has been set; stop with message if not.
      SUBROUTINE check_ISIZE1OFSy_initialized()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      IF (ISIZE1OFSy_global .LT. 0) THEN
        WRITE(*,'(A)') 'Error: ISIZE1OFSy_global not set. Call set_ISIZ'
     & // 'E1OFSy before differentiated routine.'
        STOP 1
      END IF
      RETURN
      END

C     Check that ISIZE1OFX_global has been set; stop with message if not.
      SUBROUTINE check_ISIZE1OFX_initialized()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      IF (ISIZE1OFX_global .LT. 0) THEN
        WRITE(*,'(A)') 'Error: ISIZE1OFX_global not set. Call set_ISIZE'
     & // '1OFX before differentiated routine.'
        STOP 1
      END IF
      RETURN
      END

C     Check that ISIZE1OFY_global has been set; stop with message if not.
      SUBROUTINE check_ISIZE1OFY_initialized()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      IF (ISIZE1OFY_global .LT. 0) THEN
        WRITE(*,'(A)') 'Error: ISIZE1OFY_global not set. Call set_ISIZE'
     & // '1OFY before differentiated routine.'
        STOP 1
      END IF
      RETURN
      END

C     Check that ISIZE1OFZx_global has been set; stop with message if not.
      SUBROUTINE check_ISIZE1OFZx_initialized()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      IF (ISIZE1OFZx_global .LT. 0) THEN
        WRITE(*,'(A)') 'Error: ISIZE1OFZx_global not set. Call set_ISIZ'
     & // 'E1OFZx before differentiated routine.'
        STOP 1
      END IF
      RETURN
      END

C     Check that ISIZE1OFZy_global has been set; stop with message if not.
      SUBROUTINE check_ISIZE1OFZy_initialized()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      IF (ISIZE1OFZy_global .LT. 0) THEN
        WRITE(*,'(A)') 'Error: ISIZE1OFZy_global not set. Call set_ISIZ'
     & // 'E1OFZy before differentiated routine.'
        STOP 1
      END IF
      RETURN
      END

C     Check that ISIZE2OFA_global has been set; stop with message if not.
      SUBROUTINE check_ISIZE2OFA_initialized()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      IF (ISIZE2OFA_global .LT. 0) THEN
        WRITE(*,'(A)') 'Error: ISIZE2OFA_global not set. Call set_ISIZE'
     & // '2OFA before differentiated routine.'
        STOP 1
      END IF
      RETURN
      END

C     Check that ISIZE2OFB_global has been set; stop with message if not.
      SUBROUTINE check_ISIZE2OFB_initialized()
      INTEGER ISIZE1OFAp_global, ISIZE1OFCx_global, ISIZE1OFCy_global
     & ISIZE1OFDx_global, ISIZE1OFDy_global, ISIZE1OFSx_global
     & ISIZE1OFSy_global, ISIZE1OFX_global, ISIZE1OFY_global
     & ISIZE1OFZx_global, ISIZE1OFZy_global, ISIZE2OFA_global
     & ISIZE2OFB_global
      COMMON /DIFFSIZES_COMMON/ ISIZE1OFAp_global, ISIZE1OFCx_global
     & ISIZE1OFCy_global, ISIZE1OFDx_global, ISIZE1OFDy_global
     & ISIZE1OFSx_global, ISIZE1OFSy_global, ISIZE1OFX_global
     & ISIZE1OFY_global, ISIZE1OFZx_global, ISIZE1OFZy_global
     & ISIZE2OFA_global, ISIZE2OFB_global
      IF (ISIZE2OFB_global .LT. 0) THEN
        WRITE(*,'(A)') 'Error: ISIZE2OFB_global not set. Call set_ISIZE'
     & // '2OFB before differentiated routine.'
        STOP 1
      END IF
      RETURN
      END

