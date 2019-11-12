! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! $Rev: 55 $ $Date: 2014-12-3112:16:59 -0500 (Wed, 31 Dec 2014) $
! FORTRAN 77
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module m_tenstr_disort_errpack
  contains
      SUBROUTINE  ErrMsg( MESSAG, FATAL )

!        Print out a warning or error message;  abort if error

      LOGICAL       FATAL, MsgLim
      CHARACTER*(*) MESSAG
      INTEGER       MaxMsg, NumMsg
      SAVE          MaxMsg, NumMsg, MsgLim
      DATA NumMsg / 0 /,  MaxMsg / 10 /,  MsgLim / .FALSE. /


      IF ( FATAL )  THEN
         WRITE ( *, '(/,2A,/)' )  ' ******* ERROR >>>>>>  ', MESSAG
         STOP
      END IF

      NumMsg = NumMsg + 1
      IF( MsgLim )  RETURN

      IF ( NumMsg.LE.MaxMsg )  THEN
         WRITE ( *, '(/,2A,/)' )  ' ******* WARNING >>>>>>  ', MESSAG
      ELSE
         WRITE ( *,99 )
         MsgLim = .True.
      ENDIF

      RETURN

   99 FORMAT( //,' >>>>>>  TOO MANY WARNING MESSAGES --  ',&
     &   'They will no longer be printed  <<<<<<<', // )
      END

      LOGICAL FUNCTION  WrtBad ( VarNam )

!          Write names of erroneous variables and return 'TRUE'
!
!      INPUT :   VarNam = Name of erroneous variable to be written
!                         ( CHARACTER, any length )

      CHARACTER*(*)  VarNam
      INTEGER        MaxMsg, NumMsg
      SAVE  NumMsg, MaxMsg
      DATA  NumMsg / 0 /,  MaxMsg / 5 /


      WrtBad = .TRUE.
      NumMsg = NumMsg + 1
      WRITE ( *, '(3A)' )  ' ****  Input variable  ', VarNam,&
     &                     '  in error  ****'
      IF ( NumMsg.EQ.MaxMsg )&
     &   CALL  ErrMsg ( 'Too many input errors.  Aborting...', .TRUE. )

      RETURN
      END
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! ---------------------------------------------------------------------
      LOGICAL FUNCTION  WrtDim ( DimNam, MinVal )

!          Write name of too-small symbolic dimension and
!          the value it should be increased to;  return 'TRUE'
!
!      INPUT :  DimNam = Name of symbolic dimension which is too small
!                        ( CHARACTER, any length )
!               Minval = Value to which that dimension should be
!                        increased (at least)

      CHARACTER*(*)  DimNam
      INTEGER        MinVal


      WRITE ( *, '(/,3A,I7)' )  ' ****  Symbolic dimension  ', DimNam,&
     &                     '  should be increased to at least ', MinVal
      WrtDim = .TRUE.

      RETURN
      END
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! ---------------------------------------------------------------------
      LOGICAL FUNCTION  TstBad( VarNam, RelErr )

!       Write name (VarNam) of variable failing self-test and its
!       percent error from the correct value;  return  'FALSE'.

      CHARACTER*(*)  VarNam
      REAL           RelErr


      TstBad = .FALSE.
      WRITE( *, '(/,3A,1P,E11.2,A)' )&
     &       ' Output variable ', VarNam,' differed by ', 100.*RelErr,&
     &       ' per cent from correct value.  Self-test failed.'

      RETURN
      END
end module
