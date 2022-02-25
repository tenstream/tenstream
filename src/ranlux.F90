module m_ranlux

  !implicit none

  private
  public :: ranlux, rluxgo

  dimension SEEDS(24), ISEEDS(24)
  parameter(MAXLEV=4, LXDFLT=3)
  dimension NDSKIP(0:MAXLEV)
  dimension NEXT(24)
  parameter(TWOP12=4096., IGIGA=1000000000, JSDFLT=314159265)
  parameter(ITWO24=2**24, ICONS=2147483563)
  save NOTYET, I24, J24, CARRY, SEEDS, TWOM24, TWOM12, LUXLEV
  save NSKIP, NDSKIP, IN24, NEXT, KOUNT, MKOUNT, INSEED
  integer LUXLEV
  logical NOTYET
  data NOTYET, LUXLEV, IN24, KOUNT, MKOUNT/.true., LXDFLT, 0, 0, 0/
  data I24, J24, CARRY/24, 10, 0./
  !                               default
  !  Luxury Level   0     1     2   *3*    4
  data NDSKIP/0, 24, 73, 199, 365/
  !orresponds to p=24    48    97   223   389
  !     time factor 1     2     3     6    10   on slow workstation
  !                 1    1.5    2     3     5   on fast mainframe
  !
contains
  subroutine RANLUX(RVEC, LENV)
!         Subtract-and-borrow random number generator proposed by
!         Marsaglia and Zaman, implemented by F. James with the name
!         RCARRY in 1991, and later improved by Martin Luescher
!         in 1993 to produce "Luxury Pseudorandom Numbers".
!     Fortran 77 coded by F. James, 1993
!
!       references:
!  M. Luscher, Computer Physics Communications  79 (1994) 100
!  F. James, Computer Physics Communications 79 (1994) 111
!
!   LUXURY LEVELS.
!   ------ ------      The available luxury levels are:
!
!  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
!           and Zaman, very long period, but fails many tests.
!  level 1  (p=48): considerable improvement in quality over level 0,
!           now passes the gap test, but still fails spectral test.
!  level 2  (p=97): passes all known tests, but theoretically still
!           defective.
!  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
!           correlations have very small chance of being observed.
!  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
!
!!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!!  Calling sequences for RANLUX:                                  ++
!!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
!!!!                   32-bit random floating point numbers between  ++
!!!!                   zero (not included) and one (also not incl.). ++
!!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
!!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
!!!!               which is integer between zero and MAXLEV, or if   ++
!!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
!!!!               should be set to zero unless restarting at a break++
!!!!               point given by output of RLUXAT (see RLUXAT).     ++
!!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
!!!!               which can be used to restart the RANLUX generator ++
!!!!               at the current point by calling RLUXGO.  K1 and K2++
!!!!               specify how many numbers were generated since the ++
!!!!               initialization with LUX and INT.  The restarting  ++
!!!!               skips over  K1+K2*E9   numbers, so it can be long.++
!!!!   A more efficient but less convenient way of restarting is by: ++
!!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
!!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
!!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
!!!!                 32-bit integer seeds, to be used for restarting ++
!!!!      ISVEC must be dimensioned 25 in the calling program        ++
!!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    dimension RVEC(LENV)
!  NOTYET is .TRUE. if no initialization has been performed yet.
!              Default Initialization by Multiplicative Congruential
    if (NOTYET) then
      NOTYET = .false.
      JSEED = JSDFLT
      INSEED = JSEED
      write (6, '(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ', JSEED
      LUXLEV = LXDFLT
      NSKIP = NDSKIP(LUXLEV)
      LP = NSKIP + 24
      IN24 = 0
      KOUNT = 0
      MKOUNT = 0
      write (6, '(A,I2,A,I4)') ' RANLUX DEFAULT LUXURY LEVEL =  ', &
        +LUXLEV, '      p =', LP
      TWOM24 = 1.
      do 25 I = 1, 24
        TWOM24 = TWOM24 * 0.5
        K = JSEED / 53668
        JSEED = 40014 * (JSEED - K * 53668) - K * 12211
        if (JSEED .lt. 0) JSEED = JSEED + ICONS
        ISEEDS(I) = mod(JSEED, ITWO24)
25      continue
        TWOM12 = TWOM24 * 4096.
        do 50 I = 1, 24
          SEEDS(I) = real(ISEEDS(I)) * TWOM24
          NEXT(I) = I - 1
50        continue
          NEXT(1) = 24
          I24 = 24
          J24 = 10
          CARRY = 0.
          if (abs(SEEDS(24)) .lt. tiny(SEEDS)) CARRY = TWOM24
          end if
!
!          The Generator proper: "Subtract-with-borrow",
!          as proposed by Marsaglia and Zaman,
!          Florida State University, March, 1989
!
          do 100 IVEC = 1, LENV
            UNI = SEEDS(J24) - SEEDS(I24) - CARRY
            if (UNI .lt. 0.) then
              UNI = UNI + 1.0
              CARRY = TWOM24
            else
              CARRY = 0.
            end if
            SEEDS(I24) = UNI
            I24 = NEXT(I24)
            J24 = NEXT(J24)
            RVEC(IVEC) = UNI
!  small numbers (with less than 12 "significant" bits) are "padded".
            if (UNI .lt. TWOM12) then
              RVEC(IVEC) = RVEC(IVEC) + TWOM24 * SEEDS(J24)
!        and zero is forbidden in case someone takes a logarithm
              if (abs(RVEC(IVEC)) .lt. tiny(rvec)) RVEC(IVEC) = TWOM24 * TWOM24
            end if
!        Skipping to luxury.  As proposed by Martin Luscher.
            IN24 = IN24 + 1
            if (IN24 .eq. 24) then
              IN24 = 0
              KOUNT = KOUNT + NSKIP
              do 90 ISK = 1, NSKIP
                UNI = SEEDS(J24) - SEEDS(I24) - CARRY
                if (UNI .lt. 0.) then
                  UNI = UNI + 1.0
                  CARRY = TWOM24
                else
                  CARRY = 0.
                end if
                SEEDS(I24) = UNI
                I24 = NEXT(I24)
                J24 = NEXT(J24)
90              continue
                end if
100             continue
                KOUNT = KOUNT + LENV
                if (KOUNT .ge. IGIGA) then
                  MKOUNT = MKOUNT + 1
                  KOUNT = KOUNT - IGIGA
                end if
                return
                end subroutine RANLUX
!
!           Entry to input and float integer seeds from previous run
                subroutine RLUXIN(ISDEXT)
                  dimension ISDEXT(25)
!--------------------------added to fix `bug'---------------------------
                  if (NOTYET) then
!*****$#--1----_----2----_----3----_----4----_----5----_----6----_----7--
!         WRITE(6,'(A)')  ' PROPER RESULTS ONLY WITH INITIALISATION FROM
!     $25 INTEGERS OBTAINED WITH RLUXUT'
                    NOTYET = .false.
!-----------------------------------------------------------------------
                  end if
                  TWOM24 = 1.
                  do I = 1, 24
                    NEXT(I) = I - 1
                    TWOM24 = TWOM24 * 0.5
                  end do
                  NEXT(1) = 24
                  TWOM12 = TWOM24 * 4096.
                  write (6, '(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
                  write (6, '(5X,5I12)') ISDEXT
                  do 200 I = 1, 24
                    SEEDS(I) = real(ISDEXT(I)) * TWOM24
200                 continue
                    CARRY = 0.
                    if (ISDEXT(25) .lt. 0) CARRY = TWOM24
                    ISD = IABS(ISDEXT(25))
                    I24 = mod(ISD, 100)
                    ISD = ISD / 100
                    J24 = mod(ISD, 100)
                    ISD = ISD / 100
                    IN24 = mod(ISD, 100)
                    ISD = ISD / 100
                    LUXLEV = ISD
                    if (LUXLEV .le. MAXLEV) then
                      NSKIP = NDSKIP(LUXLEV)
                      write (6, '(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ', &
                        +LUXLEV
                    else if (LUXLEV .ge. 24) then
                      NSKIP = LUXLEV - 24
                      write (6, '(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:', LUXLEV
                    else
                      NSKIP = NDSKIP(MAXLEV)
                      write (6, '(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ', LUXLEV
                      LUXLEV = MAXLEV
                    end if
                    INSEED = -1
                    end subroutine RLUXIN
!
!                    Entry to ouput seeds as integers
                    subroutine RLUXUT(ISDEXT)
                      dimension ISDEXT(25)
                      do I = 1, 24
                        ISDEXT(I) = int(SEEDS(I) * TWOP12 * TWOP12)
                      end do
                      ISDEXT(25) = I24 + 100 * J24 + 10000 * IN24 + 1000000 * LUXLEV
                      if (CARRY .gt. 0.) ISDEXT(25) = -ISDEXT(25)
                    end subroutine RLUXUT
!
!                    Entry to output the "convenient" restart point
                    subroutine RLUXAT(LOUT, INOUT, K1, K2)
                      LOUT = LUXLEV
                      INOUT = INSEED
                      K1 = KOUNT
                      K2 = MKOUNT
                    end subroutine RLUXAT
!
!                    Entry to initialize from one or three integers
                    subroutine RLUXGO(LUX, INS, K1, K2)
                      if (LUX .lt. 0) then
                        LUXLEV = LXDFLT
                      else if (LUX .le. MAXLEV) then
                        LUXLEV = LUX
                      else if (LUX .lt. 24 .or. LUX .gt. 2000) then
                        LUXLEV = MAXLEV
                        write (6, '(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ', LUX
                      else
                        LUXLEV = LUX
                        do 310 ILX = 0, MAXLEV
                          if (LUX .eq. NDSKIP(ILX) + 24) LUXLEV = ILX
310                       continue
                          end if
                          if (LUXLEV .le. MAXLEV) then
                            NSKIP = NDSKIP(LUXLEV)
!         WRITE(6,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',
!     +        LUXLEV,'     P=', NSKIP+24
                          else
                            NSKIP = LUXLEV - 24
                            write (6, '(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:', LUXLEV
                          end if
                          IN24 = 0
                          if (INS .lt. 0) write (6, '(A)') ' Illegal initialization by RLUXGO, negative input seed'
                          if (INS .gt. 0) then
                            JSEED = INS
!        WRITE(6,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',
!     +      JSEED, K1,K2
                          else
                            JSEED = JSDFLT
                            write (6, '(A)') ' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
                          end if
                          INSEED = JSEED
                          NOTYET = .false.
                          TWOM24 = 1.
                          do 325 I = 1, 24
                            TWOM24 = TWOM24 * 0.5
                            K = JSEED / 53668
                            JSEED = 40014 * (JSEED - K * 53668) - K * 12211
                            if (JSEED .lt. 0) JSEED = JSEED + ICONS
                            ISEEDS(I) = mod(JSEED, ITWO24)
325                         continue
                            TWOM12 = TWOM24 * 4096.
                            do 350 I = 1, 24
                              SEEDS(I) = real(ISEEDS(I)) * TWOM24
                              NEXT(I) = I - 1
350                           continue
                              NEXT(1) = 24
                              I24 = 24
                              J24 = 10
                              CARRY = 0.
                              if (abs(SEEDS(24)) .lt. tiny(seeds)) CARRY = TWOM24
!        If restarting at a break point, skip K1 + IGIGA*K2
!        Note that this is the number of numbers delivered to
!        the user PLUS the number skipped (if luxury .GT. 0).
                              KOUNT = K1
                              MKOUNT = K2
                              if (K1 + K2 .ne. 0) then
                                do 500 IOUTER = 1, K2 + 1
                                  INNER = IGIGA
                                  if (IOUTER .eq. K2 + 1) INNER = K1
                                  do 450 ISK = 1, INNER
                                    UNI = SEEDS(J24) - SEEDS(I24) - CARRY
                                    if (UNI .lt. 0.) then
                                      UNI = UNI + 1.0
                                      CARRY = TWOM24
                                    else
                                      CARRY = 0.
                                    end if
                                    SEEDS(I24) = UNI
                                    I24 = NEXT(I24)
                                    J24 = NEXT(J24)
450                                 continue
500                                 continue
!         Get the right value of IN24 by direct calculation
                                    IN24 = mod(KOUNT, NSKIP + 24)
                                    if (MKOUNT .gt. 0) then
                                      IZIP = mod(IGIGA, NSKIP + 24)
                                      IZIP2 = MKOUNT * IZIP + IN24
                                      IN24 = mod(IZIP2, NSKIP + 24)
                                    end if
!       Now IN24 had better be between zero and 23 inclusive
                                    if (IN24 .gt. 23) then
                                      write (6, '(A/A,3I11,A,I5)') '  Error in RESTARTING with RLUXGO:',&
                                        & '  The values', INS, K1, K2, ' cannot occur at luxury level', LUXLEV
                                      IN24 = 0
                                    end if
                                    end if
                                    end subroutine RLUXGO
                                    end module
