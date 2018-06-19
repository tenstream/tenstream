module m_kiss_rng
  ! https://www.thecodingforums.com/threads/64-bit-kiss-rngs.673657/

  use iso_fortran_env, only: real32, real64, int64
  implicit none
  interface kiss_real
    module procedure kiss_real32, kiss_real64
  end interface
contains
  subroutine kiss_real64(r)
    real(real64), intent(out) :: r
    r = rkiss()
  end subroutine
  subroutine kiss_real32(r)
    real(real32), intent(out) :: r
    r = rkiss()
  end subroutine

  !function KISS()
  !  implicit integer*8(a-z)
  !  data x,y,z,c /1234567890987654321_8, 362436362436362436_8,&
  !    1066149217761810_8, 123456123456123456_8/
  !  save x,y,z,c
  !  m(x,k)=ieor(x,ishft(x,k)) !statement function
  !  s(x)=ishft(x,-63) !statement function
  !  t=ishft(x,58)+c
  !  if(s(x).eq.s(t)) then; c=ishft(x,-6)+s(x)
  !  else; c=ishft(x,-6)+1-s(x+t); endif
  !    x=t+x
  !    y=m(m(m(y,13_8),-17_8),43_8)
  !    z=6906969069_8*z+1234567
  !    KISS=x+y+z
  !    return
  !end function

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Random number generator KISS05 after a suggestion by George Marsaglia
  ! in "Random numbers for C: The END?" posted on sci.crypt.random-numbers
  ! in 1999
  !
  ! version as in "double precision RNGs" in  sci.math.num-analysis
  ! http://sci.tech-archive.net/Archive/sci.math.num-analysis/2005-11/msg00352.html
  !
  ! The  KISS (Keep It Simple Stupid) random number generator. Combines:
  ! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
  ! (2) A 3-shift shift-register generator, period 2^32-1,
  ! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
  ! Overall period > 2^123
  !
  !
  ! A call to rkiss05() gives one random real in the interval [0,1),
  ! i.e., 0 <= rkiss05 < 1
  !
  ! Before using rkiss05 call kissinit(seed) to initialize
  ! the generator by random integers produced by Park/Millers
  ! minimal standard LCG.
  ! Seed should be any positive integer.
  !
  ! FORTRAN implementation by Thomas Vojta, vojta@mst.edu
  ! built on a module found at www.fortran.com
  !
  !
  ! History:
  !        v0.9     Dec 11, 2010    first implementation
  !        V0.91    Dec 11, 2010    inlined internal function for the SR component
  !        v0.92    Dec 13, 2010    extra shuffle of seed in kissinit
  !        v093     Aug 13, 2012    changed inter representation test to avoid data statements
  !
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine kiss_init(iinit)
    implicit none
    integer,parameter :: r8b= selected_real_kind(p=14,r=99)   ! 8-byte reals
    integer,parameter :: i4b= selected_int_kind(8)            ! 4-byte integers

    integer(i4b) idum,ia,im,iq,ir,iinit
    integer(i4b) k,x,y,z,w,c1,c2,c3,c4
    real(r8b)    rkiss05,rdum
    parameter (ia=16807,im=2147483647,iq=127773,ir=2836)
    common /kisscom/x,y,z,w

    !!! test integer representation !!!
    c1=-8
    c1=ishftc(c1,-3)
    !     print *,c1
    if (c1.ne.536870911) then
      print *,'nonstandard integer representation. stopped.'
      stop
    endif

    idum=iinit
    idum= abs(1099087573 * idum)               ! 32-bit lcg to shuffle seeds
    if (idum.eq.0) idum=1
    if (idum.ge.im) idum=im-1

    k=(idum)/iq
    idum=ia*(idum-k*iq)-ir*k
    if (idum.lt.0) idum = idum + im
    if (idum.lt.1) then
      x=idum+1
    else
      x=idum
    endif
    k=(idum)/iq
    idum=ia*(idum-k*iq)-ir*k
    if (idum.lt.0) idum = idum + im
    if (idum.lt.1) then
      y=idum+1
    else
      y=idum
    endif
    k=(idum)/iq
    idum=ia*(idum-k*iq)-ir*k
    if (idum.lt.0) idum = idum + im
    if (idum.lt.1) then
      z=idum+1
    else
      z=idum
    endif
    k=(idum)/iq
    idum=ia*(idum-k*iq)-ir*k
    if (idum.lt.0) idum = idum + im
    if (idum.lt.1) then
      w=idum+1
    else
      w=idum
    endif

    rdum = rkiss()

    return
  end subroutine

  function rkiss()
    implicit none

    integer,parameter      :: r8b= selected_real_kind(p=14,r=99)   ! 8-byte reals
    integer,parameter      :: i4b= selected_int_kind(8)            ! 4-byte integers
    real(r8b),parameter    :: am = 4.656612873077392578d-10        ! multiplier 1/2^31

    real(r8b)             :: rkiss
    integer(i4b)          :: kiss
    integer(i4b)          :: x,y,z,w              ! working variables for the four generators
    common /kisscom/x,y,z,w

    x = 69069 * x + 1327217885
    y= ieor (y, ishft (y, 13)); y= ieor (y, ishft (y, -17)); y= ieor (y, ishft (y, 5))
    z = 18000 * iand (z, 65535) + ishft (z, - 16)
    w = 30903 * iand (w, 65535) + ishft (w, - 16)
    kiss = ishft(x + y + ishft (z, 16) + w , -1)
    rkiss=kiss*am
  end function
end module
