!=====================================================================
!
!  Routines from "Numerical Recipes: the Art of Scientific Computing"
!  W. H. Press et al., Cambridge University Press
!
!=====================================================================

! double precision routines

  double precision function erf(x)

  implicit none

  double precision x

! this routine uses routine gammp
  double precision gammp

  if(x<0.)then
    erf=-gammp(0.5d0,x**2)
  else
    erf=gammp(0.5d0,x**2)
  endif

  end function erf

! ---------------------------------

  double precision function gammp(a,x)

  implicit none

  double precision a,x

! this routine uses routines gcf and gser
  double precision gammcf,gamser,gln

  if(x<0.d0 .or. a <= 0.d0) stop 'bad arguments in gammp'

  if(x<a+1.d0)then
    call gser(gamser,a,x,gln)
    gammp=gamser
  else
    call gcf(gammcf,a,x,gln)
    gammp=1.d0-gammcf
  endif

  end function gammp

! ---------------------------------

  subroutine gcf(gammcf,a,x,gln)

  implicit none

  double precision a,gammcf,gln,x

  double precision, parameter :: EPS=3.d-7,FPMIN=1.d-30
  integer, parameter :: ITMAX=100

! this routine uses routine gammln

  integer i
  double precision an,b,c,d,del,h

  double precision, external :: gammln

  gln=gammln(a)
  b=x+1.d0-a
  c=1.d0/FPMIN
  d=1.d0/b
  h=d
  do i=1,ITMAX
    an=-i*(i-a)
    b=b+2.d0
    d=an*d+b
    if(dabs(d)<FPMIN)d=FPMIN
    c=b+an/c
    if(dabs(c)<FPMIN)c=FPMIN
    d=1.d0/d
    del=d*c
    h=h*del
    if(dabs(del-1.d0)<EPS) then
      gammcf=exp(-x+a*log(x)-gln)*h
      return
    endif
  enddo

  stop 'a too large, ITMAX too small in gcf'

  end subroutine gcf

! ---------------------------------

  subroutine gser(gamser,a,x,gln)

  implicit none

  double precision a,gamser,gln,x

  integer, parameter :: ITMAX=100
  double precision, parameter :: EPS=3.d-7

! this routine uses routine gammln

  integer n
  double precision ap,del,sumval

  double precision, external :: gammln

  gln=gammln(a)

  if(x <= 0.d0)then
    if(x<0.d0) stop 'x < 0 in gser'
    gamser=0.d0
    return
  endif

  ap=a
  sumval=1.d0/a
  del=sumval

  do n=1,ITMAX
    ap=ap+1.d0
    del=del*x/ap
    sumval=sumval+del
    if(dabs(del)<dabs(sumval)*EPS) then
      gamser=sumval*exp(-x+a*log(x)-gln)
      return
    endif
  enddo

  stop 'a too large, ITMAX too small in gser'

  end subroutine gser

! ---------------------------------

  double precision function gammln(xx)

  implicit none

  double precision xx

  integer j
  double precision ser,stp,tmp,x,y,cof(6)

  cof(1) = 76.18009172947146d0
  cof(2) = -86.50532032941677d0
  cof(3) = 24.01409824083091d0
  cof(4) = -1.231739572450155d0
  cof(5) = 0.1208650973866179d-2
  cof(6) = -0.5395239384953d-5

  stp = 2.5066282746310005d0

  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
  do j=1,6
    y=y+1.d0
    ser=ser+cof(j)/y
  enddo
  gammln=tmp+log(stp*ser/x)

  end function gammln

