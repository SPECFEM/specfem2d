
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.1
!                   ------------------------------
!
!                         Dimitri Komatitsch
!          Universite de Pau et des Pays de l'Adour, France
!
!                          (c) December 2004
!
!========================================================================

  subroutine positrec(coord,posrec,npoin,nrec)

!
!---- calculer la position reelle des recepteurs
!

  implicit none

  include "constants.h"

  integer npoin,nrec
  double precision coord(NDIME,npoin)
  double precision posrec(NDIME,nrec)

  double precision dminmax,dmin,xs,zs,xp,zp,dist
  integer n,ip,ipoint

  write(iout,200)

  dminmax = -HUGEVAL

  do n=1,nrec

      dmin = +HUGEVAL

! coordonnees demandees
  xs = posrec(1,n)
  zs = posrec(2,n)

    do ip=1,npoin

! coordonnees du point de grille
      xp = coord(1,ip)
      zp = coord(2,ip)

      dist = dsqrt((xp-xs)**2 + (zp-zs)**2)

! retenir le point pour lequel l'ecart est minimal
      if (dist < dmin) then
        dmin = dist
        ipoint = ip
      endif

    enddo

    dminmax = dmax1(dmin,dminmax)

  write(iout,150) n,xs,zs,coord(1,ipoint),coord(2,ipoint),dmin

! stocker numero global dans premiere coordonnee
  posrec(1,n) = dble(ipoint)

  enddo

  write(iout,160) dminmax

 150   format(1x,i7,1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3,f12.3)
 160   format(/2x,'Maximum distance between asked and real =',f12.3)
 200  format(//1x,51('=')/,' =  R e c e i v e r s  ', &
  'r e a l  p o s i t i o n s  ='/1x,51('=')// &
  '  Receiver    x-asked      z-asked     ', &
  'x-obtain     z-obtain       dist'/)

  end subroutine positrec

