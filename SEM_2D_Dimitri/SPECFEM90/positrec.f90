!=====================================================================
!
!                 S p e c f e m  V e r s i o n  4 . 2
!                 -----------------------------------
!
!                         Dimitri Komatitsch
!    Department of Earth and Planetary Sciences - Harvard University
!                         Jean-Pierre Vilotte
!                 Departement de Sismologie - IPGP - Paris
!                           (c) June 1998
!
!=====================================================================

  subroutine positrec(coord,posrec,ndime,npoin,nrec)

!
!---- calculer la position reelle des recepteurs
!

  use iounit

  implicit none

  integer ndime,npoin,nrec
  double precision coord(ndime,npoin)
  double precision posrec(ndime,nrec)

  double precision dminmax,dmin,xs,zs,xp,zp,dist
  integer n,ip,ipoint

  write(iout,200)

  dminmax = -1.d30

  do n=1,nrec

      dmin = +1.d30

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

  return
  end subroutine positrec
