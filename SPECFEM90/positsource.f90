
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

  subroutine positsource(coord,ibool,gltfu,npoin,nspec)

!
!----- calculer la position reelle de la source
!

  implicit none

  include "constants.h"

  integer npoin,nspec
  double precision coord(NDIME,npoin)
  double precision gltfu(20)
  integer ibool(NGLLX,NGLLZ,nspec)

  double precision dminmax,dmin,xs,zs,xp,zp,dist
  integer ip,ipoint,ix,iy,numelem,ilowx,ilowy,ihighx,ihighy

  write(iout,200)

  dminmax = -HUGEVAL

      dmin = +HUGEVAL

! coordonnees demandees pour la source
      xs = gltfu(3)
      zs = gltfu(4)

      ilowx = 1
      ilowy = 1
      ihighx = NGLLX
      ihighy = NGLLZ

! on ne fait la recherche que sur l'interieur de l'element si source explosive
  if(nint(gltfu(2)) == 2) then
    ilowx = 2
    ilowy = 2
    ihighx = NGLLX-1
    ihighy = NGLLZ-1
  endif

! recherche du point de grille le plus proche
      do numelem=1,nspec
      do ix=ilowx,ihighx
      do iy=ilowy,ihighy

! numero global du point
        ip=ibool(ix,iy,numelem)

! coordonnees du point de grille
            xp = coord(1,ip)
            zp = coord(2,ip)

            dist = dsqrt((xp-xs)**2 + (zp-zs)**2)

! retenir le point pour lequel l'ecart est minimal
            if (dist < dmin) then
              dmin = dist
              gltfu(9) = ip
              gltfu(10) = ix
              gltfu(11) = iy
              gltfu(12) = numelem
            endif

      enddo
      enddo
      enddo

  ipoint = nint(gltfu(9))

  dminmax = dmax1(dmin,dminmax)

  write(iout,150) xs,zs,coord(1,ipoint),coord(2,ipoint),dmin
  write(iout,160) dminmax

 150 format(1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3,f12.3)
 160 format(/2x,'Maximum distance between asked and real =',f12.3)
 200 format(//1x,48('=')/,' =  S o u r c e s  ', &
  'r e a l  p o s i t i o n s  ='/1x,48('=')// &
  '    Source    x-asked      z-asked     x-obtain     z-obtain       dist'/)

  end subroutine positsource

