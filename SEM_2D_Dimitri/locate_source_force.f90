
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.1
!                   ------------------------------
!
!                         Dimitri Komatitsch
!          Universite de Pau et des Pays de l'Adour, France
!
!                          (c) January 2005
!
!========================================================================

  subroutine locate_source_force(coord,ibool,npoin,nspec,x_source,z_source,source_type,ix_source,iz_source,ispec_source,iglob_source)

!
!----- calculer la position reelle de la source
!

  implicit none

  include "constants.h"

  integer npoin,nspec,source_type
  integer ibool(NGLLX,NGLLZ,nspec)

  double precision x_source,z_source
  double precision coord(NDIM,npoin)

  integer ip,ix,iz,numelem,ilowx,ilowz,ihighx,ihighz,ix_source,iz_source,ispec_source,iglob_source

  double precision distminmax,distmin,xp,zp,dist

  write(iout,200)

  distminmax = -HUGEVAL

      distmin = +HUGEVAL

      ilowx = 1
      ilowz = 1
      ihighx = NGLLX
      ihighz = NGLLZ

! on ne fait la recherche que sur l'interieur de l'element si source explosive
  if(source_type == 2) then
    ilowx = 2
    ilowz = 2
    ihighx = NGLLX-1
    ihighz = NGLLZ-1
  endif

! recherche du point de grille le plus proche
      do numelem=1,nspec
      do ix=ilowx,ihighx
      do iz=ilowz,ihighz

! numero global du point
        ip=ibool(ix,iz,numelem)

! coordonnees du point de grille
            xp = coord(1,ip)
            zp = coord(2,ip)

            dist = sqrt((xp-x_source)**2 + (zp-z_source)**2)

! retenir le point pour lequel l'ecart est minimal
            if(dist < distmin) then
              distmin = dist
              iglob_source = ip
              ix_source = ix
              iz_source = iz
              ispec_source = numelem
            endif

      enddo
      enddo
      enddo

  distminmax = max(distmin,distminmax)

  write(iout,150) x_source,z_source,coord(1,iglob_source),coord(2,iglob_source),distmin
  write(iout,160) distminmax

 150 format(1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3,f12.3)
 160 format(/2x,'Maximum distance between asked and real =',f12.3)
 200 format(//1x,48('=')/,' =  S o u r c e s  ', &
  'r e a l  p o s i t i o n s  ='/1x,48('=')// &
  '    Source    x-asked      z-asked     x-obtain     z-obtain       dist'/)

  end subroutine locate_source_force

