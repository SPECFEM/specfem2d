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

  subroutine positsource(coord,ibool,gltfu,ndime,npoin,nltfl,nxgll,nygll,nspec)

!
!----- calculer la position reelle des sources
!

  use iounit

  implicit none

  integer ndime,npoin,nltfl,nxgll,nygll,nspec
  double precision coord(ndime,npoin)
  double precision gltfu(20,nltfl)
  integer ibool(0:nxgll-1,0:nygll-1,nspec)

  double precision dminmax,dmin,xs,zs,xp,zp,dist
  integer n,ip,ipoint,ix,iy,numelem,ilowx,ilowy,ihighx,ihighy

  write(iout,200)

  dminmax = -1.d30

  do n=1,nltfl

      dmin = +1.d30

! coordonnees demandees pour la source
      xs = gltfu(3,n)
      zs = gltfu(4,n)

      ilowx = 0
      ilowy = 0
      ihighx = nxgll-1
      ihighy = nygll-1

! on ne fait la recherche que sur l'interieur de l'element si source explosive
  if(nint(gltfu(2,n)) == 2) then
      ilowx = 1
      ilowy = 1
      ihighx = nxgll-2
      ihighy = nygll-2
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
                  gltfu(9,n) = ip
                  gltfu(10,n) = ix
                  gltfu(11,n) = iy
                  gltfu(12,n) = numelem
            endif

      enddo
      enddo
      enddo

  ipoint = nint(gltfu(9,n))

  dminmax = dmax1(dmin,dminmax)

  write(iout,150) n,xs,zs,coord(1,ipoint),coord(2,ipoint),dmin

  enddo

  write(iout,160) dminmax

 150   format(1x,i7,1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3,f12.3)
 160   format(/2x,'Maximum distance between asked and real =',f12.3)
 200  format(//1x,48('=')/,' =  S o u r c e s  ', &
  'r e a l  p o s i t i o n s  ='/1x,48('=')// &
  '    Source    x-asked      z-asked     ', &
  'x-obtain     z-obtain       dist'/)

  return
  end subroutine positsource
