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

  subroutine qmasspec(rhoext,wx,wy,ibool,dvolu,rmass,density,kmato,npoin)
!
!=======================================================================
!
!     "q m a s s p e c" : Build the mass matrix for the spectral
!                         elements
!
!=======================================================================
!
  use spela202
  use constspec

  implicit none

  integer npoin

  double precision wx(0:nxgll-1),wy(0:nygll-1),rmass(npoin), &
          dvolu(nspec,0:nxgll-1,0:nxgll-1),density(numat)
  double precision rhoext(npoin)

  integer kmato(nspec),ibool(0:nxgll-1,0:nxgll-1,nspec)

  integer numelem,material,i,j,iglobnum
  double precision denst

  double precision, parameter :: zero=0.d0, one=1.d0

!
!----  compute the mass matrix by summing the contribution of each point
!

  rmass = zero

  do numelem = 1,nspec

  material = kmato(numelem)
  denst    = density(material)

  do i=0,nxgll-1
       do j=0,nygll-1

  iglobnum = ibool(i,j,numelem)

!--- si formulation heterogene pour un modele de densite externe
  if(ireadmodel) denst = rhoext(iglobnum)

  rmass(iglobnum) = rmass(iglobnum) + &
                         denst * wx(i) * wy(j) * dvolu(numelem,i,j)

      enddo
   enddo

  enddo

!----  in case of periodic boundary conditions, fill the mass matrix
  where(rmass == zero) rmass = one

  return
  end subroutine qmasspec
