
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.0
!                   ------------------------------
!
!                         Dimitri Komatitsch
!          Universite de Pau et des Pays de l'Adour, France
!
!                          (c) May 2004
!
!========================================================================

  subroutine qmasspec(rhoext,wxgll,wygll,ibool,dvolu,rmass,density,kmato,npoin,ireadmodel,nspec,numat)

! build the mass matrix

  implicit none

  include "constants.h"

  integer npoin,nspec,numat

  double precision wxgll(NGLLX),wygll(NGLLY),rmass(npoin),dvolu(nspec,NGLLX,NGLLX),density(numat)
  double precision rhoext(npoin)

  integer kmato(nspec),ibool(NGLLX,NGLLX,nspec)

  integer numelem,material,i,j,iglobnum
  logical ireadmodel
  double precision denst

!
!----  compute the mass matrix by summing the contribution of each point
!

  rmass = zero

  do numelem = 1,nspec

  material = kmato(numelem)
  denst    = density(material)

  do i=1,NGLLX
       do j=1,NGLLY

  iglobnum = ibool(i,j,numelem)

!--- si formulation heterogene pour un modele de densite externe
  if(ireadmodel) denst = rhoext(iglobnum)

  rmass(iglobnum) = rmass(iglobnum) + denst * wxgll(i) * wygll(j) * dvolu(numelem,i,j)

      enddo
   enddo

  enddo

  end subroutine qmasspec

