
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

  subroutine defarrays(vpext,vsext,rhoext,density,elastcoef, &
          ibool,kmato,coord,npoin,rsizemin,rsizemax, &
          cpoverdxmax,lambdaSmin,lambdaSmax,lambdaPmin,lambdaPmax, &
          vpmin,vpmax,read_external_model,nspec,numat)

! define all the arrays for the variational formulation

  implicit none

  include "constants.h"

  integer i,j,ispec,material,ipointnum,npoin,nspec,numat

  integer kmato(nspec),ibool(NGLLX,NGLLX,nspec)

  double precision density(numat),elastcoef(4,numat)

  double precision coord(NDIM,npoin)

  double precision vpext(npoin)
  double precision vsext(npoin)
  double precision rhoext(npoin)

  double precision vsmin,vsmax,densmin,densmax
  double precision lambdaplus2mu,lambda,mu,denst
  double precision kappa,cploc,csloc,x0,z0
  double precision x1,z1,x2,z2,rdist1,rdist2,rapportmin,rapportmax
  double precision lambdamin,lambdamax

  double precision rsizemin,rsizemax,cpoverdxmax, &
    lambdaSmin,lambdaSmax,lambdaPmin,lambdaPmax,vpmin,vpmax

  logical read_external_model

!
!-----------------------------------------------------------------------
!

!---- compute parameters for the spectral elements

  vpmin = HUGEVAL
  vsmin = HUGEVAL
  vpmax = -HUGEVAL
  vsmax = -HUGEVAL
  densmin = HUGEVAL
  densmax = -HUGEVAL

  rsizemin = HUGEVAL
  rsizemax = -HUGEVAL

  cpoverdxmax = -HUGEVAL

  lambdaPmin = HUGEVAL
  lambdaSmin = HUGEVAL
  lambdaPmax = -HUGEVAL
  lambdaSmax = -HUGEVAL

  do ispec=1,nspec

    material = kmato(ispec)

    lambda = elastcoef(1,material)
    mu = elastcoef(2,material)
    lambdaplus2mu  = elastcoef(3,material)
    denst = density(material)

    kappa = lambda + 2.d0*mu/3.d0

    cploc = sqrt((kappa + 4.d0*mu/3.d0)/denst)
    csloc = sqrt(mu/denst)

  do j=1,NGLLZ
    do i=1,NGLLX

!--- si formulation heterogene pour un modele de vitesse externe
  if(read_external_model) then
    ipointnum = ibool(i,j,ispec)
    cploc = vpext(ipointnum)
    csloc = vsext(ipointnum)
    denst = rhoext(ipointnum)
    mu = denst*csloc*csloc
    lambda = denst*cploc*cploc - 2.d0*mu
    lambdaplus2mu = lambda + 2.d0*mu
  endif

!--- calculer min et max du modele de vitesse et densite
  vpmin = min(vpmin,cploc)
  vpmax = max(vpmax,cploc)

  vsmin = min(vsmin,csloc)
  vsmax = max(vsmax,csloc)

  densmin = min(densmin,denst)
  densmax = max(densmax,denst)

!--- stocker parametres pour verifs diverses
  if(i < NGLLX .and. j < NGLLZ) then

    x0 = coord(1,ibool(i,j,ispec))
    z0 = coord(2,ibool(i,j,ispec))
    x1 = coord(1,ibool(i+1,j,ispec))
    z1 = coord(2,ibool(i+1,j,ispec))
    x2 = coord(1,ibool(i,j+1,ispec))
    z2 = coord(2,ibool(i,j+1,ispec))

    rdist1 = sqrt((x1-x0)**2 + (z1-z0)**2)
    rdist2 = sqrt((x2-x0)**2 + (z2-z0)**2)

    rsizemin = min(rsizemin,rdist1)
    rsizemin = min(rsizemin,rdist2)
    rsizemax = max(rsizemax,rdist1)
    rsizemax = max(rsizemax,rdist2)

    rapportmin = cploc / max(rdist1,rdist2)
    rapportmax = cploc / min(rdist1,rdist2)
    cpoverdxmax = max(cpoverdxmax,rapportmax)

    x0 = coord(1,ibool(1,1,ispec))
    z0 = coord(2,ibool(1,1,ispec))
    x1 = coord(1,ibool(NGLLX,1,ispec))
    z1 = coord(2,ibool(NGLLX,1,ispec))
    x2 = coord(1,ibool(1,NGLLZ,ispec))
    z2 = coord(2,ibool(1,NGLLZ,ispec))

    rdist1 = sqrt((x1-x0)**2 + (z1-z0)**2)
    rdist2 = sqrt((x2-x0)**2 + (z2-z0)**2)

    lambdamin = cploc/max(rdist1,rdist2)
    lambdamax = cploc/min(rdist1,rdist2)
    lambdaPmin = min(lambdaPmin,lambdamin)
    lambdaPmax = max(lambdaPmax,lambdamax)

    lambdamin = csloc/max(rdist1,rdist2)
    lambdamax = csloc/min(rdist1,rdist2)
    lambdaSmin = min(lambdaSmin,lambdamin)
    lambdaSmax = max(lambdaSmax,lambdamax)

  endif

    enddo
 enddo
  enddo

  write(IOUT,*)
  write(IOUT,*) '********'
  write(IOUT,*) 'Modele : vitesse P min,max = ',vpmin,vpmax
  write(IOUT,*) 'Modele : vitesse S min,max = ',vsmin,vsmax
  write(IOUT,*) 'Modele : densite min,max = ',densmin,densmax
  write(IOUT,*) '********'
  write(IOUT,*)

  end subroutine defarrays

