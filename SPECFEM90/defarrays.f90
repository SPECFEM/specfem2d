
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
          xigll,zigll,xix,xiz,gammax,gammaz,a11,a12, &
          ibool,kmato,coord,npoin,rsizemin,rsizemax, &
          cpoverdxmin,cpoverdxmax,lambdaSmin,lambdaSmax,lambdaPmin,lambdaPmax, &
          vpmin,vpmax,readmodel,nspec,numat,source_type,ix_source,iz_source,ispec_source)

! define all the arrays for the variational formulation

  implicit none

  include "constants.h"

  integer i,j,ispec,material,ipointnum,npoin,nspec,numat
  integer ix_source,iz_source,ispec_source,ir,is,source_type

  integer kmato(nspec),ibool(NGLLX,NGLLX,nspec)

  double precision xix(NGLLX,NGLLZ,nspec)
  double precision xiz(NGLLX,NGLLZ,nspec)
  double precision gammax(NGLLX,NGLLZ,nspec)
  double precision gammaz(NGLLX,NGLLZ,nspec)

  double precision density(numat),elastcoef(4,numat)

  double precision coord(NDIME,npoin)

  double precision a11(NGLLX,NGLLX),a12(NGLLX,NGLLX)

  double precision xigll(NGLLX),zigll(NGLLZ)

  double precision vpext(npoin)
  double precision vsext(npoin)
  double precision rhoext(npoin)

  double precision vsmin,vsmax,densmin,densmax
  double precision lambdaplus2mu,lambda,mu,denst
  double precision kappa,cploc,csloc,x0,z0
  double precision x1,z1,x2,z2,rdist1,rdist2,rapportmin,rapportmax
  double precision lambdamin,lambdamax
  double precision flagxprime,flagzprime,sig0

  double precision rsizemin,rsizemax,cpoverdxmin,cpoverdxmax, &
    lambdaSmin,lambdaSmax,lambdaPmin,lambdaPmax,vpmin,vpmax

  logical readmodel

  double precision, external :: lagrange_deriv_GLL

!
!-----------------------------------------------------------------------
!

!---- compute parameters for the spectral elements

  a11 = zero
  a12 = zero

  vpmin = HUGEVAL
  vsmin = HUGEVAL
  vpmax = -HUGEVAL
  vsmax = -HUGEVAL
  densmin = HUGEVAL
  densmax = -HUGEVAL

  rsizemin = HUGEVAL
  rsizemax = -HUGEVAL

  cpoverdxmin = HUGEVAL
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
  if(readmodel) then
    ipointnum = ibool(i,j,ispec)
    cploc = vpext(ipointnum)
    csloc = vsext(ipointnum)
    denst = rhoext(ipointnum)
    mu   = denst*csloc*csloc
    lambda  = denst*cploc*cploc - 2.d0*mu
    lambdaplus2mu  = lambda + 2.d0*mu
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
    cpoverdxmin = min(cpoverdxmin,rapportmin)
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

  print *
  print *,'********'
  print *,'Modele : vitesse P min,max = ',vpmin,vpmax
  print *,'Modele : vitesse S min,max = ',vsmin,vsmax
  print *,'Modele : densite min,max = ',densmin,densmax
  print *,'********'
  print *

! seulement si source explosive
  if(source_type == 2) then

  if(ix_source == 1 .or. ix_source == NGLLX .or. iz_source == 1 .or. iz_source == NGLLX) &
        stop 'Explosive source on element edge'

!---- definir a11 et a12 - dirac (schema en croix)

  sig0 = one

  do ir=1,NGLLX
    flagxprime = lagrange_deriv_GLL(ir-1,ix_source-1,xigll,NGLLX)
    a11(ir,iz_source) = a11(ir,iz_source) + sig0*xix(ix_source,iz_source,ispec_source)*flagxprime
    a12(ir,iz_source) = a12(ir,iz_source) + sig0*xiz(ix_source,iz_source,ispec_source)*flagxprime
  enddo

  do is=1,NGLLZ
    flagzprime = lagrange_deriv_GLL(is-1,iz_source-1,zigll,NGLLZ)
    a11(ix_source,is) = a11(ix_source,is) + sig0*gammax(ix_source,iz_source,ispec_source)*flagzprime
    a12(ix_source,is) = a12(ix_source,is) + sig0*gammaz(ix_source,iz_source,ispec_source)*flagzprime
  enddo

  endif

  end subroutine defarrays

