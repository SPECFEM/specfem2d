
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

  subroutine defarrays(vpext,vsext,rhoext,density,elastcoef, &
          xigll,zigll,xix,xiz,gammax,gammaz,a11,a12, &
          ibool,kmato,coord,gltfu,npoin,rsizemin,rsizemax, &
          cpoverdxmin,cpoverdxmax,rlambdaSmin,rlambdaSmax, &
          rlambdaPmin,rlambdaPmax,vpmin,vpmax,ireadmodel,nspec,numat)

! define all the arrays for the variational formulation

  implicit none

  include "constants.h"

  integer i,j,ispec,material,ipointnum,npoin,nspec,numat
  integer isourx,isourz,ielems,ir,is

  integer kmato(nspec),ibool(NGLLX,NGLLX,nspec)

  double precision xix(NGLLX,NGLLZ,nspec)
  double precision xiz(NGLLX,NGLLZ,nspec)
  double precision gammax(NGLLX,NGLLZ,nspec)
  double precision gammaz(NGLLX,NGLLZ,nspec)

  double precision density(numat),elastcoef(4,numat)

  double precision coord(NDIME,npoin)

  double precision a11(NGLLX,NGLLX),a12(NGLLX,NGLLX)

  double precision xigll(NGLLX),zigll(NGLLZ)

  double precision gltfu(20)
  double precision vpext(npoin)
  double precision vsext(npoin)
  double precision rhoext(npoin)

  double precision vsmin,vsmax,densmin,densmax
  double precision rKmod,rlambda,rmu,denst
  double precision rKvol,cploc,csloc,x0,z0
  double precision x1,z1,x2,z2,rdist1,rdist2,rapportmin,rapportmax
  double precision rlambmin,rlambmax
  double precision flagxprime,flagzprime,sig0

  double precision rsizemin,rsizemax,cpoverdxmin,cpoverdxmax, &
    rlambdaSmin,rlambdaSmax,rlambdaPmin,rlambdaPmax,vpmin,vpmax

  logical ireadmodel

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

  rlambdaPmin = HUGEVAL
  rlambdaSmin = HUGEVAL
  rlambdaPmax = -HUGEVAL
  rlambdaSmax = -HUGEVAL

  do ispec=1,nspec

 material = kmato(ispec)

 rlambda = elastcoef(1,material)
 rmu    = elastcoef(2,material)
 rKmod  = elastcoef(3,material)
 denst  = density(material)

 rKvol  = rlambda + 2.d0*rmu/3.d0
 cploc = dsqrt((rKvol + 4.d0*rmu/3.d0)/denst)
 csloc = dsqrt(rmu/denst)

  do j=1,NGLLZ
    do i=1,NGLLX

!--- si formulation heterogene pour un modele de vitesse externe
  if(ireadmodel) then
    ipointnum = ibool(i,j,ispec)
    cploc = vpext(ipointnum)
    csloc = vsext(ipointnum)
    denst = rhoext(ipointnum)
    rmu   = denst*csloc*csloc
    rlambda  = denst*cploc*cploc - 2.d0*rmu
    rKmod  = rlambda + 2.d0*rmu
  endif

!--- calculer min et max du modele de vitesse et densite
  vpmin = dmin1(vpmin,cploc)
  vpmax = dmax1(vpmax,cploc)

  vsmin = dmin1(vsmin,csloc)
  vsmax = dmax1(vsmax,csloc)

  densmin = dmin1(densmin,denst)
  densmax = dmax1(densmax,denst)

!--- stocker parametres pour verifs diverses
  if(i < NGLLX .and. j < NGLLZ) then

    x0 = coord(1,ibool(i,j,ispec))
    z0 = coord(2,ibool(i,j,ispec))
    x1 = coord(1,ibool(i+1,j,ispec))
    z1 = coord(2,ibool(i+1,j,ispec))
    x2 = coord(1,ibool(i,j+1,ispec))
    z2 = coord(2,ibool(i,j+1,ispec))

    rdist1 = dsqrt((x1-x0)**2 + (z1-z0)**2)
    rdist2 = dsqrt((x2-x0)**2 + (z2-z0)**2)

    rsizemin = dmin1(rsizemin,rdist1)
    rsizemin = dmin1(rsizemin,rdist2)
    rsizemax = dmax1(rsizemax,rdist1)
    rsizemax = dmax1(rsizemax,rdist2)

    rapportmin = cploc / dmax1(rdist1,rdist2)
    rapportmax = cploc / dmin1(rdist1,rdist2)
    cpoverdxmin = dmin1(cpoverdxmin,rapportmin)
    cpoverdxmax = dmax1(cpoverdxmax,rapportmax)

    x0 = coord(1,ibool(1,1,ispec))
    z0 = coord(2,ibool(1,1,ispec))
    x1 = coord(1,ibool(NGLLX,1,ispec))
    z1 = coord(2,ibool(NGLLX,1,ispec))
    x2 = coord(1,ibool(1,NGLLZ,ispec))
    z2 = coord(2,ibool(1,NGLLZ,ispec))

    rdist1 = dsqrt((x1-x0)**2 + (z1-z0)**2)
    rdist2 = dsqrt((x2-x0)**2 + (z2-z0)**2)

    rlambmin = cploc/dmax1(rdist1,rdist2)
    rlambmax = cploc/dmin1(rdist1,rdist2)
    rlambdaPmin = dmin1(rlambdaPmin,rlambmin)
    rlambdaPmax = dmax1(rlambdaPmax,rlambmax)

    rlambmin = csloc/dmax1(rdist1,rdist2)
    rlambmax = csloc/dmin1(rdist1,rdist2)
    rlambdaSmin = dmin1(rlambdaSmin,rlambmin)
    rlambdaSmax = dmax1(rlambdaSmax,rlambmax)

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
  if(nint(gltfu(2)) == 2) then

  isourx = nint(gltfu(10))
  isourz = nint(gltfu(11))
  ielems = nint(gltfu(12))

  if(isourx == 1 .or. isourx == NGLLX .or. isourz == 1 .or. isourz == NGLLX) &
        stop 'Explosive source on element edge'

!---- definir a11 et a12 - dirac (schema en croix)

  sig0 = one

  do ir=1,NGLLX
    flagxprime = lagrange_deriv_GLL(ir-1,isourx-1,xigll,NGLLX)
    a11(ir,isourz) = a11(ir,isourz) + sig0*xix(isourx,isourz,ielems)*flagxprime
    a12(ir,isourz) = a12(ir,isourz) + sig0*xiz(isourx,isourz,ielems)*flagxprime
  enddo

  do is=1,NGLLZ
    flagzprime = lagrange_deriv_GLL(is-1,isourz-1,zigll,NGLLZ)
    a11(isourx,is) = a11(isourx,is) + sig0*gammax(isourx,isourz,ielems)*flagzprime
    a12(isourx,is) = a12(isourx,is) + sig0*gammaz(isourx,isourz,ielems)*flagzprime
  enddo

  endif

  end subroutine defarrays

