
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
          xigll,yigll,wxgll,wygll,hprime,hTprime, &
          a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13x,a13z, &
          ibool,kmato,dvolu,xjaci,coord,gltfu, &
          numabs,codeabs,anyabs,npoin,rsizemin,rsizemax, &
          cpoverdxmin,cpoverdxmax,rlamdaSmin,rlamdaSmax, &
          rlamdaPmin,rlamdaPmax,vpmin,vpmax,ireadmodel,nelemabs,nspec,numat)


! define all the arrays for the variational formulation

  implicit none

  include "constants.h"

  integer i,j,ispec,material,ipointnum,npoin,nelemabs,nspec,numat
  integer isourx,isourz,ielems,ir,is,ip

  integer kmato(nspec),ibool(NGLLX,NGLLX,nspec)

  double precision density(numat),elastcoef(4,numat), &
              xigll(NGLLX),yigll(NGLLY),wxgll(NGLLX),wygll(NGLLX), &
              dvolu(nspec,NGLLX,NGLLX), &
              xjaci(nspec,NDIME,NDIME,NGLLX,NGLLX), &
              hprime(NGLLX,NGLLX), hTprime(NGLLX,NGLLX)

  double precision coord(NDIME,npoin)
  double precision a1(NGLLX,NGLLX,nspec),a2(NGLLX,NGLLX,nspec), &
     a3(NGLLX,NGLLX,nspec),a4(NGLLX,NGLLX,nspec), &
     a5(NGLLX,NGLLX,nspec),a6(NGLLX,NGLLX,nspec), &
     a7(NGLLX,NGLLX,nspec),a8(NGLLX,NGLLX,nspec), &
     a9(NGLLX,NGLLX,nspec),a10(NGLLX,NGLLX,nspec)
  double precision a13x(NGLLX,NGLLX,nelemabs),a13z(NGLLX,NGLLX,nelemabs)
  double precision a11(NGLLX,NGLLX),a12(NGLLX,NGLLX)

  double precision gltfu(20)
  double precision vpext(npoin)
  double precision vsext(npoin)
  double precision rhoext(npoin)

  double precision vsmin,vsmax,densmin,densmax
  double precision rKmod,rlamda,rmu,xix,xiz,etax,etaz,denst,rjacob
  double precision rKvol,cploc,csloc,xxi,zeta,rwgll,x0,z0
  double precision c11,c13,c33,c44
  double precision x1,z1,x2,z2,rdist1,rdist2,rapportmin,rapportmax
  double precision rlambmin,rlambmax,coefintegr
  double precision flagxprime,flagzprime,sig0

  double precision rsizemin,rsizemax,cpoverdxmin,cpoverdxmax, &
    rlamdaSmin,rlamdaSmax,rlamdaPmin,rlamdaPmax,vpmin,vpmax

  logical anyabs,anisotrope,ireadmodel

  integer numabs(nelemabs),codeabs(4,nelemabs)

  double precision, external :: lagrange_deriv_GLL

!
!-----------------------------------------------------------------------
!

!---- compute parameters a1 to a13 for the spectral elements

  a11 = zero
  a12 = zero

  a13x = zero
  a13z = zero

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

  rlamdaPmin = HUGEVAL
  rlamdaSmin = HUGEVAL
  rlamdaPmax = -HUGEVAL
  rlamdaSmax = -HUGEVAL

  do ispec=1,nspec


 material = kmato(ispec)

 rlamda = elastcoef(1,material)
 rmu    = elastcoef(2,material)
 rKmod  = elastcoef(3,material)
 denst  = density(material)

 rKvol  = rlamda + 2.d0*rmu/3.d0
 cploc = dsqrt((rKvol + 4.d0*rmu/3.d0)/denst)
 csloc = dsqrt(rmu/denst)

! determiner si le materiau est anisotrope ou non
  if(elastcoef(4,material)  >  0.00001d0) then
    anisotrope = .true.
    c11 = elastcoef(1,material)
    c13 = elastcoef(2,material)
    c33 = elastcoef(3,material)
    c44 = elastcoef(4,material)
  else
    anisotrope = .false.
  endif

  do j=1,NGLLY
    do i=1,NGLLX

  xix = xjaci(ispec,1,1,i,j)
  xiz = xjaci(ispec,1,2,i,j)
  etax = xjaci(ispec,2,1,i,j)
  etaz = xjaci(ispec,2,2,i,j)
  rjacob = dvolu(ispec,i,j)

  xxi = etaz * rjacob
  zeta = xix * rjacob

  rwgll = - wxgll(i)*wygll(j)

!--- si formulation heterogene pour un modele de vitesse externe
  if(ireadmodel) then
    ipointnum = ibool(i,j,ispec)
    cploc = vpext(ipointnum)
    csloc = vsext(ipointnum)
    denst = rhoext(ipointnum)
    rmu   = denst*csloc*csloc
    rlamda  = denst*cploc*cploc - 2.d0*rmu
    rKmod  = rlamda + 2.d0*rmu
  endif

!--- si materiau transverse isotrope, donner une idee des proprietes
  if(anisotrope) then
    cploc = sqrt(c11/denst)
    csloc = sqrt(c44/denst)
  endif

!--- calculer min et max du modele de vitesse et densite
  vpmin = dmin1(vpmin,cploc)
  vpmax = dmax1(vpmax,cploc)

  vsmin = dmin1(vsmin,csloc)
  vsmax = dmax1(vsmax,csloc)

  densmin = dmin1(densmin,denst)
  densmax = dmax1(densmax,denst)

!--- stocker parametres pour verifs diverses
  if(i < NGLLX .and. j < NGLLY) then

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
    x2 = coord(1,ibool(1,NGLLY,ispec))
    z2 = coord(2,ibool(1,NGLLY,ispec))

    rdist1 = dsqrt((x1-x0)**2 + (z1-z0)**2)
    rdist2 = dsqrt((x2-x0)**2 + (z2-z0)**2)

    rlambmin = cploc/dmax1(rdist1,rdist2)
    rlambmax = cploc/dmin1(rdist1,rdist2)
    rlamdaPmin = dmin1(rlamdaPmin,rlambmin)
    rlamdaPmax = dmax1(rlamdaPmax,rlambmax)

    rlambmin = csloc/dmax1(rdist1,rdist2)
    rlambmax = csloc/dmin1(rdist1,rdist2)
    rlamdaSmin = dmin1(rlamdaSmin,rlambmin)
    rlamdaSmax = dmax1(rlamdaSmax,rlambmax)

  endif

!--- definir tableaux
  if(.not. anisotrope) then
    a1(i,j,ispec) = rwgll*(rKmod*xix*xix + rmu*xiz*xiz)*rjacob
    a2(i,j,ispec) = rwgll*(rKmod*etax*xix + rmu*etaz*xiz)*rjacob
    a3(i,j,ispec) = rwgll*(rlamda+rmu)*xiz*xix*rjacob
    a4(i,j,ispec) = rwgll*(rlamda*etaz*xix + rmu*etax*xiz)*rjacob
    a5(i,j,ispec) = rwgll*(rKmod*etaz*etaz + rmu*etax*etax)*rjacob
    a6(i,j,ispec) = rwgll*(rKmod*etax*etax + rmu*etaz*etaz)*rjacob
    a7(i,j,ispec) = rwgll*(rlamda*etax*xiz + rmu*etaz*xix)*rjacob
    a8(i,j,ispec) = rwgll*(rlamda+rmu)*etax*etaz*rjacob
    a9(i,j,ispec) = rwgll*(rKmod*xiz*xiz + rmu*xix*xix)*rjacob
    a10(i,j,ispec) = rwgll*(rKmod*etaz*xiz + rmu*etax*xix)*rjacob
  else
    a3(i,j,ispec) = rwgll*(c13+c44)*xiz*xix*rjacob
    a4(i,j,ispec) = rwgll*(c13*etaz*xix + c44*etax*xiz)*rjacob
    a7(i,j,ispec) = rwgll*(c13*etax*xiz + c44*etaz*xix)*rjacob
    a8(i,j,ispec) = rwgll*(c13+c44)*etax*etaz*rjacob

    a1(i,j,ispec) = rwgll*(c11*xix*xix + c44*xiz*xiz)*rjacob
    a2(i,j,ispec) = rwgll*(c11*etax*xix + c44*etaz*xiz)*rjacob
    a6(i,j,ispec) = rwgll*(c11*etax*etax + c44*etaz*etaz)*rjacob

    a5(i,j,ispec) = rwgll*(c33*etaz*etaz + c44*etax*etax)*rjacob
    a9(i,j,ispec) = rwgll*(c33*xiz*xiz + c44*xix*xix)*rjacob
    a10(i,j,ispec) = rwgll*(c33*etaz*xiz + c44*etax*xix)*rjacob
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

!
!--- definition coefficients pour bords absorbants
!

  if(anyabs) then

  do ispec=1,nelemabs

 material = kmato(numabs(ispec))

 rlamda = elastcoef(1,material)
 rmu    = elastcoef(2,material)
 rKmod  = elastcoef(3,material)
 denst  = density(material)

 rKvol  = rlamda + 2.d0*rmu/3.d0
 cploc = dsqrt((rKvol + 4.d0*rmu/3.d0)/denst)
 csloc = dsqrt(rmu/denst)

  do i=1,NGLLX
    do j=1,NGLLY

!--- si formulation heterogene pour un modele de vitesse externe
  if(ireadmodel) then
     ipointnum = ibool(i,j,numabs(ispec))
     cploc = vpext(ipointnum)
     csloc = vsext(ipointnum)
     denst = rhoext(ipointnum)
     rmu   = denst*csloc*csloc
     rlamda  = denst*cploc*cploc - 2.d0*rmu
     rKmod  = rlamda + 2.d0*rmu
  endif

  xix = xjaci(numabs(ispec),1,1,i,j)
  xiz = xjaci(numabs(ispec),1,2,i,j)
  etax = xjaci(numabs(ispec),2,1,i,j)
  etaz = xjaci(numabs(ispec),2,2,i,j)
  rjacob = dvolu(numabs(ispec),i,j)

  xxi = etaz * rjacob
  zeta = xix * rjacob

  rwgll = - wxgll(i)*wygll(j)

!---- sommer les contributions dans les coins pour l'ancienne formulation
!---- ne pas sommer les contributions dans les coins pour la nouvelle

! bord absorbant du bas
  if(codeabs(ibas,ispec) /= 1 .and. j == 1) then
    coefintegr = wxgll(i)*xxi
    a13x(i,j,ispec) = denst*csloc*coefintegr
    a13z(i,j,ispec) = denst*cploc*coefintegr
  endif

! bord absorbant du haut (signe moins)
  if(codeabs(ihaut,ispec) /= 1 .and. j == NGLLY) then
    coefintegr = wxgll(i)*xxi
    a13x(i,j,ispec) = denst*csloc*coefintegr
    a13z(i,j,ispec) = denst*cploc*coefintegr
  endif

! bord absorbant de gauche
  if(codeabs(igauche,ispec) /= 1 .and. i == 1) then
    coefintegr = wygll(j)*zeta
    a13x(i,j,ispec) = denst*cploc*coefintegr
    a13z(i,j,ispec) = denst*csloc*coefintegr
  endif

! bord absorbant de droite
  if(codeabs(idroite,ispec) /= 1 .and. i == NGLLX) then
    coefintegr = wygll(j)*zeta
    a13x(i,j,ispec) = denst*cploc*coefintegr
    a13z(i,j,ispec) = denst*csloc*coefintegr
  endif

    enddo
 enddo
  enddo

  endif

! seulement si source explosive
  if(nint(gltfu(2)) == 2) then

  isourx = nint(gltfu(10))
  isourz = nint(gltfu(11))
  ielems = nint(gltfu(12))

  if(isourx == 1.or.isourx == NGLLX.or.isourz == 1 .or.isourz == NGLLX) &
        stop 'Explosive source on element edge'

!---- definir a11 et a12 - dirac (schema en croix)

  xix = xjaci(ielems,1,1,isourx,isourz)
  xiz = xjaci(ielems,1,2,isourx,isourz)
  etax = xjaci(ielems,2,1,isourx,isourz)
  etaz = xjaci(ielems,2,2,isourx,isourz)

  sig0 = one

  do ir=1,NGLLX
    flagxprime = lagrange_deriv_GLL(ir-1,isourx-1,xigll,NGLLX)
    a11(ir,isourz) = a11(ir,isourz) + sig0*xix*flagxprime
    a12(ir,isourz) = a12(ir,isourz) + sig0*xiz*flagxprime
  enddo

  do is=1,NGLLY
    flagzprime = lagrange_deriv_GLL(is-1,isourz-1,yigll,NGLLY)
    a11(isourx,is) = a11(isourx,is) + sig0*etax*flagzprime
    a12(isourx,is) = a12(isourx,is) + sig0*etaz*flagzprime
  enddo

  endif

!---- compute hprime coefficients (derivatives of Lagrange polynomials)
!---- (works only if NGLLX = NGLLY)
  do ip=1,NGLLX
    do i=1,NGLLX
      hprime(ip,i) = lagrange_deriv_GLL(ip-1,i-1,xigll,NGLLX)
      hTprime(i,ip) = hprime(ip,i)
    enddo
  enddo

  end subroutine defarrays

