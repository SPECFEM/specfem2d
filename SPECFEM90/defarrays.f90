
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

  subroutine defarrays(vpext,vsext,rhoext,density,elastcoef, &
                      xi,yi,wx,wy,hprime,hTprime, &
                      a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13x,a13z, &
                      ibool,iboolori,kmato,dvolu,xjaci,coord,gltfu, &
                      numabs,codeabs,anyabs,anyperio)
!
!=======================================================================
!
!     "d e f a r r a y s" : Define arrays a1 to a13 for the spectral
!                           elements solver
!
!=======================================================================
!

  use loadft
  use iounit
  use infos
  use mesh01
  use spela202
  use constspec
  use vparams
  use verifs
  use energie
  use codebord

  implicit none

  integer kmato(nspec),ibool(0:nxgll-1,0:nxgll-1,nspec)
  integer iboolori(0:nxgll-1,0:nxgll-1,nspec)

  double precision density(numat),elastcoef(4,numat), &
              xi(0:nxgll-1),yi(0:nygll-1),wx(0:nxgll-1),wy(0:nxgll-1), &
              dvolu(nspec,0:nxgll-1,0:nxgll-1), &
              xjaci(nspec,ndime,ndime,0:nxgll-1,0:nxgll-1), &
              hprime(0:nxgll-1,0:nxgll-1), hTprime(0:nxgll-1,0:nxgll-1)

  double precision coord(ndime,npoin)
  double precision a1(0:nxgll-1,0:nxgll-1,nspec), &
     a2(0:nxgll-1,0:nxgll-1,nspec), &
     a3(0:nxgll-1,0:nxgll-1,nspec),a4(0:nxgll-1,0:nxgll-1,nspec), &
     a5(0:nxgll-1,0:nxgll-1,nspec),a6(0:nxgll-1,0:nxgll-1,nspec), &
     a7(0:nxgll-1,0:nxgll-1,nspec),a8(0:nxgll-1,0:nxgll-1,nspec), &
     a9(0:nxgll-1,0:nxgll-1,nspec),a10(0:nxgll-1,0:nxgll-1,nspec)
  double precision a13x(0:nxgll-1,0:nxgll-1,nelemabs), &
     a13z(0:nxgll-1,0:nxgll-1,nelemabs)
  double precision a11(0:nxgll-1,0:nxgll-1,nltfl), &
     a12(0:nxgll-1,0:nxgll-1,nltfl)

  double precision gltfu(20,nltfl)
  double precision vpext(npoin)
  double precision vsext(npoin)
  double precision rhoext(npoin)

  integer numabs(nelemabs),codeabs(4,nelemabs)

  double precision, external :: hdgll

  double precision, parameter :: zero=0.d0,one=1.d0

  integer i,j,k
  integer numelem,material
  integer ipointnum,n
  integer isourx,isourz,ielems,ir,is,ip,noffsetelem
  double precision vsmin,vsmax,densmin,densmax
  double precision rKmod,rlamda,rmu,xix,xiz,etax,etaz,denst,rjacob
  double precision rKvol,cploc,csloc,xxi,zeta,rwx,x0,z0
  double precision c11,c13,c33,c44
  double precision x1,z1,x2,z2,rdist1,rdist2,rapportmin,rapportmax
  double precision rlambmin,rlambmax,coefintegr
  double precision flagxprime,flagzprime,sig0
  logical anyabs,anyperio,anisotrope

!
!-----------------------------------------------------------------------
!

!---- compute parameters a1 to a13 for the spectral elements

  a11 = zero
  a12 = zero

  a13x = zero
  a13z = zero

  vpmin = 1.d30
  vsmin = 1.d30
  vpmax = -1.d30
  vsmax = -1.d30
  densmin = 1.d30
  densmax = -1.d30

  rsizemin = 1.d30
  rsizemax = -1.d30

  cpoverdxmin = 1.d30
  cpoverdxmax = -1.d30

  rlamdaPmin = 1.d30
  rlamdaSmin = 1.d30
  rlamdaPmax = -1.d30
  rlamdaSmax = -1.d30

  do numelem=1,nspec

 material = kmato(numelem)

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

  do i=0,nxgll-1
    do j=0,nygll-1

  xix = xjaci(numelem,1,1,i,j)
  xiz = xjaci(numelem,1,2,i,j)
  etax = xjaci(numelem,2,1,i,j)
  etaz = xjaci(numelem,2,2,i,j)
  rjacob = dvolu(numelem,i,j)

  xxi = etaz * rjacob
  zeta = xix * rjacob

  rwx = - wx(i)*wy(j)

!--- si formulation heterogene pour un modele de vitesse externe
  if(ireadmodel) then
    ipointnum = ibool(i,j,numelem)
    cploc = vpext(ipointnum)
    csloc = vsext(ipointnum)
    denst = rhoext(ipointnum)
    rmu   = denst*csloc*csloc
    rlamda  = denst*cploc*cploc - 2.d0*rmu
    rKmod  = rlamda + 2.d0*rmu
  endif

!--- si materiau transverse isotrope, donner une idee des proprietes
  if (anisotrope) then
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
  if(i /= nxgll-1 .and. j /= nygll-1) then

  if(anyperio) then
    x0 = coord(1,iboolori(i,j,numelem))
    z0 = coord(2,iboolori(i,j,numelem))
    x1 = coord(1,iboolori(i+1,j,numelem))
    z1 = coord(2,iboolori(i+1,j,numelem))
    x2 = coord(1,iboolori(i,j+1,numelem))
    z2 = coord(2,iboolori(i,j+1,numelem))
  else
    x0 = coord(1,ibool(i,j,numelem))
    z0 = coord(2,ibool(i,j,numelem))
    x1 = coord(1,ibool(i+1,j,numelem))
    z1 = coord(2,ibool(i+1,j,numelem))
    x2 = coord(1,ibool(i,j+1,numelem))
    z2 = coord(2,ibool(i,j+1,numelem))
  endif

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

  if(anyperio) then
    x0 = coord(1,iboolori(0,0,numelem))
    z0 = coord(2,iboolori(0,0,numelem))
    x1 = coord(1,iboolori(nxgll-1,0,numelem))
    z1 = coord(2,iboolori(nxgll-1,0,numelem))
    x2 = coord(1,iboolori(0,nygll-1,numelem))
    z2 = coord(2,iboolori(0,nygll-1,numelem))
  else
    x0 = coord(1,ibool(0,0,numelem))
    z0 = coord(2,ibool(0,0,numelem))
    x1 = coord(1,ibool(nxgll-1,0,numelem))
    z1 = coord(2,ibool(nxgll-1,0,numelem))
    x2 = coord(1,ibool(0,nygll-1,numelem))
    z2 = coord(2,ibool(0,nygll-1,numelem))
  endif

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
    a1(i,j,numelem) = rwx*(rKmod*xix*xix + rmu*xiz*xiz)*rjacob
    a2(i,j,numelem) = rwx*(rKmod*etax*xix + rmu*etaz*xiz)*rjacob
    a3(i,j,numelem) = rwx*(rlamda+rmu)*xiz*xix*rjacob
    a4(i,j,numelem) = rwx*(rlamda*etaz*xix + rmu*etax*xiz)*rjacob
    a5(i,j,numelem) = rwx*(rKmod*etaz*etaz + rmu*etax*etax)*rjacob
    a6(i,j,numelem) = rwx*(rKmod*etax*etax + rmu*etaz*etaz)*rjacob
    a7(i,j,numelem) = rwx*(rlamda*etax*xiz + rmu*etaz*xix)*rjacob
    a8(i,j,numelem) = rwx*(rlamda+rmu)*etax*etaz*rjacob
    a9(i,j,numelem) = rwx*(rKmod*xiz*xiz + rmu*xix*xix)*rjacob
    a10(i,j,numelem) = rwx*(rKmod*etaz*xiz + rmu*etax*xix)*rjacob
  else
    a3(i,j,numelem) = rwx*(c13+c44)*xiz*xix*rjacob
    a4(i,j,numelem) = rwx*(c13*etaz*xix + c44*etax*xiz)*rjacob
    a7(i,j,numelem) = rwx*(c13*etax*xiz + c44*etaz*xix)*rjacob
    a8(i,j,numelem) = rwx*(c13+c44)*etax*etaz*rjacob

    a1(i,j,numelem) = rwx*(c11*xix*xix + c44*xiz*xiz)*rjacob
    a2(i,j,numelem) = rwx*(c11*etax*xix + c44*etaz*xiz)*rjacob
    a6(i,j,numelem) = rwx*(c11*etax*etax + c44*etaz*etaz)*rjacob

    a5(i,j,numelem) = rwx*(c33*etaz*etaz + c44*etax*etax)*rjacob
    a9(i,j,numelem) = rwx*(c33*xiz*xiz + c44*xix*xix)*rjacob
    a10(i,j,numelem) = rwx*(c33*etaz*xiz + c44*etax*xix)*rjacob
  endif

!--- valeurs pour solution analytique (recuperer deux points de topo)
  noffsetelem = 20
  if(numelem  ==  nspec-noffsetelem.and.i  ==  0.and.j  ==  nygll-1) then
    cp1 = cploc
    cs1 = csloc
    rho1 = denst
  if(anyperio) then
    xt1 = coord(1,iboolori(i,j,numelem))
    zt1 = coord(2,iboolori(i,j,numelem))
  else
    xt1 = coord(1,ibool(i,j,numelem))
    zt1 = coord(2,ibool(i,j,numelem))
  endif
  else if(numelem  ==  nspec.and.i  ==  nxgll-1.and. j  ==  nygll-1) then
 if(anyperio) then
    xt2 = coord(1,iboolori(i,j,numelem))
    zt2 = coord(2,iboolori(i,j,numelem))
  else
    xt2 = coord(1,ibool(i,j,numelem))
    zt2 = coord(2,ibool(i,j,numelem))
  endif
  else if(numelem  ==  1) then
    cp2 = cploc
    cs2 = csloc
    rho2 = denst
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

  do numelem=1,nelemabs

 material = kmato(numabs(numelem))

 rlamda = elastcoef(1,material)
 rmu    = elastcoef(2,material)
 rKmod  = elastcoef(3,material)
 denst  = density(material)

 rKvol  = rlamda + 2.d0*rmu/3.d0
 cploc = dsqrt((rKvol + 4.d0*rmu/3.d0)/denst)
 csloc = dsqrt(rmu/denst)

  do i=0,nxgll-1
    do j=0,nygll-1

!--- si formulation heterogene pour un modele de vitesse externe
  if(ireadmodel) then
     ipointnum = ibool(i,j,numabs(numelem))
     cploc = vpext(ipointnum)
     csloc = vsext(ipointnum)
     denst = rhoext(ipointnum)
     rmu   = denst*csloc*csloc
     rlamda  = denst*cploc*cploc - 2.d0*rmu
     rKmod  = rlamda + 2.d0*rmu
  endif

  xix = xjaci(numabs(numelem),1,1,i,j)
  xiz = xjaci(numabs(numelem),1,2,i,j)
  etax = xjaci(numabs(numelem),2,1,i,j)
  etaz = xjaci(numabs(numelem),2,2,i,j)
  rjacob = dvolu(numabs(numelem),i,j)

  xxi = etaz * rjacob
  zeta = xix * rjacob

  rwx = - wx(i)*wy(j)

!---- sommer les contributions dans les coins pour l'ancienne formulation
!---- ne pas sommer les contributions dans les coins pour la nouvelle

! bord absorbant du bas
  if(codeabs(ibas,numelem) /= 0 .and. j == 0) then
    coefintegr = wx(i)*xxi
    a13x(i,j,numelem) = denst*csloc*coefintegr
    a13z(i,j,numelem) = denst*cploc*coefintegr
  endif

! bord absorbant du haut (signe moins)
  if(codeabs(ihaut,numelem) /= 0 .and. j == nygll-1) then
    coefintegr = wx(i)*xxi
    a13x(i,j,numelem) = denst*csloc*coefintegr
    a13z(i,j,numelem) = denst*cploc*coefintegr
  endif

! bord absorbant de gauche
  if(codeabs(igauche,numelem) /= 0 .and. i == 0) then
    coefintegr = wy(j)*zeta
    a13x(i,j,numelem) = denst*cploc*coefintegr
    a13z(i,j,numelem) = denst*csloc*coefintegr
  endif

! bord absorbant de droite
  if(codeabs(idroite,numelem) /= 0 .and. i == nxgll-1) then
    coefintegr = wy(j)*zeta
    a13x(i,j,numelem) = denst*cploc*coefintegr
    a13z(i,j,numelem) = denst*csloc*coefintegr
  endif

    enddo
 enddo
  enddo

  endif

! pour source explosive
  do n=1,nltfl

! seulement si source explosive
  if(nint(gltfu(2,n)) == 2) then

  isourx = nint(gltfu(10,n))
  isourz = nint(gltfu(11,n))
  ielems = nint(gltfu(12,n))

  if(isourx == 0.or.isourx == nxgll-1.or.isourz == 0 .or.isourz == nxgll-1) &
        stop 'Explosive source on element edge'

!---- definir a11 et a12 - dirac (schema en croix)

  xix = xjaci(ielems,1,1,isourx,isourz)
  xiz = xjaci(ielems,1,2,isourx,isourz)
  etax = xjaci(ielems,2,1,isourx,isourz)
  etaz = xjaci(ielems,2,2,isourx,isourz)

  sig0 = one

  do ir=0,nxgll-1
    flagxprime = hdgll(ir,isourx,xi,nxgll)
    a11(ir,isourz,n) = a11(ir,isourz,n) + sig0*xix*flagxprime
    a12(ir,isourz,n) = a12(ir,isourz,n) + sig0*xiz*flagxprime
  enddo

  do is=0,nygll-1
    flagzprime = hdgll(is,isourz,yi,nygll)
    a11(isourx,is,n) = a11(isourx,is,n) + sig0*etax*flagzprime
    a12(isourx,is,n) = a12(isourx,is,n) + sig0*etaz*flagzprime
  enddo

  endif

  enddo

!---- compute hprime coefficients (derivatives of Lagrange polynomials)
!---- (works only if nxgll = nygll)
  do ip=0,nxgll-1
    do i=0,nxgll-1
      hprime(ip,i) = hdgll(ip,i,xi,nxgll)
      hTprime(i,ip) = hprime(ip,i)
    enddo
  enddo

  return
  end subroutine defarrays

