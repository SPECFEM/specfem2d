
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

  subroutine qinpspec(density,elastcoef,xi,yi,wx,wy,knods, &
                      ibool,kmato,shape,shapeint,dershape,dvolu,xjaci, &
                      coorg,xirec,etarec,flagrange, &
          numabs,codeabs,codeperio,anyabs,anyperio)
!
!=======================================================================
!
!     "q i n p s p e c" : Read, generate and write data for the spectral
!                         elements
!
!=======================================================================
!

  use iounit
  use infos
  use mesh01
  use spela202

  implicit none

  double precision, parameter :: zero=0.d0
  double precision, parameter :: gaussalpha=zero,gaussbeta=zero

! choix entre version lente et rapide pour la numerotation globale
  logical, parameter :: fast_numbering = .false.

  integer knods(ngnod,nspec),kmato(nspec),ibool(0:nxgll-1,0:nxgll-1,nspec)

  double precision density(numat),elastcoef(4,numat), &
              xi(0:nxgll-1),yi(0:nygll-1),wx(0:nxgll-1),wy(0:nxgll-1), &
              dvolu(nspec,nxgll,nxgll),xjaci(nspec,ndime,ndime,nxgll,nxgll), &
              coorg(ndime,npgeo)
  double precision shape(ngnod,nxgll,nxgll)
  double precision shapeint(ngnod,iptsdisp,iptsdisp)
  double precision dershape(ndime,ngnod,nxgll,nxgll)
  double precision xirec(iptsdisp),etarec(iptsdisp)
  double precision flagrange(0:nxgll-1,iptsdisp)

  integer numabs(nelemabs),codeabs(4,nelemabs)
  integer codeperio(4,nelemperio)
  logical anyabs,anyperio

  integer nelemabs2,nelemperio2

!
!-----------------------------------------------------------------------
!

! check that numbering is fine (no fast numbering if periodic conditions)
  if(fast_numbering .and. anyperio) stop 'no fast numbering if periodic conditions'

!
!---- print element group main parameters
!
  nelemabs2 = nelemabs
  nelemperio2 = nelemperio
  if(.not. anyabs) nelemabs2 = 0
  if(.not. anyperio) nelemperio2 = 0
  if(iecho /= 0) then
    write(iout,100)
    write(iout,200) nspec,ngnod,nxgll, &
           nygll,nxgll*nygll,iptsdisp,numat,nelemabs2,nelemperio2
  endif

!
!----    set up coordinates of the Gauss-Lobatto-Legendre points
!
 call zwgljd(xi,wx,nxgll,gaussalpha,gaussbeta)
 call zwgljd(yi,wy,nygll,gaussalpha,gaussbeta)

!
!---- if nb of points is odd, the middle abscissa is exactly zero
!
  if(mod(nxgll,2)  /=  0) xi((nxgll-1)/2) = zero
  if(mod(nygll,2)  /=  0) yi((nygll-1)/2) = zero

!
!---- read the material properties
!
  call gmat01(density,elastcoef,numat)

!
!---- read topology and material number for spectral elements
!
  call getelspec(knods,kmato,numabs,codeabs,codeperio,anyabs,anyperio)

!
!---- compute the spectral element shape functions and their local derivatives
!
  call q49shape(shape,dershape,xi,yi,ngnod,nxgll,nygll,ndime)

!
!---- generate the global numbering
!

! version "propre mais lente" ou version "sale mais rapide"
  if(fast_numbering) then
      call createnum_fast(knods,ibool,shape,coorg,npoin,ndime,npgeo)
  else
      call createnum_slow(knods,ibool,npoin)
  endif

!
!---- compute the spectral element jacobian matrix
!

  call q49spec(shapeint,dershape,dvolu,xjaci,xi,coorg, &
                 knods,ngnod,nxgll,nygll,ndime,nspec,npgeo, &
                 xirec,etarec,flagrange,iptsdisp)

!
!---- formats
!
  100   format(/5x,'--> Isoparametric Spectral Elements <--',//)
  200   format(5x, &
           'Number of spectral elements . . . . .  (nspec) =',i7,/5x, &
           'Number of control nodes per element .  (ngnod) =',i7,/5x, &
           'Number of points in X-direction . . .  (nxgll) =',i7,/5x, &
           'Number of points in Y-direction . . .  (nygll) =',i7,/5x, &
           'Number of points per element. . .(nxgll*nygll) =',i7,/5x, &
           'Number of points for display . . . .(iptsdisp) =',i7,/5x, &
           'Number of element material sets . . .  (numat) =',i7,/5x, &
           'Number of absorbing elements . . . .(nelemabs) =',i7,/5x, &
           'Number of periodic elements. . . .(nelemperio) =',i7)

  end subroutine qinpspec

