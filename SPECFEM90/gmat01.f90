
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

  subroutine gmat01(density,elastcoef,numat)

! read properties of a 2D isotropic or anisotropic linear elastic element

  implicit none

  include "constants.h"

  character(len=80) datlin
  double precision Kmod,Kvol

  integer numat
  double precision density(numat),elastcoef(4,numat)

  integer in,n,indic
  double precision young,poisson,denst,cp,cs,amu,a2mu,alam
  double precision val1,val2,val3,val4
  double precision c11,c13,c33,c44

!
!---- loop over the different material sets
!
  density(:) = zero
  elastcoef(:,:) = zero

  write(iout,100) numat

  read(iin ,40) datlin
  do in = 1,numat

   read(iin ,*) n,indic,denst,val1,val2,val3,val4

   if(n<1 .or. n>numat) stop 'Wrong material set number'

!---- materiau isotrope, vitesse P et vitesse S donnees
   if(indic == 0) then
      cp = val1
      cs = val2
      amu   = denst*cs*cs
      a2mu  = 2.d0*amu
      alam  = denst*cp*cp - a2mu
      Kmod  = alam + a2mu
      Kvol  = alam + a2mu/3.d0
      young = 9.d0*Kvol*amu/(3.d0*Kvol + amu)
      poisson = half*(3.d0*Kvol-a2mu)/(3.d0*Kvol+amu)
! Poisson's ratio must be between -1 and +1/2
      if (poisson < -1.d0 .or. poisson > 0.5d0) stop 'Poisson''s ratio out of range'

!---- materiau isotrope, module de Young et coefficient de Poisson donnes
   else if(indic == 1) then
      young = val1
      poisson = val2
! Poisson's ratio must be between -1 and +1/2
      if (poisson < -1.d0 .or. poisson > 0.5d0) stop 'Poisson''s ratio out of range'
      a2mu  = young/(one+poisson)
      amu   = half*a2mu
      alam  = a2mu*poisson/(one-two*poisson)
      Kmod  = alam + a2mu
      Kvol  = alam + a2mu/3.d0
      cp    = dsqrt((Kvol + 4.d0*amu/3.d0)/denst)
      cs    = dsqrt(amu/denst)

!---- materiau anisotrope, c11, c13, c33 et c44 donnes en Pascal
   else if(indic == 2) then
      c11 = val1
      c13 = val2
      c33 = val3
      c44 = val4

   else
      stop 'Improper value while reading material sets'
   endif

!
!----  set elastic coefficients and density
!
!  Isotropic              :  lambda, mu, K (= lambda + 2*mu), zero
!  Transverse anisotropic :  c11, c13, c33, c44
!
  if(indic == 0 .or. indic == 1) then
    elastcoef(1,n) = alam
    elastcoef(2,n) = amu
    elastcoef(3,n) = Kmod
    elastcoef(4,n) = zero
  else
    elastcoef(1,n) = c11
    elastcoef(2,n) = c13
    elastcoef(3,n) = c33
    elastcoef(4,n) = c44
  endif

  density(n) = denst

!
!----    check the input
!
  if(indic == 0 .or. indic == 1) then
    write(iout,200) n,cp,cs,denst,poisson,alam,amu,Kvol,young
  else
    write(iout,300) n,c11,c13,c33,c44,denst, &
        sqrt(c33/denst),sqrt(c11/denst),sqrt(c44/denst),sqrt(c44/denst)
  endif

  enddo

!
!---- formats
!
  40    format(a80)
  100   format(//,' M a t e r i a l   s e t s :  ', &
         ' 2 D  e l a s t i c i t y', &
         /1x,54('='),//5x,'Number of material sets . . . . . . (numat) =',i5)
  200   format(//5x,'------------------------',/5x, &
         '-- Isotropic material --',/5x, &
         '------------------------',/5x, &
         'Material set number. . . . . . . . (jmat) =',i5,/5x, &
         'P-wave velocity. . . . . . . . . . . (cp) =',1pe15.8,/5x, &
         'S-wave velocity. . . . . . . . . . . (cs) =',1pe15.8,/5x, &
         'Mass density. . . . . . . . . . . (denst) =',1pe15.8,/5x, &
         'Poisson''s ratio . . . . . . . .(poisson) =',1pe15.8,/5x, &
         'First Lame parameter Lambda. . . . (alam) =',1pe15.8,/5x, &
         'Second Lame parameter Mu. . . . . . (amu) =',1pe15.8,/5x, &
         'Bulk modulus K . . . . . . . . . . (Kvol) =',1pe15.8,/5x, &
         'Young''s modulus E . . . . . . . . (young) =',1pe15.8)
  300   format(//5x,'-------------------------------------',/5x, &
         '-- Transverse anisotropic material --',/5x, &
         '-------------------------------------',/5x, &
         'Material set number. . . . . . . . (jmat) =',i5,/5x, &
         'c11 coefficient (Pascal). . . . . . (c11) =',1pe15.8,/5x, &
         'c13 coefficient (Pascal). . . . . . (c13) =',1pe15.8,/5x, &
         'c33 coefficient (Pascal). . . . . . (c33) =',1pe15.8,/5x, &
         'c44 coefficient (Pascal). . . . . . (c44) =',1pe15.8,/5x, &
         'Mass density. . . . . . . . . . . (denst) =',1pe15.8,/5x, &
         'Velocity of qP along vertical axis. . . . =',1pe15.8,/5x, &
         'Velocity of qP along horizontal axis. . . =',1pe15.8,/5x, &
         'Velocity of qSV along vertical axis . . . =',1pe15.8,/5x, &
         'Velocity of qSV along horizontal axis . . =',1pe15.8)

  end subroutine gmat01

