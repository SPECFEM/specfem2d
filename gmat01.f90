
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
! Copyright Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

  subroutine gmat01(density_array,elastcoef,numat)

! read properties of a 2D isotropic or anisotropic linear elastic element

  implicit none

  include "constants.h"

  character(len=80) datlin
  double precision lambdaplus2mu,kappa

  integer numat
  double precision density_array(numat),elastcoef(4,numat)

  integer in,n,indic
  double precision young,poisson,density,cp,cs,mu,two_mu,lambda
  double precision val1,val2,val3,val4
  double precision c11,c13,c33,c44

!
!---- loop over the different material sets
!
  density_array(:) = zero
  elastcoef(:,:) = zero

  write(iout,100) numat

  read(iin ,"(a80)") datlin
  do in = 1,numat

   read(iin ,*) n,indic,density,val1,val2,val3,val4

   if(n<1 .or. n>numat) call exit_MPI('Wrong material set number')

!---- isotropic material, P and S velocities given
   if(indic == 1) then

! P and S velocity
      cp = val1
      cs = val2

! Lam'e parameters
      lambdaplus2mu = density*cp*cp
      mu = density*cs*cs
      two_mu = 2.d0*mu
      lambda = lambdaplus2mu - two_mu

! bulk modulus Kappa
      kappa = lambda + two_mu/3.d0

! Young modulus
      young = 9.d0*kappa*mu/(3.d0*kappa + mu)

! Poisson's ratio
      poisson = half*(3.d0*kappa-two_mu)/(3.d0*kappa+mu)

! Poisson's ratio must be between -1 and +1/2
      if (poisson < -1.d0 .or. poisson > 0.5d0) call exit_MPI('Poisson''s ratio out of range')

!---- anisotropic material, c11, c13, c33 and c44 given in Pascal
   else if (indic == 2) then
      c11 = val1
      c13 = val2
      c33 = val3
      c44 = val4

   else
     call exit_MPI('wrong model flag read')

   endif

!
!----  set elastic coefficients and density
!
!  Isotropic              :  lambda, mu, K (= lambda + 2*mu), zero
!  Transverse anisotropic :  c11, c13, c33, c44
!
  if(indic == 1) then
    elastcoef(1,n) = lambda
    elastcoef(2,n) = mu
    elastcoef(3,n) = lambdaplus2mu
    elastcoef(4,n) = zero
  else
    elastcoef(1,n) = c11
    elastcoef(2,n) = c13
    elastcoef(3,n) = c33
    elastcoef(4,n) = c44
  endif

  density_array(n) = density

!
!----    check the input
!
  if(indic == 1) then
! material can be acoustic (fluid) or elastic (solid)
    if(elastcoef(2,n) > TINYVAL) then
      write(iout,200) n,cp,cs,density,poisson,lambda,mu,kappa,young
    else
      write(iout,300) n,cp,density,kappa
    endif
  else
    write(iout,400) n,c11,c13,c33,c44,density,sqrt(c33/density),sqrt(c11/density),sqrt(c44/density),sqrt(c44/density)
  endif

  enddo

!
!---- formats
!
  100   format(//,' M a t e r i a l   s e t s :  ', &
         ' 2 D  e l a s t i c i t y', &
         /1x,54('='),//5x,'Number of material sets . . . . . . (numat) =',i6)

  200   format(//5x,'----------------------------------------',/5x, &
         '-- Elastic (solid) isotropic material --',/5x, &
         '----------------------------------------',/5x, &
         'Material set number. . . . . . . . (jmat) =',i6,/5x, &
         'P-wave velocity. . . . . . . . . . . (cp) =',1pe15.8,/5x, &
         'S-wave velocity. . . . . . . . . . . (cs) =',1pe15.8,/5x, &
         'Mass density. . . . . . . . . . (density) =',1pe15.8,/5x, &
         'Poisson''s ratio. . . . . . . . .(poisson) =',1pe15.8,/5x, &
         'First Lame parameter Lambda. . . (lambda) =',1pe15.8,/5x, &
         'Second Lame parameter Mu. . . . . . .(mu) =',1pe15.8,/5x, &
         'Bulk modulus Kappa . . . . . . . .(kappa) =',1pe15.8,/5x, &
         'Young''s modulus E. . . . . . . . .(young) =',1pe15.8)

  300   format(//5x,'-------------------------------',/5x, &
         '-- Acoustic (fluid) material --',/5x, &
         '-------------------------------',/5x, &
         'Material set number. . . . . . . . (jmat) =',i6,/5x, &
         'P-wave velocity. . . . . . . . . . . (cp) =',1pe15.8,/5x, &
         'Mass density. . . . . . . . . . (density) =',1pe15.8,/5x, &
         'Bulk modulus Kappa . . . . . . . .(kappa) =',1pe15.8)

  400   format(//5x,'-------------------------------------',/5x, &
         '-- Transverse anisotropic material --',/5x, &
         '-------------------------------------',/5x, &
         'Material set number. . . . . . . . (jmat) =',i6,/5x, &
         'c11 coefficient (Pascal). . . . . . (c11) =',1pe15.8,/5x, &
         'c13 coefficient (Pascal). . . . . . (c13) =',1pe15.8,/5x, &
         'c33 coefficient (Pascal). . . . . . (c33) =',1pe15.8,/5x, &
         'c44 coefficient (Pascal). . . . . . (c44) =',1pe15.8,/5x, &
         'Mass density. . . . . . . . . . (density) =',1pe15.8,/5x, &
         'Velocity of qP along vertical axis. . . . =',1pe15.8,/5x, &
         'Velocity of qP along horizontal axis. . . =',1pe15.8,/5x, &
         'Velocity of qSV along vertical axis . . . =',1pe15.8,/5x, &
         'Velocity of qSV along horizontal axis . . =',1pe15.8)

  end subroutine gmat01

