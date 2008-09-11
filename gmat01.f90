
!========================================================================
!
!                   S P E C F E M 2 D  Version 6.3
!                   ------------------------------
!
! Copyright Universite de Pau et des Pays de l'Adour and CNRS, France.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT gps DOT caltech DOT edu
!               Jeroen Tromp, jtromp aT gps DOT caltech DOT edu
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic and poroelastic wave equations
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

  subroutine gmat01(density_array,porosity_array,tortuosity_array,permeability,poroelastcoef,numat)

! read properties of a 2D isotropic or anisotropic (to be defined) linear elastic element
! velocities cpI, cpII, and cs are calculated using solid, fluid, and frame properties 

  implicit none

  include "constants.h"

  character(len=80) datlin

  integer numat
  double precision density_array(2,numat),poroelastcoef(4,3,numat),porosity_array(numat)
  double precision tortuosity_array(numat),permeability(3,numat)

  integer in,n,indic
  double precision lambdaplus2mu_s,lambdaplus2mu_fr,kappa_s,kappa_f,kappa_fr
  double precision young_s,poisson_s,density(2),phi,tortuosity,permxx,permzz,permxz
  double precision cpIsquare,cpIIsquare,cssquare,mu_s,mu_fr,eta_f,lambda_s,lambda_fr
  double precision vals(2),valf(2),valfr(2)
  double precision c11,c13,c33,c44

  double precision afactor,bfactor,cfactor,D_biot,H_biot,C_biot,M_biot,density_bar
!
!---- loop over the different material sets
!
  density_array(:,:) = zero
  porosity_array(:) = zero
  tortuosity_array(:) = zero
  permeability(:,:) = zero
  poroelastcoef(:,:,:) = zero

  write(iout,100) numat

  read(iin ,"(a80)") datlin
  do in = 1,numat

   read(iin ,*) n,indic,density(1),density(2),phi,tortuosity,permxx,permxz,permzz,vals(1),valf(1),valfr(1),vals(2),valf(2),valfr(2)

   if(n<1 .or. n>numat) call exit_MPI('Wrong material set number')

!---- isotropic material, kappa and mu/eta, for solid, fluid, and frame given
   if(indic == 1) then

! Solid properties 
      kappa_s = vals(1)
      mu_s = vals(2)
! Fluid properties 
      kappa_f = valf(1)
      eta_f = valf(2)
! Frame properties 
      kappa_fr = valfr(1)
      mu_fr = valfr(2)
! Lam'e parameters for the solid phase and the frame
      lambdaplus2mu_s = kappa_s + FOUR_THIRDS*mu_s
      lambda_s = lambdaplus2mu_s - 2.d0*mu_s
      lambdaplus2mu_fr = kappa_fr + FOUR_THIRDS*mu_fr
      lambda_fr = lambdaplus2mu_fr - 2.d0*mu_fr

! Biot coefficients for the input phi
      D_biot = kappa_s*(1.d0 + phi*(kappa_s/kappa_f - 1.d0))
      H_biot = (kappa_s - kappa_fr)*(kappa_s - kappa_fr)/(D_biot - kappa_fr) + kappa_fr + FOUR_THIRDS*mu_fr
      C_biot = kappa_s*(kappa_s - kappa_fr)/(D_biot - kappa_fr)
      M_biot = kappa_s*kappa_s/(D_biot - kappa_fr)
! Approximated velocities (no viscous dissipation)
      density_bar = (1.d0 - phi)*density(1) + phi*density(2)
      afactor = density_bar - phi/tortuosity*density(2)
      bfactor = H_biot + phi*density_bar/(tortuosity*density(2))*M_biot - 2.d0*phi/tortuosity*C_biot
      cfactor = phi/(tortuosity*density(2))*(H_biot*M_biot - C_biot*C_biot)
      cpIsquare = (bfactor + sqrt(bfactor*bfactor - 4.d0*afactor*cfactor))/(2.d0*afactor)
      cpIIsquare = (bfactor - sqrt(bfactor*bfactor - 4.d0*afactor*cfactor))/(2.d0*afactor)
      cssquare = mu_fr/afactor

! Young modulus for the solid phase 
      young_s = 9.d0*kappa_s*mu_s/(3.d0*kappa_s + mu_s)

! Poisson's ratio for the solid phase 
      poisson_s = HALF*(3.d0*kappa_s- 2.d0*mu_s)/(3.d0*kappa_s+mu_s)

! Poisson's ratio must be between -1 and +1/2
      if (poisson_s < -1.d0 .or. poisson_s > 0.5d0) stop 'Poisson''s ratio for the solid phase out of range'

!---- anisotropic material, c11, c13, c33 and c44 given in Pascal
   else if (indic == 2) then
      stop 'Attention, anisotropic still needs to be defined'
!      c11 = val1
!      c13 = val2
!      c33 = val3
!      c44 = val4

   else
     call exit_MPI('wrong model flag read')

   endif

!
!----  set elastic coefficients and density
!
!  Isotropic              :  lambda, mu, K (= lambda + 2*mu), zero for the solid phase (1) and the frame (3)
!  Transverse anisotropic :  c11, c13, c33, c44
!
  if(indic == 1) then
    poroelastcoef(1,1,n) = lambda_s
    poroelastcoef(2,1,n) = mu_s
    poroelastcoef(3,1,n) = lambdaplus2mu_s
    poroelastcoef(4,1,n) = zero

    poroelastcoef(1,2,n) = kappa_f
    poroelastcoef(2,2,n) = eta_f
    poroelastcoef(3,2,n) = zero
    poroelastcoef(4,2,n) = zero

    poroelastcoef(1,3,n) = lambda_fr
    poroelastcoef(2,3,n) = mu_fr
    poroelastcoef(3,3,n) = lambdaplus2mu_fr
    poroelastcoef(4,3,n) = zero
  else
     stop 'Attention, anisotropic still needs to be defined'
    poroelastcoef(1,1,n) = c11
    poroelastcoef(2,1,n) = c13
    poroelastcoef(3,1,n) = c33
    poroelastcoef(4,1,n) = c44
  endif

  density_array(1,n) = density(1)
  density_array(2,n) = density(2)
  porosity_array(n) = phi
  tortuosity_array(n) = tortuosity
  permeability(1,n) = permxx
  permeability(2,n) = permxz
  permeability(3,n) = permzz

!
!----    check the input
!
  if(indic == 1) then
! material can be acoustic (fluid) or poroelastic (solid/fluid) or elastic (solid)
    if(phi < TINYVAL) then ! material is elastic
      write(iout,800) n,sqrt(cpIsquare),sqrt(cssquare),density(1),poisson_s,lambda_s,mu_s,kappa_s,young_s
    elseif(phi >=1.d0)then ! material is acoustic
      write(iout,900)  n,sqrt(kappa_f/density(2)),density(2),kappa_f
    else ! material is poroelastic
      write(iout,200) n,sqrt(cpIsquare),sqrt(cpIIsquare),sqrt(cssquare)
      write(iout,300) density(1),poisson_s,lambda_s,mu_s,kappa_s,young_s
      write(iout,400) density(2),kappa_f,eta_f
      write(iout,500) lambda_fr,mu_fr,kappa_fr,phi,tortuosity,permxx,permxz,permzz
      write(iout,600) D_biot,H_biot,C_biot,M_biot
    endif
  else
     stop 'Attention, anisotropic still needs to be defined'
    write(iout,700) n,c11,c13,c33,c44,density(1),sqrt(c33/density(1)),sqrt(c11/density(1)),sqrt(c44/density(1)),sqrt(c44/density(1))
  endif

  enddo

!
!---- formats
!
  100   format(//,' M a t e r i a l   s e t s :  ', &
         ' 2 D  (p o r o) e l a s t i c i t y', &
         /1x,54('='),//5x,'Number of material sets . . . . . . (numat) =',i6)

  200   format(//5x,'----------------------------------------',/5x, &
         '-- Poroelastic isotropic material --',/5x, &
         '----------------------------------------',/5x, &
         'Material set number. . . . . . . . (jmat) =',i6,/5x, &
         'First P-wave velocity. . . . . . . . . . . (cpI) =',1pe15.8,/5x, &
         'Second P-wave velocity. . . . . . . . . . . (cpII) =',1pe15.8,/5x, &
         'S-wave velocity. . . . . . . . . . . (cs) =',1pe15.8)

  300   format(//5x,'-------------------------------',/5x, &
         '-- Solid phase properties --',/5x, &
         'Mass density. . . . . . . . . . (density_s) =',1pe15.8,/5x, &
         'Poisson''s ratio. . . . . . . . .(poisson_s) =',1pe15.8,/5x, &
         'First Lame parameter Lambda. . . (lambda_s) =',1pe15.8,/5x, &
         'Second Lame parameter Mu. . . . . . .(mu_s) =',1pe15.8,/5x, &
         'Solid bulk modulus Kappa . . . . . . . .(kappa_s) =',1pe15.8,/5x, &
         'Young''s modulus E. . . . . . . . .(young_s) =',1pe15.8)

  400   format(//5x,'-------------------------------',/5x, &
         '-- Fluid phase properties --',/5x, &
         'Mass density. . . . . . . . . . (density_f) =',1pe15.8,/5x, &
         'Fluid bulk modulus Kappa . . . . . . . .(kappa_f) =',1pe15.8,/5x, &
         'Fluid viscosity Eta . . . . . . . .(eta_f) =',1pe15.8)

  500   format(//5x,'-------------------------------',/5x, &
         '-- Frame properties --',/5x, &
         'First Lame parameter Lambda. . . (lambda_fr) =',1pe15.8,/5x, &
         'Second Lame parameter Mu. . . . . . .(mu_fr) =',1pe15.8,/5x, &
         'Frame bulk modulus Kappa . . . . . . . .(kappa_fr) =',1pe15.8,/5x, &
         'Porosity. . . . . . . . . . . . . . . . .(phi) =',1pe15.8,/5x,&
         'Tortuosity. . . . . . . . . . . . . . . . . . =',1pe15.8,/5x,&
         'Permeability xx component. . . . . . . . . . =',1pe15.8,/5x,&
         'Permeability zx component. . . . . . . . . . =',1pe15.8,/5x,&
         'Permeability zz component. . . . . . . . . . =',1pe15.8)

  600   format(//5x,'-------------------------------',/5x, &
         '-- Biot coefficients --',/5x, &
         '-------------------------------',/5x, &
         'D. . . . . . . . =',1pe15.8,/5x, &
         'H. . . . . . . . =',1pe15.8,/5x, &
         'C. . . . . . . . =',1pe15.8,/5x, &
         'M. . . . . . . . =',1pe15.8)

  700   format(//5x,'-------------------------------------',/5x, &
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

  800   format(//5x,'--------------------------------------------------',/5x, &
         '-- Elastic (solid - phi = 0) isotropic material --',/5x, &
         '--------------------------------------------------',/5x, &
         'Material set number. . . . . . . . (jmat) =',i6,/5x, &
         'P-wave velocity. . . . . . . . . . . (cp) =',1pe15.8,/5x, &
         'S-wave velocity. . . . . . . . . . . (cs) =',1pe15.8,/5x, &
         'Mass density. . . . . . . . . . (density) =',1pe15.8,/5x, &
         'Poisson''s ratio. . . . . . . . .(poisson) =',1pe15.8,/5x, &
         'First Lame parameter Lambda. . . (lambda) =',1pe15.8,/5x, &
         'Second Lame parameter Mu. . . . . . .(mu) =',1pe15.8,/5x, &
         'Bulk modulus Kappa . . . . . . . .(kappa) =',1pe15.8,/5x, &
         'Young''s modulus E. . . . . . . . .(young) =',1pe15.8)

  900   format(//5x,'-------------------------------',/5x, &
         '-- Acoustic (fluid) material --',/5x, &
         '-------------------------------',/5x, &
         'Material set number. . . . . . . . (jmat) =',i6,/5x, &
         'P-wave velocity. . . . . . . . . . . (cp) =',1pe15.8,/5x, &
         'Mass density. . . . . . . . . . (density) =',1pe15.8,/5x, &
         'Bulk modulus Kappa . . . . . . . .(kappa) =',1pe15.8)

  end subroutine gmat01

