
!========================================================================
!
!                   S P E C F E M 2 D  Version 6.0
!                   ------------------------------
!
! Copyright Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France,
! and Princeton University, USA.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
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

  subroutine gmat01(density_array,porosity_array,tortuosity_array,permeability,poroelastcoef,&
                    numat,myrank,ipass,Qp_array,Qs_array)

! read properties of a 2D isotropic or anisotropic linear elastic element

  implicit none

  include "constants.h"

  character(len=80) datlin
  double precision lambdaplus2mu,kappa

  integer numat,myrank,ipass
  double precision density_array(2,numat),poroelastcoef(4,3,numat),porosity_array(numat)
  double precision tortuosity_array(numat),permeability(3,numat)
  double precision Qp_array(numat),Qs_array(numat)

  integer in,n,indic
  double precision young,poisson,cp,cs,mu,two_mu,lambda,Qp,Qs
  double precision lambdaplus2mu_s,lambdaplus2mu_fr,kappa_s,kappa_f,kappa_fr
  double precision young_s,poisson_s,density(2),phi,tortuosity,permxx,permzz,permxz
  double precision cpIsquare,cpIIsquare,cssquare,mu_s,mu_fr,eta_f,lambda_s,lambda_fr
  double precision val1,val2,val3,val4,val5,val6
  double precision val7,val8,val9,val10,val11,val12,val0
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
  Qp_array(:) = zero
  Qs_array(:) = zero

  if(myrank == 0 .and. ipass == 1) write(IOUT,100) numat

  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  do in = 1,numat

   read(IIN,*) n,indic,val0,val1,val2,val3,val4,val5,val6,val7,val8,val9,val10,val11,val12

   if(n<1 .or. n>numat) call exit_MPI('Wrong material set number')

!---- isotropic material, P and S velocities given, allows for declaration of elastic/acoustic material
!---- elastic (cs/=0) and acoustic (cs=0)
   if(indic == 1) then
      density(1) = val0

! P and S velocity
      cp = val1
      cs = val2

! Qp and Qs values
      Qp = val5
      Qs = val6

! Lam'e parameters
      lambdaplus2mu = density(1)*cp*cp
      mu = density(1)*cs*cs
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

      density(1) =val0
      c11 = val1
      c13 = val2
      c33 = val3
      c44 = val4

!---- isotropic material, moduli are given, allows for declaration of poroelastic material
!---- poroelastic (<0phi<1)
   else if (indic == 3) then
! Qs values
      Qs = val12

      density(1) =val0
      density(2) =val1

! Solid properties
      kappa_s = val7
      mu_s = val11
! Fluid properties
      kappa_f = val8
      eta_f = val10
! Frame properties
      kappa_fr = val9
      mu_fr = val11
! Lam'e parameters for the solid phase and the frame
      lambdaplus2mu_s = kappa_s + FOUR_THIRDS*mu_s
      lambda_s = lambdaplus2mu_s - 2.d0*mu_s
      lambdaplus2mu_fr = kappa_fr + FOUR_THIRDS*mu_fr
      lambda_fr = lambdaplus2mu_fr - 2.d0*mu_fr
      phi = val2
      tortuosity = val3

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
      cssquare = val11/afactor

  porosity_array(n) = val2
  tortuosity_array(n) = val3
  permeability(1,n) = val4
  permeability(2,n) = val5
  permeability(3,n) = val6

! Young modulus for the solid phase
      young_s = 9.d0*kappa_s*mu_s/(3.d0*kappa_s + mu_s)

! Poisson's ratio for the solid phase
      poisson_s = HALF*(3.d0*kappa_s- 2.d0*mu_s)/(3.d0*kappa_s+mu_s)

! Poisson's ratio must be between -1 and +1/2
      if (poisson_s < -1.d0 .or. poisson_s > 0.5d0) stop 'Poisson''s ratio for the solid phase out of range'

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
    density_array(1,n) = density(1)
    poroelastcoef(1,1,n) = lambda
    poroelastcoef(2,1,n) = mu
    poroelastcoef(3,1,n) = lambdaplus2mu
    poroelastcoef(4,1,n) = zero
  Qp_array(n) = Qp
  Qs_array(n) = Qs
     if(mu > TINYVAL) then
  porosity_array(n) = 0.d0
     else
  porosity_array(n) = 1.d0
     endif
  elseif(indic == 2) then
    density_array(1,n) = density(1)
    poroelastcoef(1,1,n) = c11
    poroelastcoef(2,1,n) = c13
    poroelastcoef(3,1,n) = c33
    poroelastcoef(4,1,n) = c44
  Qp_array(n) = Qp
  Qs_array(n) = Qs
  porosity_array(n) = 0.d0
  else
    density_array(1,n) = density(1)
    density_array(2,n) = density(2)
    poroelastcoef(1,1,n) = lambda_s
    poroelastcoef(2,1,n) = mu_s    ! = mu_fr
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
  Qp_array(n) = 10.d0 ! dummy for attenuation_model
  Qs_array(n) = Qs
  endif

!
!----    check what has been read
!
  if(myrank == 0 .and. ipass == 1) then
  if(indic == 1) then
! material can be acoustic (fluid) or elastic (solid)
    if(poroelastcoef(2,1,n) > TINYVAL) then    ! elastic
      write(IOUT,200) n,cp,cs,density(1),poisson,lambda,mu,kappa,young,Qp,Qs
    else                                       ! acoustic
      write(IOUT,300) n,cp,density(1),kappa,Qp,Qs
    endif
  elseif(indic == 2) then                      ! elastic (anisotropic)
    write(IOUT,400) n,c11,c13,c33,c44,density,sqrt(c33/density),sqrt(c11/density),sqrt(c44/density),sqrt(c44/density),Qp,Qs
  else
! material is poroelastic (solid/fluid)
      write(iout,500) n,sqrt(cpIsquare),sqrt(cpIIsquare),sqrt(cssquare)
      write(iout,600) density(1),poisson_s,lambda_s,mu_s,kappa_s,young_s
      write(iout,700) density(2),kappa_f,eta_f
      write(iout,800) lambda_fr,mu_fr,kappa_fr,porosity_array(n),tortuosity_array(n),&
                              permeability(1,n),permeability(2,n),permeability(3,n),Qs
      write(iout,900) D_biot,H_biot,C_biot,M_biot
  endif
  endif

  enddo

!
!---- formats
!
  100   format(//,' M a t e r i a l   s e t s :  ', &
         ' 2 D  (p o r o) e l a s t i c i t y', &
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
         'Young''s modulus E. . . . . . . . .(young) =',1pe15.8,/5x, &
         'Qp_attenuation. . . . . . . . . . . .(Qp) =',1pe15.8,/5x, &
         'Qs_attenuation. . . . . . . . . . . .(Qs) =',1pe15.8)

  300   format(//5x,'-------------------------------',/5x, &
         '-- Acoustic (fluid) material --',/5x, &
         '-------------------------------',/5x, &
         'Material set number. . . . . . . . (jmat) =',i6,/5x, &
         'P-wave velocity. . . . . . . . . . . (cp) =',1pe15.8,/5x, &
         'Mass density. . . . . . . . . . (density) =',1pe15.8,/5x, &
         'Bulk modulus Kappa . . . . . . . .(kappa) =',1pe15.8,/5x, &
         'Qp_attenuation. . . . . . . . . . . .(Qp) =',1pe15.8,/5x, &
         'Qs_attenuation. . . . . . . . . . . .(Qs) =',1pe15.8)

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
         'Velocity of qSV along horizontal axis . . =',1pe15.8,/5x, &
         'Qp_attenuation. . . . . . . . . . . .(Qp) =',1pe15.8,/5x, &
         'Qs_attenuation. . . . . . . . . . . .(Qs) =',1pe15.8)

  500   format(//5x,'----------------------------------------',/5x, &
         '-- Poroelastic isotropic material --',/5x, &
         '----------------------------------------',/5x, &
         'Material set number. . . . . . . . (jmat) =',i6,/5x, &
         'First P-wave velocity. . . . . . . . . . . (cpI) =',1pe15.8,/5x, &
         'Second P-wave velocity. . . . . . . . . . . (cpII) =',1pe15.8,/5x, &
         'S-wave velocity. . . . . . . . . . . (cs) =',1pe15.8)

  600   format(//5x,'-------------------------------',/5x, &
         '-- Solid phase properties --',/5x, &
         'Mass density. . . . . . . . . . (density_s) =',1pe15.8,/5x, &
         'Poisson''s ratio. . . . . . . . .(poisson_s) =',1pe15.8,/5x, &
         'First Lame parameter Lambda. . . (lambda_s) =',1pe15.8,/5x, &
         'Second Lame parameter Mu. . . . . . .(mu_s) =',1pe15.8,/5x, &
         'Solid bulk modulus Kappa . . . . . . . .(kappa_s) =',1pe15.8,/5x, &
         'Young''s modulus E. . . . . . . . .(young_s) =',1pe15.8)

  700   format(//5x,'-------------------------------',/5x, &
         '-- Fluid phase properties --',/5x, &
         'Mass density. . . . . . . . . . (density_f) =',1pe15.8,/5x, &
         'Fluid bulk modulus Kappa . . . . . . . .(kappa_f) =',1pe15.8,/5x, &
         'Fluid viscosity Eta . . . . . . . .(eta_f) =',1pe15.8)

  800   format(//5x,'-------------------------------',/5x, &
         '-- Frame properties --',/5x, &
         'First Lame parameter Lambda. . . (lambda_fr) =',1pe15.8,/5x, &
         'Second Lame parameter Mu. . . . . . .(mu_fr) =',1pe15.8,/5x, &
         'Frame bulk modulus Kappa . . . . . . . .(kappa_fr) =',1pe15.8,/5x, &
         'Porosity. . . . . . . . . . . . . . . . .(phi) =',1pe15.8,/5x,&
         'Tortuosity. . . . . . . . . . . . . . . . .(c) =',1pe15.8,/5x,&
         'Permeability xx component. . . . . . . . . . =',1pe15.8,/5x,&
         'Permeability zx component. . . . . . . . . . =',1pe15.8,/5x,&
         'Permeability zz component. . . . . . . . . . =',1pe15.8,/5x,&
         'Qs_attenuation. . . . . . . . . . . .(Qs) =',1pe15.8)

  900   format(//5x,'-------------------------------',/5x, &
         '-- Biot coefficients --',/5x, &
         '-------------------------------',/5x, &
         'D. . . . . . . . =',1pe15.8,/5x, &
         'H. . . . . . . . =',1pe15.8,/5x, &
         'C. . . . . . . . =',1pe15.8,/5x, &
         'M. . . . . . . . =',1pe15.8)

  end subroutine gmat01

