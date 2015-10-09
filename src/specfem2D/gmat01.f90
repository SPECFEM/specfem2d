
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and Inria at the following URL
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

  subroutine gmat01(f0)

! reads properties of a 2D isotropic or anisotropic linear elastic element

  use specfem_par, only : AXISYM,ANY_ANISOTROPY,density,porosity,tortuosity, &
                          anisotropy,permeability,poroelastcoef, &
                          numat,myrank,QKappa_attenuation,Qmu_attenuation, &
                          freq0,Q0,ATTENUATION_PORO_FLUID_PART,assign_external_model, &
                          tomo_material

  implicit none
  include "constants.h"

  double precision f0

  ! local parameters
  double precision :: lambdaplus2mu,kappa
  double precision :: young,poisson,cp,cs,mu,two_mu,lambda,QKappa,Qmu
  double precision :: lambdaplus2mu_s,lambdaplus2mu_fr,kappa_s,kappa_f,kappa_fr
  double precision :: young_s,poisson_s,density_mat(2),phi,tortuosity_mat
  double precision :: cpIsquare,cpIIsquare,cssquare,mu_s,mu_fr,eta_f,lambda_s,lambda_fr
  double precision :: val1,val2,val3,val4,val5,val6
  double precision :: val7,val8,val9,val10,val11,val12,val0
  double precision ::  c11,c13,c15,c33,c35,c55,c12,c23,c25,c22
  double precision :: D_biot,H_biot,C_biot,M_biot
  double precision :: w_c
  integer in,n,indic
  character(len=80) datlin


  ANY_ANISOTROPY = .false.

  !
  !---- loop over the different material sets
  !
  density(:,:) = zero
  porosity(:) = zero
  tortuosity(:) = zero
  anisotropy(:,:) = zero
  permeability(:,:) = zero
  poroelastcoef(:,:,:) = zero
  QKappa_attenuation(:) = 9999.
  Qmu_attenuation(:) = 9999.
  tomo_material = 0 ! Index of the material that will be defined by an external tomo file if needed (TOMOGRAPHY_FILE)

  if(myrank == 0) write(IOUT,100) numat

  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  do in = 1,numat

     read(IIN,*) n,indic,val0,val1,val2,val3,val4,val5,val6,val7,val8,val9,val10,val11,val12

     if(n<1 .or. n>numat) call exit_MPI('Wrong material set number')

     !---- isotropic material, P and S velocities given, allows for declaration of elastic/acoustic material
     !---- elastic (cs/=0) and acoustic (cs=0)
     if(indic == 1) then
        density_mat(1) = val0

        ! P and S velocity
        cp = val1
        cs = val2

        ! QKappa and Qmu values
        QKappa = val5
        Qmu = val6
        if(QKappa <= 0.0000001 .or. Qmu <= 0.0000001) &
          stop 'negative or null values of Q attenuation factor not allowed; set them equal to 9999 to indicate no attenuation'

        ! Lame parameters
        lambdaplus2mu = density_mat(1)*cp*cp
        mu = density_mat(1)*cs*cs
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

        density_mat(1) =val0

        ! Anisotropy parameters
        c11 = val1
        c13 = val2
        c15 = val3
        c33 = val4
        c35 = val5
        c55 = val6
        c12 = val7
        c23 = val8
        c25 = val9
        c22 = val10  ! This value is used for AXISYM only

        ! P and S velocity
        cp = sqrt(c33/density_mat(1))
        cs = sqrt(c55/density_mat(1))

        ! Lame parameters
        lambdaplus2mu = density_mat(1)*cp*cp
        mu = density_mat(1)*cs*cs
        two_mu = 2.d0*mu
        lambda = lambdaplus2mu - two_mu

        ! bulk modulus Kappa
        kappa = lambda + two_mu/3.d0

        ! Young modulus
        young = 9.d0*kappa*mu/(3.d0*kappa + mu)

        ! Poisson's ratio
        poisson = half*(3.d0*kappa-two_mu)/(3.d0*kappa+mu)

        ANY_ANISOTROPY = .true.


        !---- isotropic material, moduli are given, allows for declaration of poroelastic material
        !---- poroelastic (0<phi<1)
     else if (indic == 3) then
        ! Qmu values
        Qmu = val12

        density_mat(1) =val0
        density_mat(2) =val1

        ! Solid properties
        kappa_s = val7
        mu_s = val11
        ! Fluid properties
        kappa_f = val8
        eta_f = val10
        ! Frame properties
        kappa_fr = val9
        mu_fr = val11
        ! Lame parameters for the solid phase and the frame
        lambdaplus2mu_s = kappa_s + FOUR_THIRDS*mu_s
        lambda_s = lambdaplus2mu_s - 2.d0*mu_s
        lambdaplus2mu_fr = kappa_fr + FOUR_THIRDS*mu_fr
        lambda_fr = lambdaplus2mu_fr - 2.d0*mu_fr
        phi = val2
        tortuosity_mat = val3

        ! Biot coefficients for the input phi
        D_biot = kappa_s*(1.d0 + phi*(kappa_s/kappa_f - 1.d0))
        H_biot = (kappa_s - kappa_fr)*(kappa_s - kappa_fr)/(D_biot - kappa_fr) + kappa_fr + FOUR_THIRDS*mu_fr
        C_biot = kappa_s*(kappa_s - kappa_fr)/(D_biot - kappa_fr)
        M_biot = kappa_s*kappa_s/(D_biot - kappa_fr)

        call get_poroelastic_velocities(cpIsquare,cpIIsquare,cssquare, &
                                  H_biot,C_biot,M_biot,mu_fr,phi, &
                                  tortuosity_mat,density_mat(1),density_mat(2),eta_f, &
                                  val4,f0,freq0,Q0,w_c,ATTENUATION_PORO_FLUID_PART)

        porosity(n) = val2
        tortuosity(n) = val3
        permeability(1,n) = val4
        permeability(2,n) = val5
        permeability(3,n) = val6

        ! Young modulus for the solid phase
        young_s = 9.d0*kappa_s*mu_s/(3.d0*kappa_s + mu_s)

        ! Poisson's ratio for the solid phase
        poisson_s = HALF*(3.d0*kappa_s- 2.d0*mu_s)/(3.d0*kappa_s+mu_s)

        ! Poisson's ratio must be between -1 and +1/2
        if (poisson_s < -1.d0 .or. poisson_s > 0.5d0) stop 'Poisson''s ratio for the solid phase out of range'
     else if (indic <= 0) then
       assign_external_model = .true.
       tomo_material = n
       mu = val2 ! for acoustic medium vs must be 0 anyway
     else
        call exit_MPI('wrong model flag read')

     endif

     !
     !----  set elastic coefficients and density
     !
     if(indic == 1) then
        density(1,n) = density_mat(1)
        poroelastcoef(1,1,n) = lambda
        poroelastcoef(2,1,n) = mu
        poroelastcoef(3,1,n) = lambdaplus2mu
        poroelastcoef(4,1,n) = zero
        QKappa_attenuation(n) = QKappa
        Qmu_attenuation(n) = Qmu
        if(mu > TINYVAL) then
           porosity(n) = 0.d0
        else
           porosity(n) = 1.d0
        endif
     else if (indic == 2) then
        density(1,n) = density_mat(1)
! dummy poroelastcoef values, trick to avoid floating invalid
        poroelastcoef(1,1,n) = lambda
        poroelastcoef(2,1,n) = mu
        poroelastcoef(3,1,n) = lambdaplus2mu
        poroelastcoef(4,1,n) = zero
        anisotropy(1,n) = c11
        anisotropy(2,n) = c13
        anisotropy(3,n) = c15
        anisotropy(4,n) = c33
        anisotropy(5,n) = c35
        anisotropy(6,n) = c55
        anisotropy(7,n) = c12
        anisotropy(8,n) = c23
        anisotropy(9,n) = c25
        anisotropy(10,n) = c22 ! This value is used for AXISYM only
        porosity(n) = 0.d0
     else if (indic == 3) then
        density(1,n) = density_mat(1)
        density(2,n) = density_mat(2)
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
     else if (indic <= 0) then
        ! Assign dummy values for now (for acoustic medium vs must be 0 anyway), these values will be read in read_external_model
        density(1,n) = -1.0d0
        poroelastcoef(1,1,n) = -1.0d0
        poroelastcoef(2,1,n) = mu
        poroelastcoef(3,1,n) = -1.0d0
        poroelastcoef(4,1,n) = zero
        QKappa_attenuation(n) = 9999.
        Qmu_attenuation(n) = 9999.
        if(mu > TINYVAL) then
           porosity(n) = 0.d0
        else
           porosity(n) = 1.d0
        endif
     else
        call exit_MPI('wrong model flag read')
     endif

     !
     !----    check what has been read
     !
     if(myrank == 0) then
        if(indic == 1) then
           ! material can be acoustic (fluid) or elastic (solid)
           if(poroelastcoef(2,1,n) > TINYVAL) then    ! elastic
              write(IOUT,200) n,cp,cs,density_mat(1),poisson,lambda,mu,kappa,young,QKappa,Qmu
              if(poisson < 0.d0) then
                write(IOUT,*)
                write(IOUT,*) 'Materials with a negative Poisson''s ratio can exist,'
                write(IOUT,*) 'see e.g. R. Lakes, "Science" vol. 235, p. 1038-1040 (1987),'
                write(IOUT,*) 'but are extremely rare.'
                write(IOUT,*) 'Hope you know what you are doing...'
                write(IOUT,*)
              endif
           else                                       ! acoustic
              write(IOUT,300) n,cp,density_mat(1),kappa,QKappa,Qmu
           endif
        else if(indic == 2) then                      ! elastic (anisotropic)
          if (AXISYM) then
            write(IOUT,450) n,density_mat(1),c11,c13,c15,c33,c35,c55,c12,c23,c25,c22
          else
            write(IOUT,400) n,density_mat(1),c11,c13,c15,c33,c35,c55,c12,c23,c25
          endif
        else if(indic == 3) then
           ! material is poroelastic (solid/fluid)
           write(iout,500) n,sqrt(cpIsquare),sqrt(cpIIsquare),sqrt(cssquare)
           write(iout,600) density_mat(1),poisson_s,lambda_s,mu_s,kappa_s,young_s
           if(poisson_s < 0.d0) then
             write(IOUT,*)
             write(IOUT,*) 'Materials with a negative Poisson''s ratio can exist,'
             write(IOUT,*) 'see e.g. R. Lakes, "Science" vol. 235, p. 1038-1040 (1987),'
             write(IOUT,*) 'but are extremely rare.'
             write(IOUT,*) 'Hope you know what you are doing...'
             write(IOUT,*)
           endif
           write(iout,700) density_mat(2),kappa_f,eta_f
           write(iout,800) lambda_fr,mu_fr,kappa_fr,porosity(n),tortuosity(n),&
                permeability(1,n),permeability(2,n),permeability(3,n),Qmu
           write(iout,900) D_biot,H_biot,C_biot,M_biot,w_c
        else if (indic <= 0) then
          write(IOUT,1000) n
        endif
     endif

  enddo

  !
  !---- formats
  !

100 format(//,' M a t e r i a l   s e t s :  ', &
       ' 2 D  (p o r o) e l a s t i c i t y', &
       /1x,54('='),//5x,'Number of material sets . . . . . . (numat) =',i6)

200 format(//5x,'----------------------------------------',/5x, &
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
       'QKappa_attenuation .  . . . . . .(QKappa) =',1pe15.8,/5x, &
       'Qmu_attenuation . . . . . . . . . . (Qmu) =',1pe15.8)

300 format(//5x,'-------------------------------',/5x, &
       '-- Acoustic (fluid) material --',/5x, &
       '-------------------------------',/5x, &
       'Material set number. . . . . . . . (jmat) =',i6,/5x, &
       'P-wave velocity. . . . . . . . . . . (cp) =',1pe15.8,/5x, &
       'Mass density. . . . . . . . . . (density) =',1pe15.8,/5x, &
       'Bulk modulus Kappa . . . . . . . .(kappa) =',1pe15.8,/5x, &
       'QKappa_attenuation. . . . . . . .(QKappa) =',1pe15.8,/5x, &
       'Qmu_attenuation. . . . . . . . . . .(Qmu) =',1pe15.8)

400 format(//5x,'-------------------------------------',/5x, &
       '-- Anisotropic material --',/5x, &
       '-------------------------------------',/5x, &
       'Material set number. . . . . . . . (jmat) =',i6,/5x, &
       'Mass density. . . . . . . . . . (density) =',1pe15.8,/5x, &
       'c11 coefficient (Pascal). . . . . . (c11) =',1pe15.8,/5x, &
       'c13 coefficient (Pascal). . . . . . (c13) =',1pe15.8,/5x, &
       'c15 coefficient (Pascal). . . . . . (c15) =',1pe15.8,/5x, &
       'c33 coefficient (Pascal). . . . . . (c33) =',1pe15.8,/5x, &
       'c35 coefficient (Pascal). . . . . . (c35) =',1pe15.8,/5x, &
       'c55 coefficient (Pascal). . . . . . (c55) =',1pe15.8,/5x, &
       'c12 coefficient (Pascal). . . . . . (c12) =',1pe15.8,/5x, &
       'c23 coefficient (Pascal). . . . . . (c23) =',1pe15.8,/5x, &
       'c25 coefficient (Pascal). . . . . . (c25) =',1pe15.8,/5x)

450 format(//5x,'-------------------------------------',/5x, &
       '-- Axisymetrical anisotropic material --',/5x, &
       '-------------------------------------',/5x, &
       'Material set number. . . . . . . . (jmat) =',i6,/5x, &
       'Mass density. . . . . . . . . . (density) =',1pe15.8,/5x, &
       'c11 coefficient (Pascal). . . . . . (c11) =',1pe15.8,/5x, &
       'c13 coefficient (Pascal). . . . . . (c13) =',1pe15.8,/5x, &
       'c15 coefficient (Pascal). . . . . . (c15) =',1pe15.8,/5x, &
       'c33 coefficient (Pascal). . . . . . (c33) =',1pe15.8,/5x, &
       'c35 coefficient (Pascal). . . . . . (c35) =',1pe15.8,/5x, &
       'c55 coefficient (Pascal). . . . . . (c55) =',1pe15.8,/5x, &
       'c12 coefficient (Pascal). . . . . . (c12) =',1pe15.8,/5x, &
       'c23 coefficient (Pascal). . . . . . (c23) =',1pe15.8,/5x, &
       'c25 coefficient (Pascal). . . . . . (c25) =',1pe15.8,/5x, &
       'c22 coefficient (Pascal). . . . . . (c22) =',1pe15.8,/5x)

500 format(//5x,'----------------------------------------',/5x, &
       '-- Poroelastic isotropic material --',/5x, &
       '----------------------------------------',/5x, &
       'Material set number. . . . . . . . (jmat) =',i6,/5x, &
       'First P-wave velocity. . . . . . . .(cpI) =',1pe15.8,/5x, &
       'Second P-wave velocity. . . . . . .(cpII) =',1pe15.8,/5x, &
       'S-wave velocity. . . . . . . . . . . (cs) =',1pe15.8)

600 format(//5x,'-------------------------------',/5x, &
       '-- Solid phase properties --',/5x, &
       'Mass density. . . . . . . . . . (density_s) =',1pe15.8,/5x, &
       'Poisson''s ratio . . . . . . . . (poisson_s) =',1pe15.8,/5x, &
       'First Lame parameter Lambda. . . (lambda_s) =',1pe15.8,/5x, &
       'Second Lame parameter Mu. . . . . . .(mu_s) =',1pe15.8,/5x, &
       'Solid bulk modulus Kappa . . . . .(kappa_s) =',1pe15.8,/5x, &
       'Young''s modulus E . . . . . . .  .(young_s) =',1pe15.8)

700 format(//5x,'-------------------------------',/5x, &
       '-- Fluid phase properties --',/5x, &
       'Mass density. . . . . . . . . . (density_f) =',1pe15.8,/5x, &
       'Fluid bulk modulus Kappa . . . . .(kappa_f) =',1pe15.8,/5x, &
       'Fluid viscosity Eta . . . . . . . . (eta_f) =',1pe15.8)

800 format(//5x,'-------------------------------',/5x, &
       '-- Frame properties --',/5x, &
       'First Lame parameter Lambda. . . (lambda_fr) =',1pe15.8,/5x, &
       'Second Lame parameter Mu. . . . . . .(mu_fr) =',1pe15.8,/5x, &
       'Frame bulk modulus Kappa . . . . .(kappa_fr) =',1pe15.8,/5x, &
       'Porosity. . . . . . . . . . . . . . . .(phi) =',1pe15.8,/5x,&
       'Tortuosity. . . . . . . . . . . . . . . .(c) =',1pe15.8,/5x,&
       'Permeability xx component. . . . . . . . . . =',1pe15.8,/5x,&
       'Permeability zx component. . . . . . . . . . =',1pe15.8,/5x,&
       'Permeability zz component. . . . . . . . . . =',1pe15.8,/5x,&
       'Qmu_attenuation. . . . . . . . . . . . (Qmu) =',1pe15.8)

900 format(//5x,'-------------------------------',/5x, &
         '-- Biot coefficients --',/5x, &
         '-------------------------------',/5x, &
         'D. . . . . . . . . .=',1pe15.8,/5x, &
         'H. . . . . . . . . .=',1pe15.8,/5x, &
         'C. . . . . . . . . .=',1pe15.8,/5x, &
         'M. . . . . . . . . .=',1pe15.8,/5x, &
         'Characteristic freq =',1pe15.8)

1000 format(//5x,'----------------------------------------------------',/5x, &
       '-- Material described by external tomography file --',/5x, &
       '----------------------------------------------------',/5x, &
       'Material set number. . . . . . . . (jmat) =',i6,/5x)

  end subroutine gmat01

