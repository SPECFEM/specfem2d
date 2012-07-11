
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
! Copyright CNRS, INRIA and University of Pau, France,
! and Princeton University / California Institute of Technology, USA.
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

  subroutine compute_forces_acoustic(nglob,nspec,nelemabs,numat,it,NSTEP, &
               anyabs,assign_external_model,ibool,kmato,numabs, &
               elastic,poroelastic,codeabs,potential_dot_dot_acoustic,potential_dot_acoustic, &
               potential_acoustic, &
               density,poroelastcoef,xix,xiz,gammax,gammaz,jacobian, &
               vpext,rhoext,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll, &
               ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
               jbegin_left,jend_left,jbegin_right,jend_right, &
               SIMULATION_TYPE,SAVE_FORWARD,nspec_left,nspec_right,&
               nspec_bottom,nspec_top,ib_left,ib_right,ib_bottom,ib_top, &
               b_absorb_acoustic_left,b_absorb_acoustic_right, &
               b_absorb_acoustic_bottom,b_absorb_acoustic_top,IS_BACKWARD_FIELD,&
            is_PML,nspec_PML,npoin_PML,ibool_PML,spec_to_PML,which_PML_elem,&
            K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store,&
            rmemory_potential_acoustic,&
            rmemory_acoustic_dux_dx,rmemory_acoustic_dux_dz,&
            deltat,PML_BOUNDARY_CONDITIONS)

! compute forces for the acoustic elements

  implicit none

  include "constants.h"

  integer :: nglob,nspec,nelemabs,numat,it,NSTEP,SIMULATION_TYPE

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: kmato
  integer, dimension(nelemabs) :: numabs,ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
               jbegin_left,jend_left,jbegin_right,jend_right

  logical, dimension(nspec) :: elastic,poroelastic
  logical, dimension(4,nelemabs)  :: codeabs

  real(kind=CUSTOM_REAL), dimension(nglob) :: &
    potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic

  double precision, dimension(2,numat) :: density
  double precision, dimension(4,3,numat) :: poroelastcoef
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz,jacobian
  double precision, dimension(NGLLX,NGLLZ,nspec) :: vpext,rhoext

  logical :: anyabs,assign_external_model
  logical :: SAVE_FORWARD,IS_BACKWARD_FIELD

! derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz

! Gauss-Lobatto-Legendre weights
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: wxgll
  real(kind=CUSTOM_REAL), dimension(NGLLZ) :: wzgll

  integer :: nspec_left,nspec_right,nspec_bottom,nspec_top
  integer, dimension(nelemabs) :: ib_left
  integer, dimension(nelemabs) :: ib_right
  integer, dimension(nelemabs) :: ib_bottom
  integer, dimension(nelemabs) :: ib_top

  real(kind=CUSTOM_REAL), dimension(NGLLZ,nspec_left,NSTEP) :: b_absorb_acoustic_left
  real(kind=CUSTOM_REAL), dimension(NGLLZ,nspec_right,NSTEP) :: b_absorb_acoustic_right
  real(kind=CUSTOM_REAL), dimension(NGLLX,nspec_top,NSTEP) :: b_absorb_acoustic_top
  real(kind=CUSTOM_REAL), dimension(NGLLX,nspec_bottom,NSTEP) :: b_absorb_acoustic_bottom

!---
!--- local variables
!---

  integer :: ispec,i,j,k,iglob,ispecabs,ibegin,iend,jbegin,jend

! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,dux_dxl,dux_dzl
  real(kind=CUSTOM_REAL) :: weight,xxi,zxi,xgamma,zgamma,jacobian1D
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2

! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

! material properties of the elastic medium
  real(kind=CUSTOM_REAL) :: mul_relaxed,lambdal_relaxed,kappal,cpl,rhol

  integer :: ifirstelem,ilastelem

!CPML coefficients and memory variables
  integer :: nspec_PML,npoin_PML,iPML,ispec_PML
  logical, dimension(4,nspec) :: which_PML_elem
  logical, dimension(nspec) :: is_PML
  integer, dimension(nspec) :: spec_to_PML
  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool_PML

  real(kind=CUSTOM_REAL), dimension(2,NGLLX,NGLLZ,nspec_PML) :: rmemory_potential_acoustic
  real(kind=CUSTOM_REAL), dimension(2,NGLLX,NGLLZ,nspec_PML) :: &
    rmemory_acoustic_dux_dx,rmemory_acoustic_dux_dz
  real(kind=CUSTOM_REAL), dimension(npoin_PML) :: K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_PML) :: potential_dot_dot_acoustic_PML
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_PML) ::PML_dux_dxl,PML_dux_dzl,&
                           PML_dux_dxl_new,PML_dux_dzl_new
  real(kind=CUSTOM_REAL) :: coef0_x, coef1_x, coef2_x, coef0_z, coef1_z, coef2_z,bb,aa,&
                            deltat

  logical :: PML_BOUNDARY_CONDITIONS

  ifirstelem = 1
  ilastelem = nspec

! loop over spectral elements
  do ispec = ifirstelem,ilastelem

!---
!--- acoustic spectral element
!---
    if(.not. elastic(ispec) .and. .not. poroelastic(ispec)) then

      rhol = density(1,kmato(ispec))

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX

          ! derivative along x and along z
          dux_dxi = ZERO
          dux_dgamma = ZERO

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + potential_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + potential_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

          ! derivatives of potential
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl


          ! derivative along x and along z
          if(PML_BOUNDARY_CONDITIONS .and. is_PML(ispec))then

          ispec_PML=spec_to_PML(ispec)

          PML_dux_dxl(i,j,ispec_PML) = dux_dxl
          PML_dux_dzl(i,j,ispec_PML)=dux_dzl

          dux_dxi = ZERO
          dux_dgamma = ZERO

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi &
            +(potential_acoustic(ibool(k,j,ispec))+deltat*potential_dot_acoustic(ibool(k,j,ispec)))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma &
            +(potential_acoustic(ibool(i,k,ispec))+deltat*potential_dot_acoustic(ibool(i,k,ispec)))*hprime_zz(j,k)
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

          ! derivatives of potential
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          PML_dux_dxl_new(i,j,ispec_PML) = dux_dxl
          PML_dux_dzl_new(i,j,ispec_PML) = dux_dzl
          endif


             if(PML_BOUNDARY_CONDITIONS .and. is_PML(ispec))then
               ispec_PML=spec_to_PML(ispec)
               iPML=ibool_PML(i,j,ispec)
               if ( (which_PML_elem(ILEFT,ispec) .OR. which_PML_elem(IRIGHT,ispec)) .and. &
                       .not. (which_PML_elem(ITOP,ispec) .OR. which_PML_elem(IBOTTOM,ispec)) ) then

                 bb = 0.d0 - (d_x_store(iPML) / k_x_store(iPML) + alpha_x_store(iPML))
                 if(abs(alpha_x_store(iPML)) > 1.d-3) then
                 aa = - d_x_store(iPML) / &
                        (k_x_store(iPML) * d_x_store(iPML) + k_x_store(iPML)**2 * alpha_x_store(iPML))
                 else
                 aa = - 1.d0 / k_x_store(iPML)
                 endif

                 coef0_x = exp(bb * deltat)
                 coef1_x = (1.d0 - exp(bb * deltat / 2.d0)) * aa
                 coef2_x = (1.d0 - exp(bb* deltat / 2.d0)) * exp(bb * deltat / 2.d0) * aa


                 rmemory_acoustic_dux_dx(1,i,j,ispec_PML) = 0.d0
                 rmemory_acoustic_dux_dx(2,i,j,ispec_PML) = coef0_x*rmemory_acoustic_dux_dx(2,i,j,ispec_PML) &
                 + PML_dux_dxl_new(i,j,ispec_PML) * coef1_x + PML_dux_dxl(i,j,ispec_PML) * coef2_x

                  if(abs(alpha_x_store(iPML)) > 1.d-3) then
                   bb = alpha_x_store(iPML)
                   coef0_x = exp(- bb * deltat)
                   coef1_x = (1 - exp(- bb * deltat / 2)) * d_x_store(iPML)/bb
                   coef2_x = (1 - exp(- bb * deltat / 2)) * exp(- bb * deltat / 2) * &
                             d_x_store(iPML)/bb
                 else
                   bb = alpha_x_store(iPML)
                  coef0_x = exp(- bb * deltat)
                  coef1_x = deltat / 2.0d0 * d_x_store(iPML)
                  coef2_x = deltat / 2.0d0 * d_x_store(iPML)
                  endif

                  rmemory_acoustic_dux_dz(1,i,j,ispec_PML) = coef0_x * rmemory_acoustic_dux_dz(1,i,j,ispec_PML) &
                  + PML_dux_dzl_new(i,j,ispec_PML) *coef1_x + PML_dux_dzl(i,j,ispec_PML) * coef2_x
                  rmemory_acoustic_dux_dz(2,i,j,ispec_PML) = 0.d0

                 dux_dxl = PML_dux_dxl(i,j,ispec_PML)  + rmemory_acoustic_dux_dx(1,i,j,ispec_PML) &
                                                       + rmemory_acoustic_dux_dx(2,i,j,ispec_PML)
                 dux_dzl = PML_dux_dzl(i,j,ispec_PML)  + rmemory_acoustic_dux_dz(1,i,j,ispec_PML) &
                                                       + rmemory_acoustic_dux_dz(2,i,j,ispec_PML)

                 elseif ( (which_PML_elem(ILEFT,ispec) .OR. which_PML_elem(IRIGHT,ispec)) .and. &
                       (which_PML_elem(ITOP,ispec) .OR. which_PML_elem(IBOTTOM,ispec)) ) then
                        ispec_PML=spec_to_PML(ispec)
                        iPML=ibool_PML(i,j,ispec)

                         if(abs(alpha_z_store(iPML)) .lt. 1.d-3)then
                            bb = alpha_z_store(iPML)
                            if(abs(k_x_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))+d_x_store(iPML)) .lt. 1.d-3)then
                            aa=d_z_store(iPML)/k_x_store(iPML)
                            else
                            aa = d_z_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))/&
                                     (k_x_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))+d_x_store(iPML))
                            endif
                            coef0_x = exp(- bb * deltat)
                            coef1_x = deltat / 2.0d0 * aa
                            coef2_x = deltat / 2.0d0 * aa
                         else
                            if(abs(k_x_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))+d_x_store(iPML)) .lt. 1.d-3)then
                            aa=d_z_store(iPML)/k_x_store(iPML)
                            else
                            aa = d_z_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))/&
                                     (k_x_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))+d_x_store(iPML))
                            endif
                            bb = alpha_z_store(iPML)
                            coef0_x = exp(- bb * deltat)
                            coef1_x = (1 - exp(- bb * deltat / 2)) * aa/bb
                            coef2_x = (1 - exp(- bb * deltat / 2)) * exp(- bb * deltat / 2) * aa /bb
                         endif


                  rmemory_acoustic_dux_dx(1,i,j,ispec_PML) = coef0_x*rmemory_acoustic_dux_dx(1,i,j,ispec_PML) &
                  + PML_dux_dxl_new(i,j,ispec_PML) * coef1_x + PML_dux_dxl(i,j,ispec_PML) * coef2_x


                        if(abs(d_x_store(iPML)*&
                                       (k_x_store(iPML)*k_z_store(iPML)*(alpha_z_store(iPML)-alpha_x_store(iPML))+&
                                        k_x_store(iPML)*d_z_store(iPML)-k_z_store(iPML)*d_x_store(iPML))/&
                                       (k_x_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))+d_x_store(iPML))/&
                                       k_x_store(iPML)**2) .lt. 1.d-3 .and. &
                        abs(k_x_store(iPML)*alpha_x_store(iPML)+d_x_store(iPML)) .lt. 1.d-3)then
                            coef0_z = 1.d0
                            coef1_z = 0.d0
                            coef2_z = 0.d0

                  rmemory_acoustic_dux_dx(2,i,j,ispec_PML) = coef0_z*rmemory_acoustic_dux_dx(2,i,j,ispec_PML) &
                  + PML_dux_dxl_new(i,j,ispec_PML) * coef1_z + PML_dux_dxl(i,j,ispec_PML) * coef2_z


                        elseif(abs(k_x_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))+d_x_store(iPML)) .lt. 1.d-3 .and. &
                        abs(k_x_store(iPML)*alpha_x_store(iPML)+d_x_store(iPML)) .lt. 1.d-3)then
                            bb=(k_x_store(iPML)*alpha_x_store(iPML)+d_x_store(iPML))/k_x_store(iPML)
                            aa=d_z_store(iPML)/k_x_store(iPML)
                            coef0_z = 1.d0
                            coef1_z = 0.5d0 *deltat
                            coef2_z = 0.5d0 *deltat

                  rmemory_acoustic_dux_dx(2,i,j,ispec_PML) = coef0_z*rmemory_acoustic_dux_dx(2,i,j,ispec_PML) &
                  + PML_dux_dxl_new(i,j,ispec_PML) * coef1_z*aa + PML_dux_dxl(i,j,ispec_PML) * coef2_z*aa


                        else
                            bb=(k_x_store(iPML)*alpha_x_store(iPML)+d_x_store(iPML))/k_x_store(iPML)
                            aa=d_x_store(iPML)*&
                                       (k_x_store(iPML)*k_z_store(iPML)*(alpha_z_store(iPML)-alpha_x_store(iPML))+&
                                        k_x_store(iPML)*d_z_store(iPML)-k_z_store(iPML)*d_x_store(iPML))/&
                                       (k_x_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))+d_x_store(iPML))/&
                                       k_x_store(iPML)**2
                            coef0_z = exp(- bb * deltat)
                            coef1_z = (1 - exp(- bb * deltat / 2)) * aa/bb
                            coef2_z = (1 - exp(- bb * deltat / 2)) * exp(- bb * deltat / 2) * aa /bb

                  rmemory_acoustic_dux_dx(2,i,j,ispec_PML) = coef0_z*rmemory_acoustic_dux_dx(2,i,j,ispec_PML) &
                  + PML_dux_dxl_new(i,j,ispec_PML) * coef1_z + PML_dux_dxl(i,j,ispec_PML) * coef2_z


                        endif

                         if(abs(alpha_x_store(iPML)) .lt. 1.d-3)then
                            bb = alpha_x_store(iPML)
                            if(abs(k_z_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))-d_z_store(iPML)) .lt. 1.d-3)then
                            aa=d_x_store(iPML)/k_z_store(iPML)
                            else
                            aa = d_x_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))/&
                                     (k_z_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))-d_z_store(iPML))
                            endif
                            coef0_x = exp(- bb * deltat)
                            coef1_x = deltat / 2.0d0 * aa
                            coef2_x = deltat / 2.0d0 * aa
                         else
                            if(abs(k_z_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))-d_z_store(iPML)) .lt. 1.d-3)then
                            aa=d_x_store(iPML)/k_z_store(iPML)
                            else
                            aa = d_x_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))/&
                                     (k_z_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))-d_z_store(iPML))
                            endif
                            bb = alpha_x_store(iPML)
                            coef0_x = exp(- bb * deltat)
                            coef1_x = (1 - exp(- bb * deltat / 2)) * aa/bb
                            coef2_x = (1 - exp(- bb * deltat / 2)) * exp(- bb * deltat / 2) * aa /bb
                         endif

                  rmemory_acoustic_dux_dz(1,i,j,ispec_PML) = coef0_x * rmemory_acoustic_dux_dz(1,i,j,ispec_PML) &
                  + PML_dux_dzl_new(i,j,ispec_PML) *coef1_x + PML_dux_dzl(i,j,ispec_PML) * coef2_x


                        if(abs(d_z_store(iPML)*&
                                       (k_x_store(iPML)*k_z_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))+&
                                        k_z_store(iPML)*d_x_store(iPML)-k_x_store(iPML)*d_z_store(iPML))/&
                                       (k_z_store(iPML)*(alpha_z_store(iPML)-alpha_x_store(iPML))+d_z_store(iPML))/&
                                       k_z_store(iPML)**2) .lt. 1.d-3 .and. &
                        abs(k_z_store(iPML)*alpha_z_store(iPML)+d_z_store(iPML)) .lt. 1.d-3)then
                            coef0_z = 1.d0
                            coef1_z = 0.d0
                            coef2_z = 0.d0

                  rmemory_acoustic_dux_dz(2,i,j,ispec_PML) = coef0_z * rmemory_acoustic_dux_dz(2,i,j,ispec_PML) &
                  + PML_dux_dzl_new(i,j,ispec_PML) *coef1_z + PML_dux_dzl(i,j,ispec_PML) * coef2_z

                        elseif(abs(k_z_store(iPML)*(alpha_z_store(iPML)-alpha_x_store(iPML))+d_z_store(iPML)) .lt. 1.d-3 .and. &
                        abs(k_z_store(iPML)*alpha_z_store(iPML)+d_z_store(iPML)) .lt. 1.d-3)then
                            bb=(k_z_store(iPML)*alpha_z_store(iPML)+d_z_store(iPML))/k_z_store(iPML)
                            aa=d_x_store(iPML)/k_z_store(iPML)
                            coef0_z = 1.d0
                            coef1_z = 0.5d0 *deltat
                            coef2_z = 0.5d0 *deltat

                  rmemory_acoustic_dux_dz(2,i,j,ispec_PML) = coef0_z * rmemory_acoustic_dux_dz(2,i,j,ispec_PML) &
                  + PML_dux_dzl_new(i,j,ispec_PML) *coef1_z*aa + PML_dux_dzl(i,j,ispec_PML) * coef2_z*aa

                        else
                            bb=(k_z_store(iPML)*alpha_z_store(iPML)+d_z_store(iPML))/k_z_store(iPML)
                            aa=d_z_store(iPML)*&
                                       (k_x_store(iPML)*k_z_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))+&
                                        k_z_store(iPML)*d_x_store(iPML)-k_x_store(iPML)*d_z_store(iPML))/&
                                       (k_z_store(iPML)*(alpha_z_store(iPML)-alpha_x_store(iPML))+d_z_store(iPML))/&
                                       k_z_store(iPML)**2
                            coef0_z = exp(- bb * deltat)
                            coef1_z = (1 - exp(- bb * deltat / 2)) * aa/bb
                            coef2_z = (1 - exp(- bb * deltat / 2)) * exp(- bb * deltat / 2) * aa /bb

                  rmemory_acoustic_dux_dz(2,i,j,ispec_PML) = coef0_z * rmemory_acoustic_dux_dz(2,i,j,ispec_PML) &
                  + PML_dux_dzl_new(i,j,ispec_PML) *coef1_z + PML_dux_dzl(i,j,ispec_PML) * coef2_z

                        endif

                 dux_dxl = PML_dux_dxl(i,j,ispec_PML)  + rmemory_acoustic_dux_dx(1,i,j,ispec_PML) +&
                           rmemory_acoustic_dux_dx(2,i,j,ispec_PML)
                 dux_dzl = PML_dux_dzl(i,j,ispec_PML)  + rmemory_acoustic_dux_dz(1,i,j,ispec_PML) + &
                           rmemory_acoustic_dux_dz(2,i,j,ispec_PML)

               else

                  if(abs(alpha_z_store(iPML)) > 1.d-3) then
                   bb = alpha_z_store(iPML)
                   coef0_x = exp(- bb * deltat)
                   coef1_x = (1 - exp(- bb * deltat / 2)) * d_z_store(iPML)/bb
                   coef2_x = (1 - exp(- bb * deltat / 2)) * exp(- bb * deltat / 2) * &
                             d_z_store(iPML) / bb
                 else
                   bb = alpha_z_store(iPML)
                  coef0_x = exp(- bb * deltat)
                  coef1_x = deltat / 2.0d0 * d_z_store(iPML)
                  coef2_x = deltat / 2.0d0 * d_z_store(iPML)
                  endif


                  rmemory_acoustic_dux_dx(1,i,j,ispec_PML) = coef0_x*rmemory_acoustic_dux_dx(1,i,j,ispec_PML) &
                  + PML_dux_dxl_new(i,j,ispec_PML) * coef1_x + PML_dux_dxl(i,j,ispec_PML) * coef2_x
                  rmemory_acoustic_dux_dx(2,i,j,ispec_PML) = 0.d0

                 bb = 0.d0 - (d_z_store(iPML) / k_z_store(iPML) + alpha_z_store(iPML))
                 if(abs(alpha_z_store(iPML)) > 1.d-3) then
                 aa = - d_z_store(iPML) * k_x_store(iPML) / &
                        (k_z_store(iPML) * d_z_store(iPML) + k_z_store(iPML)**2 * alpha_z_store(iPML))
                 else
                 aa = - k_x_store(iPML) / k_z_store(iPML)
                 endif
                 coef0_x = exp(bb * deltat)
                 coef1_x = (1.d0 - exp(bb * deltat / 2.d0)) * aa
                 coef2_x = (1.d0 - exp(bb* deltat / 2.d0)) * exp(bb * deltat / 2.d0) * aa

                  rmemory_acoustic_dux_dz(1,i,j,ispec_PML) = 0.d0
                  rmemory_acoustic_dux_dz(2,i,j,ispec_PML) = coef0_x * rmemory_acoustic_dux_dz(2,i,j,ispec_PML) &
                  + PML_dux_dzl_new(i,j,ispec_PML) *coef1_x + PML_dux_dzl(i,j,ispec_PML) * coef2_x

                 dux_dxl = dux_dxl  + rmemory_acoustic_dux_dx(1,i,j,ispec_PML) + rmemory_acoustic_dux_dx(2,i,j,ispec_PML)
                 dux_dzl = dux_dzl  + rmemory_acoustic_dux_dz(1,i,j,ispec_PML) + rmemory_acoustic_dux_dz(2,i,j,ispec_PML)

               endif
             endif

          jacobianl = jacobian(i,j,ispec)

          ! if external density model
          if(assign_external_model) rhol = rhoext(i,j,ispec)

          ! for acoustic medium
          ! also add GLL integration weights
          tempx1(i,j) = wzgll(j)*jacobianl*(xixl*dux_dxl + xizl*dux_dzl) / rhol
          tempx2(i,j) = wxgll(i)*jacobianl*(gammaxl*dux_dxl + gammazl*dux_dzl) / rhol
        enddo
      enddo


        ! first double loop over GLL points to compute and store gradients
        do j = 1,NGLLZ
           do i = 1,NGLLX
            if(assign_external_model) then
              rhol = rhoext(i,j,ispec)
            else
!              rhol = density(1,kmato(ispec))
        lambdal_relaxed = poroelastcoef(1,1,kmato(ispec))
        mul_relaxed = poroelastcoef(2,1,kmato(ispec))
        kappal  = lambdal_relaxed + TWO*mul_relaxed/3._CUSTOM_REAL
        rhol = density(1,kmato(ispec))
            endif

             if(is_PML(ispec) .and. PML_BOUNDARY_CONDITIONS)then
                        iPML=ibool_PML(i,j,ispec)
                        iglob=ibool(i,j,ispec)

               if ( (which_PML_elem(ILEFT,ispec) .OR. which_PML_elem(IRIGHT,ispec)) .and. &
                       .not. (which_PML_elem(ITOP,ispec) .OR. which_PML_elem(IBOTTOM,ispec)) ) then
                        ispec_PML=spec_to_PML(ispec)
                        iPML=ibool_PML(i,j,ispec)
                        iglob=ibool(i,j,ispec)
                        bb = alpha_x_store(iPML)
                        coef0_x = exp(- bb * deltat)
                        coef1_x = (1 - exp(- bb * deltat / 2)) * (d_x_store(iPML)*(bb))
                        coef2_x = (1 - exp(- bb * deltat / 2)) * exp(- bb * deltat / 2) * &
                                  (d_x_store(iPML)*(bb))

                        rmemory_potential_acoustic(1,i,j,ispec_PML)=coef0_x * rmemory_potential_acoustic(1,i,j,ispec_PML) &
                         + (potential_acoustic(iglob)+deltat*potential_dot_acoustic(iglob)) * coef1_x &
                         + potential_acoustic(iglob) * coef2_x

                        rmemory_potential_acoustic(2,i,j,ispec_PML)=0.0

               elseif ((which_PML_elem(ILEFT,ispec) .OR. which_PML_elem(IRIGHT,ispec)) .and. &
                       (which_PML_elem(ITOP,ispec) .OR. which_PML_elem(IBOTTOM,ispec)) ) then
                        ispec_PML=spec_to_PML(ispec)
                        iPML=ibool_PML(i,j,ispec)

                        if(abs(alpha_x_store(iPML)-alpha_z_store(iPML)).lt. 1.d-3)then
                          bb = alpha_x_store(iPML)
                          aa = (k_z_store(iPML)*d_x_store(iPML)*(bb))
                          coef0_x = exp(- bb * deltat)
                          coef1_x = (1 - exp(- bb * deltat / 2)) * aa
                          coef2_x = (1 - exp(- bb * deltat / 2)) * exp(- bb * deltat / 2) * aa

                         rmemory_potential_acoustic(1,i,j,ispec_PML)=&
                            coef0_x * rmemory_potential_acoustic(1,i,j,ispec_PML) &
                          + (potential_acoustic(iglob)+deltat*potential_dot_acoustic(iglob)) * coef1_x &
                          + potential_acoustic(iglob) * coef2_x
                        else
                          bb = alpha_x_store(iPML)
                          aa = bb*d_x_store(iPML)*(k_z_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))-d_z_store(iPML))/&
                               (alpha_x_store(iPML)-alpha_z_store(iPML))
                          coef0_x = exp(- bb * deltat)
                          coef1_x = (1 - exp(- bb * deltat / 2)) * aa
                          coef2_x = (1 - exp(- bb * deltat / 2)) * exp(- bb * deltat / 2) * aa

                         rmemory_potential_acoustic(1,i,j,ispec_PML)=&
                            coef0_x * rmemory_potential_acoustic(1,i,j,ispec_PML) &
                          + (potential_acoustic(iglob)+deltat*potential_dot_acoustic(iglob)) * coef1_x &
                          + potential_acoustic(iglob) * coef2_x

                        endif

                        if(abs(alpha_x_store(iPML)-alpha_z_store(iPML)).lt. 1.d-3)then
                          bb = alpha_z_store(iPML)
                          aa = k_x_store(iPML)*d_z_store(iPML)*(bb)
                          coef0_x = exp(- bb * deltat)
                          coef1_x = (1 - exp(- bb * deltat / 2)) * aa
                          coef2_x = (1 - exp(- bb * deltat / 2)) * exp(- bb * deltat / 2) * aa

                        rmemory_potential_acoustic(2,i,j,ispec_PML)=&
                           coef0_x * rmemory_potential_acoustic(2,i,j,ispec_PML) &
                         + (potential_acoustic(iglob)+deltat*potential_dot_acoustic(iglob)) * coef1_x &
                         + potential_acoustic(iglob) * coef2_x

                        else
                          bb = alpha_z_store(iPML)
                          aa = bb*d_z_store(iPML)*(k_x_store(iPML)*(alpha_x_store(iPML)-alpha_z_store(iPML))+d_x_store(iPML))/&
                               (alpha_x_store(iPML)-alpha_z_store(iPML))
                          coef0_x = exp(- bb * deltat)
                          coef1_x = (1 - exp(- bb * deltat / 2)) * aa
                          coef2_x = (1 - exp(- bb * deltat / 2)) * exp(- bb * deltat / 2) * aa

                         rmemory_potential_acoustic(2,i,j,ispec_PML)=&
                            coef0_x * rmemory_potential_acoustic(2,i,j,ispec_PML) &
                          + (potential_acoustic(iglob)+deltat*potential_dot_acoustic(iglob)) * coef1_x &
                          + potential_acoustic(iglob) * coef2_x
                        endif

               else

                        ispec_PML=spec_to_PML(ispec)
                        iPML=ibool_PML(i,j,ispec)
                        iglob=ibool(i,j,ispec)
                        bb = alpha_z_store(iPML)
                        coef0_x = exp(- bb * deltat)
                        coef1_x = (1 - exp(- bb * deltat / 2)) * (d_z_store(iPML)*bb)
                        coef2_x = (1 - exp(- bb * deltat / 2)) * exp(- bb * deltat / 2) * &
                                  (d_z_store(iPML)* bb )


                        rmemory_potential_acoustic(1,i,j,ispec_PML)=0.d0
                        rmemory_potential_acoustic(2,i,j,ispec_PML)=coef0_x * rmemory_potential_acoustic(2,i,j,ispec_PML) &
                         + (potential_acoustic(iglob)+deltat*potential_dot_acoustic(iglob)) * coef1_x &
                         + potential_acoustic(iglob) * coef2_x

               endif


               if ( (which_PML_elem(ILEFT,ispec) .OR. which_PML_elem(IRIGHT,ispec)) .and. &
                       .not. (which_PML_elem(ITOP,ispec) .OR. which_PML_elem(IBOTTOM,ispec)) ) then
                        ispec_PML=spec_to_PML(ispec)
                        iPML=ibool_PML(i,j,ispec)
                        iglob=ibool(i,j,ispec)
                    potential_dot_dot_acoustic_PML(i,j,ispec_PML)= wxgll(i)*wzgll(j)/ kappal*jacobian(i,j,ispec) * &
                      (d_x_store(iPML)*potential_dot_acoustic(iglob)+ &
                      rmemory_potential_acoustic(1,i,j,ispec_PML)+rmemory_potential_acoustic(2,i,j,ispec_PML)&
                      -d_x_store(iPML)*alpha_x_store(iPML)*potential_acoustic(iglob))


                elseif ((which_PML_elem(ILEFT,ispec) .OR. which_PML_elem(IRIGHT,ispec)) .and. &
                       (which_PML_elem(ITOP,ispec) .OR. which_PML_elem(IBOTTOM,ispec))) then
                        ispec_PML=spec_to_PML(ispec)
                        iPML=ibool_PML(i,j,ispec)

                    potential_dot_dot_acoustic_PML(i,j,ispec_PML)= wxgll(i)*wzgll(j)/ kappal*jacobian(i,j,ispec) * &
                      ((d_x_store(iPML)*k_z_store(iPML)+d_z_store(iPML)*k_x_store(iPML))*potential_dot_acoustic(iglob)+ &
                      rmemory_potential_acoustic(1,i,j,ispec_PML)+rmemory_potential_acoustic(2,i,j,ispec_PML)&
                      -k_z_store(iPML)*d_x_store(iPML)*alpha_x_store(iPML)*potential_acoustic(iglob)-&
                       k_x_store(iPML)*d_z_store(iPML)*alpha_z_store(iPML)*potential_acoustic(iglob)+&
                       d_z_store(iPML)*d_x_store(iPML)*potential_acoustic(iglob))


               else
                        ispec_PML=spec_to_PML(ispec)
                        iPML=ibool_PML(i,j,ispec)
                        iglob=ibool(i,j,ispec)
                    potential_dot_dot_acoustic_PML(i,j,ispec_PML)= wxgll(i)*wzgll(j)/ kappal*jacobian(i,j,ispec) * &
                      (d_z_store(iPML)*potential_dot_acoustic(iglob)+ &
                      rmemory_potential_acoustic(1,i,j,ispec_PML)+rmemory_potential_acoustic(2,i,j,ispec_PML)&
                      -d_z_store(iPML)*alpha_z_store(iPML)*potential_acoustic(iglob))
               endif
             endif

           enddo
        enddo



!
! second double-loop over GLL to compute all the terms
!
      do j = 1,NGLLZ
        do i = 1,NGLLX

          iglob = ibool(i,j,ispec)

          ! along x direction and z direction
          ! and assemble the contributions
          do k = 1,NGLLX
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
                           (tempx1(k,j)*hprimewgll_xx(k,i) + tempx2(i,k)*hprimewgll_zz(k,j))
          enddo

            if(is_PML(ispec) .and. PML_BOUNDARY_CONDITIONS)then
                if ((which_PML_elem(ILEFT,ispec) .OR. which_PML_elem(IRIGHT,ispec)) .and. &
                       .not. (which_PML_elem(ITOP,ispec) .OR. which_PML_elem(IBOTTOM,ispec))) then
                        ispec_PML=spec_to_PML(ispec)
                        iPML=ibool_PML(i,j,ispec)
                      potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                                         - potential_dot_dot_acoustic_PML(i,j,ispec_PML)

                 elseif ((which_PML_elem(ILEFT,ispec) .OR. which_PML_elem(IRIGHT,ispec)) .and. &
                       (which_PML_elem(ITOP,ispec) .OR. which_PML_elem(IBOTTOM,ispec))) then
                        ispec_PML=spec_to_PML(ispec)
                        iPML=ibool_PML(i,j,ispec)
                      potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                                          + potential_dot_dot_acoustic_PML(i,j,ispec_PML)

                      potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                                          - potential_dot_dot_acoustic_PML(i,j,ispec_PML)

               else
                      ispec_PML=spec_to_PML(ispec)
                      iPML=ibool_PML(i,j,ispec)
                      potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                                         - potential_dot_dot_acoustic_PML(i,j,ispec_PML)

               endif
             endif

        enddo ! second loop over the GLL points
      enddo

    endif ! end of test if acoustic element

  enddo ! end of loop over all spectral elements

 if(PML_BOUNDARY_CONDITIONS .and. anyabs) then

! we have to put Dirichlet on the boundary of the PML
     do ispecabs = 1,nelemabs

       ispec = numabs(ispecabs)

       if (is_PML(ispec)) then
      ispec_PML=spec_to_PML(ispec)
!--- left absorbing boundary
      if(codeabs(ILEFT,ispecabs)) then
        i = 1
        do j = 1,NGLLZ
          iglob = ibool(i,j,ispec)
          potential_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
       enddo
      endif
!--- right absorbing boundary
      if(codeabs(IRIGHT,ispecabs)) then
        i = NGLLX
        do j = 1,NGLLZ
          iglob = ibool(i,j,ispec)
          potential_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
        enddo

      endif
!--- bottom absorbing boundary
      if(codeabs(IBOTTOM,ispecabs)) then
        j = 1
! exclude corners to make sure there is no contradiction on the normal
        ibegin = 1
        iend = NGLLX
        do i = ibegin,iend
          iglob = ibool(i,j,ispec)
          potential_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
        enddo
      endif
!--- top absorbing boundary
      if(codeabs(ITOP,ispecabs)) then
        j = NGLLZ
! exclude corners to make sure there is no contradiction on the normal
        ibegin = 1
        iend = NGLLX
        do i = ibegin,iend
          iglob = ibool(i,j,ispec)
          potential_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
        enddo
      endif  !  end of top absorbing boundary
      endif ! end of is_PML
    enddo ! end specabs loop
  endif  ! end of absorbing boundaries


!
!--- absorbing boundaries
!
  if(.not. PML_BOUNDARY_CONDITIONS .and. anyabs) then

    do ispecabs=1,nelemabs

      ispec = numabs(ispecabs)

      ! Sommerfeld condition if acoustic
      if(.not. elastic(ispec) .and. .not. poroelastic(ispec)) then

        ! get elastic parameters of current spectral element
        lambdal_relaxed = poroelastcoef(1,1,kmato(ispec))
        mul_relaxed = poroelastcoef(2,1,kmato(ispec))
        kappal  = lambdal_relaxed + TWO*mul_relaxed/3._CUSTOM_REAL
        rhol = density(1,kmato(ispec))

        cpl = sqrt(kappal/rhol)

        !--- left absorbing boundary
        if(codeabs(ILEFT,ispecabs)) then
          i = 1
          jbegin = jbegin_left(ispecabs)
          jend = jend_left(ispecabs)
          do j = jbegin,jend
            iglob = ibool(i,j,ispec)
            ! external velocity model
            if(assign_external_model) then
              cpl = vpext(i,j,ispec)
              rhol = rhoext(i,j,ispec)
            endif
            xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
            zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            weight = jacobian1D * wzgll(j)

            if( SIMULATION_TYPE == 1 ) then
              ! adds absorbing boundary contribution
              potential_dot_dot_acoustic(iglob) = &
                  potential_dot_dot_acoustic(iglob) &
                  - potential_dot_acoustic(iglob)*weight/cpl/rhol
            elseif(SIMULATION_TYPE == 2) then
              if(IS_BACKWARD_FIELD) then
                ! adds (previously) stored contribution
                potential_dot_dot_acoustic(iglob) = &
                    potential_dot_dot_acoustic(iglob) &
                    - b_absorb_acoustic_left(j,ib_left(ispecabs),NSTEP-it+1)
              else
                ! adds absorbing boundary contribution
                potential_dot_dot_acoustic(iglob) = &
                    potential_dot_dot_acoustic(iglob) &
                    - potential_dot_acoustic(iglob)*weight/cpl/rhol
              endif
            endif

            if(SAVE_FORWARD .and. SIMULATION_TYPE ==1) then
              ! saves contribution
              b_absorb_acoustic_left(j,ib_left(ispecabs),it) = &
                  potential_dot_acoustic(iglob)*weight/cpl/rhol
            endif

          enddo

        endif  !  end of left absorbing boundary

        !--- right absorbing boundary
        if(codeabs(IRIGHT,ispecabs)) then
          i = NGLLX
          jbegin = jbegin_right(ispecabs)
          jend = jend_right(ispecabs)
          do j = jbegin,jend
            iglob = ibool(i,j,ispec)
            ! external velocity model
            if(assign_external_model) then
              cpl = vpext(i,j,ispec)
              rhol = rhoext(i,j,ispec)
            endif
            xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
            zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            weight = jacobian1D * wzgll(j)

            if( SIMULATION_TYPE == 1 ) then
              ! adds absorbing boundary contribution
              potential_dot_dot_acoustic(iglob) = &
                  potential_dot_dot_acoustic(iglob) &
                  - potential_dot_acoustic(iglob)*weight/cpl/rhol
            elseif(SIMULATION_TYPE == 2) then
              if(IS_BACKWARD_FIELD) then
                ! adds (previously) stored contribution
                potential_dot_dot_acoustic(iglob) = &
                    potential_dot_dot_acoustic(iglob) &
                    - b_absorb_acoustic_right(j,ib_right(ispecabs),NSTEP-it+1)
              else
                ! adds absorbing boundary contribution
                potential_dot_dot_acoustic(iglob) = &
                    potential_dot_dot_acoustic(iglob) &
                    - potential_dot_acoustic(iglob)*weight/cpl/rhol
              endif
            endif

            if(SAVE_FORWARD .and. SIMULATION_TYPE ==1) then
              ! saves contribution
              b_absorb_acoustic_right(j,ib_right(ispecabs),it) = &
                    potential_dot_acoustic(iglob)*weight/cpl/rhol
            endif
          enddo
        endif  !  end of right absorbing boundary

        !--- bottom absorbing boundary
        if(codeabs(IBOTTOM,ispecabs)) then
          j = 1
          ibegin = ibegin_bottom(ispecabs)
          iend = iend_bottom(ispecabs)
          ! exclude corners to make sure there is no contradiction on the normal
          if(codeabs(ILEFT,ispecabs)) ibegin = 2
          if(codeabs(IRIGHT,ispecabs)) iend = NGLLX-1
          do i = ibegin,iend
            iglob = ibool(i,j,ispec)
            ! external velocity model
            if(assign_external_model) then
              cpl = vpext(i,j,ispec)
              rhol = rhoext(i,j,ispec)
            endif
            xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
            zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            weight = jacobian1D * wxgll(i)

            if( SIMULATION_TYPE == 1 ) then
              ! adds absorbing boundary contribution
              potential_dot_dot_acoustic(iglob) = &
                  potential_dot_dot_acoustic(iglob) &
                  - potential_dot_acoustic(iglob)*weight/cpl/rhol
            elseif(SIMULATION_TYPE == 2) then
              if(IS_BACKWARD_FIELD) then
                ! adds (previously) stored contribution
                potential_dot_dot_acoustic(iglob) = &
                    potential_dot_dot_acoustic(iglob) &
                    - b_absorb_acoustic_bottom(i,ib_bottom(ispecabs),NSTEP-it+1)
              else
                ! adds absorbing boundary contribution
                potential_dot_dot_acoustic(iglob) = &
                    potential_dot_dot_acoustic(iglob) &
                    - potential_dot_acoustic(iglob)*weight/cpl/rhol
              endif
            endif

            if(SAVE_FORWARD .and. SIMULATION_TYPE ==1) then
              ! saves contribution
              b_absorb_acoustic_bottom(i,ib_bottom(ispecabs),it) = &
                  potential_dot_acoustic(iglob)*weight/cpl/rhol
            endif
          enddo
        endif  !  end of bottom absorbing boundary

        !--- top absorbing boundary
        if(codeabs(ITOP,ispecabs)) then
          j = NGLLZ
          ibegin = ibegin_top(ispecabs)
          iend = iend_top(ispecabs)
          ! exclude corners to make sure there is no contradiction on the normal
          if(codeabs(ILEFT,ispecabs)) ibegin = 2
          if(codeabs(IRIGHT,ispecabs)) iend = NGLLX-1
          do i = ibegin,iend
            iglob = ibool(i,j,ispec)
            ! external velocity model
            if(assign_external_model) then
              cpl = vpext(i,j,ispec)
              rhol = rhoext(i,j,ispec)
            endif
            xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
            zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            weight = jacobian1D * wxgll(i)

            if( SIMULATION_TYPE == 1 ) then
              ! adds absorbing boundary contribution
              potential_dot_dot_acoustic(iglob) = &
                  potential_dot_dot_acoustic(iglob) &
                  - potential_dot_acoustic(iglob)*weight/cpl/rhol
            elseif(SIMULATION_TYPE == 2) then
              if(IS_BACKWARD_FIELD) then
                ! adds (previously) stored contribution
                potential_dot_dot_acoustic(iglob) = &
                    potential_dot_dot_acoustic(iglob) &
                    - b_absorb_acoustic_top(i,ib_top(ispecabs),NSTEP-it+1)
              else
                ! adds absorbing boundary contribution
                potential_dot_dot_acoustic(iglob) = &
                    potential_dot_dot_acoustic(iglob) &
                    - potential_dot_acoustic(iglob)*weight/cpl/rhol
              endif
            endif

            if(SAVE_FORWARD .and. SIMULATION_TYPE ==1) then
              ! saves contribution
              b_absorb_acoustic_top(i,ib_top(ispecabs),it) = &
                  potential_dot_acoustic(iglob)*weight/cpl/rhol
            endif
          enddo
        endif  !  end of top absorbing boundary

      endif ! acoustic ispec
    enddo
  endif  ! end of absorbing boundaries

  end subroutine compute_forces_acoustic

