
!========================================================================
!
!                  S P E C F E M 2 D  Version 7 . 0
!                  --------------------------------
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
               potential_acoustic,potential_acoustic_old,stage_time_scheme,i_stage, &
               density,poroelastcoef,xix,xiz,gammax,gammaz,jacobian, &
               vpext,rhoext,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll, &
               ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
               ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2, &
               SIMULATION_TYPE,SAVE_FORWARD,nspec_left,nspec_right,&
               nspec_bottom,nspec_top,ib_left,ib_right,ib_bottom,ib_top, &
               b_absorb_acoustic_left,b_absorb_acoustic_right, &
               b_absorb_acoustic_bottom,b_absorb_acoustic_top,IS_BACKWARD_FIELD,&
               is_PML,nspec_PML,spec_to_PML,region_CPML, &
               K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store,&
               rmemory_potential_acoustic,&
               rmemory_acoustic_dux_dx,rmemory_acoustic_dux_dz,&
               rmemory_potential_acoust_LDDRK,alpha_LDDRK,beta_LDDRK, &
               rmemory_acoustic_dux_dx_LDDRK,rmemory_acoustic_dux_dz_LDDRK,&
               deltat,PML_BOUNDARY_CONDITIONS,STACEY_BOUNDARY_CONDITIONS)


! compute forces for the acoustic elements

  implicit none

  include "constants.h"

  integer :: nglob,nspec,nelemabs,numat,it,NSTEP,SIMULATION_TYPE

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: kmato
  integer, dimension(nelemabs) :: numabs,ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
               ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2

  logical, dimension(nspec) :: elastic,poroelastic
  logical, dimension(4,nelemabs)  :: codeabs

  real(kind=CUSTOM_REAL), dimension(nglob) :: &
    potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic,potential_acoustic_old

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
  integer :: nspec_PML,ispec_PML
  integer, dimension(nspec) :: region_CPML
  logical, dimension(nspec) :: is_PML
  integer, dimension(nspec) :: spec_to_PML

  real(kind=CUSTOM_REAL), dimension(2,NGLLX,NGLLZ,nspec_PML) :: rmemory_potential_acoustic
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_PML,2) :: &
                          rmemory_acoustic_dux_dx,rmemory_acoustic_dux_dz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_PML) :: &
                          K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: potential_dot_dot_acoustic_PML
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: PML_dux_dxl,PML_dux_dzl,PML_dux_dxl_old,PML_dux_dzl_old
  real(kind=CUSTOM_REAL) :: kappa_x,kappa_z,d_x,d_z,alpha_x,alpha_z,beta_x,beta_z,time_n,time_nsub1,&                            
                            A5,A6,A7, bb_zx_1,bb_zx_2,coef0_zx_1,coef1_zx_1,coef2_zx_1,coef0_zx_2,coef1_zx_2,coef2_zx_2,&
                            A8,A9,A10,bb_xz_1,bb_xz_2,coef0_xz_1,coef1_xz_1,coef2_xz_1,coef0_xz_2,coef1_xz_2,coef2_xz_2,&
                            A0,A1,A2,A3,A4,bb_1,coef0_1,coef1_1,coef2_1,bb_2,coef0_2,coef1_2,coef2_2
  integer :: CPML_region_local,singularity_type_zx,singularity_type_xz,singularity_type
  double precision :: deltat

  logical :: PML_BOUNDARY_CONDITIONS,STACEY_BOUNDARY_CONDITIONS

!coefficients and memory variables when using CPML with LDDRK
  integer :: stage_time_scheme,i_stage
  real(kind=CUSTOM_REAL), dimension(Nstages) :: alpha_LDDRK,beta_LDDRK
  real(kind=CUSTOM_REAL), dimension(2,NGLLX,NGLLZ,nspec_PML) :: rmemory_potential_acoust_LDDRK
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_PML,2) :: rmemory_acoustic_dux_dx_LDDRK,rmemory_acoustic_dux_dz_LDDRK

  ifirstelem = 1
  ilastelem = nspec
  time_n = (it-1) * deltat
  time_nsub1 = (it-2) * deltat

  if( PML_BOUNDARY_CONDITIONS ) then
    potential_dot_dot_acoustic_PML = 0._CUSTOM_REAL
    PML_dux_dxl = 0._CUSTOM_REAL
    PML_dux_dzl = 0._CUSTOM_REAL
    PML_dux_dxl_old = 0._CUSTOM_REAL
    PML_dux_dzl_old = 0._CUSTOM_REAL
  endif

! loop over spectral elements
  do ispec = ifirstelem,ilastelem
!---
!--- acoustic spectral element
!---
    if(.not. elastic(ispec) .and. .not. poroelastic(ispec)) then
      if(CUSTOM_REAL == SIZE_REAL) then
        rhol = sngl(density(1,kmato(ispec)))
      else
        rhol = density(1,kmato(ispec))
      endif

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

          ! derivative along x and along zbb
          if(PML_BOUNDARY_CONDITIONS .and. is_PML(ispec) .and. nspec_PML > 0)then
            ispec_PML=spec_to_PML(ispec)
            CPML_region_local = region_CPML(ispec)
            kappa_x = K_x_store(i,j,ispec_PML)
            kappa_z = K_z_store(i,j,ispec_PML)
            d_x = d_x_store(i,j,ispec_PML)
            d_z = d_z_store(i,j,ispec_PML)
            alpha_x = alpha_x_store(i,j,ispec_PML)
            alpha_z = alpha_z_store(i,j,ispec_PML)
            beta_x = alpha_x + d_x / kappa_x
            beta_z = alpha_z + d_z / kappa_z

            PML_dux_dxl(i,j) = dux_dxl
            PML_dux_dzl(i,j) = dux_dzl

            dux_dxi = ZERO
            dux_dgamma = ZERO

            ! first double loop over GLL points to compute and store gradients
            ! we can merge the two loops because NGLLX == NGLLZ
            do k = 1,NGLLX
              dux_dxi = dux_dxi + potential_acoustic_old(ibool(k,j,ispec)) * hprime_xx(i,k)
              dux_dgamma = dux_dgamma + potential_acoustic_old(ibool(i,k,ispec)) * hprime_zz(j,k)
            enddo

            xixl = xix(i,j,ispec)
            xizl = xiz(i,j,ispec)
            gammaxl = gammax(i,j,ispec)
            gammazl = gammaz(i,j,ispec)
            ! derivatives of potential
            PML_dux_dxl_old(i,j) = dux_dxi*xixl + dux_dgamma*gammaxl
            PML_dux_dzl_old(i,j) = dux_dxi*xizl + dux_dgamma*gammazl

            call lik_parameter_computation(time_n,deltat,kappa_z,beta_z,alpha_z,kappa_x,beta_x,alpha_x,&
                                           CPML_region_local,31,A5,A6,A7,singularity_type_zx,bb_zx_1,bb_zx_2,&
                                           coef0_zx_1,coef1_zx_1,coef2_zx_1,coef0_zx_2,coef1_zx_2,coef2_zx_2)
            call lik_parameter_computation(time_n,deltat,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z,&
                                           CPML_region_local,13,A8,A9,A10,singularity_type_xz,bb_xz_1,bb_xz_2,&
                                           coef0_xz_1,coef1_xz_1,coef2_xz_1,coef0_xz_2,coef1_xz_2,coef2_xz_2)
            if(stage_time_scheme == 1) then
              rmemory_acoustic_dux_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_acoustic_dux_dx(i,j,ispec_PML,1) + &
                                                         coef1_zx_1 * PML_dux_dxl(i,j) + coef2_zx_1 * PML_dux_dxl_old(i,j)
              if(singularity_type_zx == 0)then

                rmemory_acoustic_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + &
                                                           coef1_zx_2 * PML_dux_dxl(i,j) + coef2_zx_2 * PML_dux_dxl_old(i,j)
              else
                rmemory_acoustic_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + &
                                                           coef1_zx_2 * time_n * PML_dux_dxl(i,j) + &
                                                           coef2_zx_2 * time_nsub1 * PML_dux_dxl_old(i,j)
              endif

              rmemory_acoustic_dux_dz(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_acoustic_dux_dz(i,j,ispec_PML,1) + &
                                                         coef1_xz_1 * PML_dux_dzl(i,j) + coef2_xz_1 * PML_dux_dzl_old(i,j)
              if(singularity_type_xz == 0)then
                rmemory_acoustic_dux_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + &
                                                           coef1_xz_2 * PML_dux_dzl(i,j) + coef2_xz_2 * PML_dux_dzl_old(i,j)
              else
                rmemory_acoustic_dux_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + &
                                                           coef1_xz_2 * time_n * PML_dux_dzl(i,j) + &
                                                           coef2_xz_2 * time_nsub1 * PML_dux_dzl_old(i,j)
              endif
            endif

            if(stage_time_scheme == 6) then
              rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,1) = &
                     alpha_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,1) + &
                     deltat * (-bb_zx_1 * rmemory_acoustic_dux_dx(i,j,ispec_PML,1) + PML_dux_dxl(i,j))
              rmemory_acoustic_dux_dx(i,j,ispec_PML,1) = rmemory_acoustic_dux_dx(i,j,ispec_PML,1) + &
                     beta_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,1)
              if(singularity_type_zx == 0)then
                rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2) = &
                       alpha_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2) + &
                       deltat * (-bb_zx_2 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + PML_dux_dxl(i,j))
                rmemory_acoustic_dux_dx(i,j,ispec_PML,2) = rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + &
                       beta_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2)
              else
                rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2) = &
                       alpha_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2) + &
                       deltat * (-bb_zx_2 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + PML_dux_dxl(i,j) * time_n)
                rmemory_acoustic_dux_dx(i,j,ispec_PML,2) = rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + &
                       beta_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2)
              endif

              rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,1) = &
                     alpha_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,1) + &
                     deltat * (-bb_xz_1 * rmemory_acoustic_dux_dz(i,j,ispec_PML,1) + PML_dux_dzl(i,j))
              rmemory_acoustic_dux_dz(i,j,ispec_PML,1) = rmemory_acoustic_dux_dz(i,j,ispec_PML,1) + &
                     beta_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,1)
              if(singularity_type_xz == 0)then
                rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2) = &
                       alpha_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2) + &
                       deltat * (-bb_xz_2 * rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + PML_dux_dzl(i,j))
                rmemory_acoustic_dux_dz(i,j,ispec_PML,2) = rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + &
                       beta_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2)
              else
                rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2) = &
                       alpha_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2) + &
                       deltat * (-bb_xz_2 * rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + PML_dux_dzl(i,j) * time_n)
                rmemory_acoustic_dux_dz(i,j,ispec_PML,2) = rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + &
                       beta_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2)
              endif

            endif

            dux_dxl = A5 * PML_dux_dxl(i,j) + A6 * rmemory_acoustic_dux_dx(i,j,ispec_PML,1) + &
                                              A7 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2)
            dux_dzl = A8 * PML_dux_dzl(i,j) + A9 *  rmemory_acoustic_dux_dz(i,j,ispec_PML,1) + &
                                              A10 * rmemory_acoustic_dux_dz(i,j,ispec_PML,2)
          endif

          jacobianl = jacobian(i,j,ispec)

          ! if external density model
          if(assign_external_model)then
            if(CUSTOM_REAL == SIZE_REAL) then
              rhol = sngl(rhoext(i,j,ispec))
            else
              rhol = rhoext(i,j,ispec)
            endif
          endif
          ! for acoustic medium
          ! also add GLL integration weights
          tempx1(i,j) = wzgll(j)*jacobianl*(xixl*dux_dxl + xizl*dux_dzl) / rhol
          tempx2(i,j) = wxgll(i)*jacobianl*(gammaxl*dux_dxl + gammazl*dux_dzl) / rhol
        enddo
      enddo

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          if(assign_external_model) then
            if(CUSTOM_REAL == SIZE_REAL) then
              rhol = sngl(rhoext(i,j,ispec))
            else
              rhol = rhoext(i,j,ispec)
            endif
            cpl = vpext(i,j,ispec)
            !assuming that in fluid(acoustic) part input cpl is defined by sqrt(kappal/rhol), &
            !which is not the same as in cpl input in elastic part
            kappal = rhol*cpl*cpl
          else
            if(CUSTOM_REAL == SIZE_REAL) then
              lambdal_relaxed = sngl(poroelastcoef(1,1,kmato(ispec)))
              mul_relaxed = sngl(poroelastcoef(2,1,kmato(ispec)))
            else
              lambdal_relaxed = poroelastcoef(1,1,kmato(ispec))
              mul_relaxed = poroelastcoef(2,1,kmato(ispec))
            endif
            kappal  = lambdal_relaxed + TWO*mul_relaxed/3._CUSTOM_REAL
            rhol = density(1,kmato(ispec))
          endif

          if(is_PML(ispec) .and. PML_BOUNDARY_CONDITIONS .and. nspec_PML > 0)then
            ispec_PML=spec_to_PML(ispec)
            CPML_region_local = region_CPML(ispec)
            kappa_x = K_x_store(i,j,ispec_PML)
            kappa_z = K_z_store(i,j,ispec_PML)
            d_x = d_x_store(i,j,ispec_PML)
            d_z = d_z_store(i,j,ispec_PML)
            alpha_x = alpha_x_store(i,j,ispec_PML)
            alpha_z = alpha_z_store(i,j,ispec_PML)
            beta_x = alpha_x + d_x / kappa_x
            beta_z = alpha_z + d_z / kappa_z
            call l_parameter_computation(time_n,deltat,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                         CPML_region_local,A0,A1,A2,A3,A4,singularity_type,&
                                         bb_1,coef0_1,coef1_1,coef2_1,bb_2,coef0_2,coef1_2,coef2_2)

            if(stage_time_scheme == 1) then
              rmemory_potential_acoustic(1,i,j,ispec_PML) = coef0_1 * rmemory_potential_acoustic(1,i,j,ispec_PML) + &
                     coef1_1 * potential_acoustic(iglob) + coef2_1 * potential_acoustic_old(iglob)

              if(singularity_type == 0)then
                rmemory_potential_acoustic(2,i,j,ispec_PML) = coef0_2 * rmemory_potential_acoustic(2,i,j,ispec_PML) + &
                       coef1_2 * potential_acoustic(iglob) + coef2_2 * potential_acoustic_old(iglob)
              else
                rmemory_potential_acoustic(2,i,j,ispec_PML) = coef0_2 * rmemory_potential_acoustic(2,i,j,ispec_PML) + &
                       coef1_2 * time_n * potential_acoustic(iglob) + coef2_2 * time_nsub1 * potential_acoustic_old(iglob)
              endif
            endif

            if(stage_time_scheme == 6) then
              rmemory_potential_acoust_LDDRK(1,i,j,ispec_PML) = &
                     alpha_LDDRK(i_stage) * rmemory_potential_acoust_LDDRK(1,i,j,ispec_PML) + &
                     deltat * (-bb_1 * rmemory_potential_acoustic(1,i,j,ispec_PML) + potential_acoustic(iglob))
              if(singularity_type == 0)then
                rmemory_potential_acoust_LDDRK(2,i,j,ispec_PML) = &
                       alpha_LDDRK(i_stage) * rmemory_potential_acoust_LDDRK(2,i,j,ispec_PML) + &
                       deltat * (-bb_2 * rmemory_potential_acoustic(2,i,j,ispec_PML) + potential_acoustic(iglob))
              else
                rmemory_potential_acoust_LDDRK(2,i,j,ispec_PML) = &
                       alpha_LDDRK(i_stage) * rmemory_potential_acoust_LDDRK(2,i,j,ispec_PML) + &
                       deltat * (-bb_2 * rmemory_potential_acoustic(2,i,j,ispec_PML) + potential_acoustic(iglob) * time_n)
              endif
            endif

            potential_dot_dot_acoustic_PML(i,j)= wxgll(i)*wzgll(j)/kappal*jacobian(i,j,ispec) * &
                     (A1 * potential_dot_acoustic(iglob) + A2 * potential_acoustic(iglob) + &
                      A3 * rmemory_potential_acoustic(1,i,j,ispec_PML) + A4 * rmemory_potential_acoustic(2,i,j,ispec_PML))
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
          if(PML_BOUNDARY_CONDITIONS .and. is_PML(ispec))then
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                              - potential_dot_dot_acoustic_PML(i,j)
          endif
        enddo ! second loop over the GLL points
      enddo
            
    endif ! end of test if acoustic element
  enddo ! end of loop over all spectral elements



!
!--- absorbing boundaries
!

! The outer boundary condition to use for PML elements in fluid layers is Neumann for the potential
! because we need Dirichlet conditions for the displacement vector, which means Neumann for the potential.
! Thus, there is nothing to enforce explicitly here.
! There is something to enforce explicitly only in the case of elastic elements, for which a Dirichlet
! condition is needed for the displacement vector, which is the vectorial unknown for these elements.

! for Stacey paraxial absorbing conditions (more precisely: Sommerfeld in the case of a fluid) we implement them here

  if(STACEY_BOUNDARY_CONDITIONS) then

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
        if(codeabs(IEDGE4,ispecabs)) then
          i = 1
          jbegin = ibegin_edge4(ispecabs)
          jend = iend_edge4(ispecabs)
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
            else if(SIMULATION_TYPE == 3) then
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
        if(codeabs(IEDGE2,ispecabs)) then
          i = NGLLX
          jbegin = ibegin_edge2(ispecabs)
          jend = iend_edge2(ispecabs)
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
            else if(SIMULATION_TYPE == 3) then
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
        if(codeabs(IEDGE1,ispecabs)) then
          j = 1
          ibegin = ibegin_edge1(ispecabs)
          iend = iend_edge1(ispecabs)
          ! exclude corners to make sure there is no contradiction on the normal
          if(codeabs(IEDGE4,ispecabs)) ibegin = 2
          if(codeabs(IEDGE2,ispecabs)) iend = NGLLX-1
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
            else if(SIMULATION_TYPE == 3) then
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
        if(codeabs(IEDGE3,ispecabs)) then
          j = NGLLZ
          ibegin = ibegin_edge3(ispecabs)
          iend = iend_edge3(ispecabs)
          ! exclude corners to make sure there is no contradiction on the normal
          if(codeabs(IEDGE4,ispecabs)) ibegin = 2
          if(codeabs(IEDGE2,ispecabs)) iend = NGLLX-1
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
            else if(SIMULATION_TYPE == 3) then
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!--- absorbing boundaries
!
 if( PML_BOUNDARY_CONDITIONS .and. anyabs) then

! we have to put Dirichlet on the boundary of the PML
   do ispecabs = 1,nelemabs

     ispec = numabs(ispecabs)

     if (is_PML(ispec)) then
      ispec_PML=spec_to_PML(ispec)
!--- left absorbing boundary
      if(codeabs(IEDGE4,ispecabs)) then
        i = 1
        do j = 1,NGLLZ
          iglob = ibool(i,j,ispec)
          potential_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
       enddo
      endif
!--- right absorbing boundary
      if(codeabs(IEDGE2,ispecabs)) then
        i = NGLLX
        do j = 1,NGLLZ
          iglob = ibool(i,j,ispec)
          potential_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
        enddo

      endif
!--- bottom absorbing boundary
      if(codeabs(IEDGE1,ispecabs)) then
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
      if(codeabs(IEDGE3,ispecabs)) then
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
  endif  ! end of absorbing boundaries PML_BOUNDARY_CONDITIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine compute_forces_acoustic





