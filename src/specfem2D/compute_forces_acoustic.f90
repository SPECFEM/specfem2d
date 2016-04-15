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
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================


  subroutine compute_forces_acoustic(potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic, &
                                     PML_BOUNDARY_CONDITIONS,potential_acoustic_old)


! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,CPML_X_ONLY,CPML_Z_ONLY,IRIGHT,ILEFT,IBOTTOM,ITOP, &
    ZERO,ONE,TWO,IEDGE1,IEDGE2,IEDGE3,IEDGE4, &
    ALPHA_LDDRK,BETA_LDDRK,C_LDDRK

  use specfem_par, only: nglob,nspec,it, &
                         assign_external_model,ibool,kmato,ispec_is_acoustic, &
                         stage_time_scheme,i_stage, &
                         density,poroelastcoef,xix,xiz,gammax,gammaz,jacobian, &
                         vpext,rhoext, &
                         hprime_xx,hprimewgll_xx, &
                         hprime_zz,hprimewgll_zz,wxgll,wzgll, &
                         AXISYM,coord, is_on_the_axis,hprimeBar_xx,hprimeBarwglj_xx,xiglj,wxglj, &
                         rmemory_potential_acoustic,&
                         rmemory_acoustic_dux_dx,rmemory_acoustic_dux_dz,&
                         rmemory_potential_acoustic_LDDRK, &
                         rmemory_acoustic_dux_dx_LDDRK,rmemory_acoustic_dux_dz_LDDRK,&
                         deltat

  ! PML arrays
  use specfem_par, only: nspec_PML,ispec_is_PML,spec_to_PML,region_CPML, &
                K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob) :: potential_dot_dot_acoustic,potential_dot_acoustic, &
                                              potential_acoustic

  logical,intent(in) :: PML_BOUNDARY_CONDITIONS
  real(kind=CUSTOM_REAL), dimension(nglob) :: potential_acoustic_old

  ! local parameters
  integer :: ispec,i,j,k,iglob
  integer :: ifirstelem,ilastelem

  ! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,dux_dxl,dux_dzl
  real(kind=CUSTOM_REAL) :: xxi

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  ! material properties of the acoustic medium
  real(kind=CUSTOM_REAL) :: mul_relaxed,lambdal_relaxed,kappal,cpl,rhol

  ! local PML parameters
  integer :: ispec_PML
  integer :: CPML_region_local,singularity_type_zx,singularity_type_xz,singularity_type
  double precision :: kappa_x,kappa_z,d_x,d_z,alpha_x,alpha_z,beta_x,beta_z,time_n,time_nsub1,&
                      A5,A6,A7, bb_zx_1,bb_zx_2,coef0_zx_1,coef1_zx_1,coef2_zx_1,coef0_zx_2,coef1_zx_2,coef2_zx_2,&
                      A8,A9,A10,bb_xz_1,bb_xz_2,coef0_xz_1,coef1_xz_1,coef2_xz_1,coef0_xz_2,coef1_xz_2,coef2_xz_2,&
                      A0,A1,A2,A3,A4,bb_1,coef0_1,coef1_1,coef2_1,bb_2,coef0_2,coef1_2,coef2_2

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: PML_dux_dxl,PML_dux_dzl,PML_dux_dxl_old,PML_dux_dzl_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: potential_dot_dot_acoustic_PML

  ifirstelem = 1
  ilastelem = nspec
  if (stage_time_scheme == 1) then
    time_n = (it-1) * deltat
    time_nsub1 = (it-2) * deltat
  else if (stage_time_scheme == 6) then
    time_n = (it-1) * deltat + C_LDDRK(i_stage) * deltat
  endif

  if (PML_BOUNDARY_CONDITIONS) then
    potential_dot_dot_acoustic_PML(:,:) = 0._CUSTOM_REAL
    PML_dux_dxl(:,:) = 0._CUSTOM_REAL
    PML_dux_dzl(:,:) = 0._CUSTOM_REAL
    PML_dux_dxl_old(:,:) = 0._CUSTOM_REAL
    PML_dux_dzl_old(:,:) = 0._CUSTOM_REAL
  endif

! loop over spectral elements
  do ispec = ifirstelem,ilastelem

    ! acoustic spectral element
    if (ispec_is_acoustic(ispec)) then

      rhol = density(1,kmato(ispec))

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX
          ! derivative along x and along z
          dux_dxi = 0._CUSTOM_REAL
          dux_dgamma = 0._CUSTOM_REAL

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + potential_acoustic(ibool(k,j,ispec)) * hprime_xx(i,k)
            dux_dgamma = dux_dgamma + potential_acoustic(ibool(i,k,ispec)) * hprime_zz(j,k)
          enddo

          ! AXISYM case overwrites dux_dxi
          if (AXISYM) then
            if (is_on_the_axis(ispec)) then
              dux_dxi = 0._CUSTOM_REAL
              do k = 1,NGLLX
                dux_dxi = dux_dxi + potential_acoustic(ibool(k,j,ispec)) * hprimeBar_xx(i,k)
              enddo
            endif
          endif

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

          ! derivatives of potential
          dux_dxl = dux_dxi * xixl + dux_dgamma * gammaxl
          dux_dzl = dux_dxi * xizl + dux_dgamma * gammazl

          ! AXISYM case
          if (AXISYM) then
            if (is_on_the_axis(ispec) .and. i == 1) then ! dchi/dr=rho * u_r=0 on the axis
              dux_dxl = ZERO
            endif
          endif

          ! derivative along x and along zbb
          if (PML_BOUNDARY_CONDITIONS) then
            if (ispec_is_PML(ispec) .and. nspec_PML > 0) then
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

              dux_dxi = 0._CUSTOM_REAL; dux_dgamma = 0._CUSTOM_REAL

              if (stage_time_scheme == 1) then
                do k = 1,NGLLX
                  if (AXISYM) then
                    if (is_on_the_axis(ispec)) then
                      dux_dxi = dux_dxi + potential_acoustic_old(ibool(k,j,ispec)) * hprimeBar_xx(i,k)
                    else
                      dux_dxi = dux_dxi + potential_acoustic_old(ibool(k,j,ispec)) * hprime_xx(i,k)
                    endif
                  else
                    dux_dxi = dux_dxi + potential_acoustic_old(ibool(k,j,ispec)) * hprime_xx(i,k)
                  endif
                  dux_dgamma = dux_dgamma + potential_acoustic_old(ibool(i,k,ispec)) * hprime_zz(j,k)
                enddo

                ! derivatives of potential
                PML_dux_dxl_old(i,j) = dux_dxi * xixl + dux_dgamma * gammaxl
                PML_dux_dzl_old(i,j) = dux_dxi * xizl + dux_dgamma * gammazl

                if (AXISYM .and. is_on_the_axis(ispec) .and. i == 1) then ! dchi/dr=rho * u_r=0 on the axis
                  PML_dux_dxl_old(i,j) = ZERO
                endif
              endif

              ! the subroutine of lik_parameter_computation is located at the end of compute_forces_viscoelastic.F90
              call lik_parameter_computation(time_n,deltat,kappa_z,beta_z,alpha_z,kappa_x,beta_x,alpha_x,&
                                             CPML_region_local,31,A5,A6,A7,singularity_type_zx,bb_zx_1,bb_zx_2,&
                                             coef0_zx_1,coef1_zx_1,coef2_zx_1,coef0_zx_2,coef1_zx_2,coef2_zx_2)
              call lik_parameter_computation(time_n,deltat,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z,&
                                             CPML_region_local,13,A8,A9,A10,singularity_type_xz,bb_xz_1,bb_xz_2,&
                                             coef0_xz_1,coef1_xz_1,coef2_xz_1,coef0_xz_2,coef1_xz_2,coef2_xz_2)

              if (stage_time_scheme == 1) then
                rmemory_acoustic_dux_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_acoustic_dux_dx(i,j,ispec_PML,1) + &
                                                           coef1_zx_1 * PML_dux_dxl(i,j) + coef2_zx_1 * PML_dux_dxl_old(i,j)
                if (singularity_type_zx == 0) then

                  rmemory_acoustic_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + &
                                                             coef1_zx_2 * PML_dux_dxl(i,j) + coef2_zx_2 * PML_dux_dxl_old(i,j)
                else
                  rmemory_acoustic_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + &
                                                             coef1_zx_2 * time_n * PML_dux_dxl(i,j) + &
                                                             coef2_zx_2 * time_nsub1 * PML_dux_dxl_old(i,j)
                endif

                rmemory_acoustic_dux_dz(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_acoustic_dux_dz(i,j,ispec_PML,1) + &
                                                           coef1_xz_1 * PML_dux_dzl(i,j) + coef2_xz_1 * PML_dux_dzl_old(i,j)
                if (singularity_type_xz == 0) then
                  rmemory_acoustic_dux_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + &
                                                             coef1_xz_2 * PML_dux_dzl(i,j) + coef2_xz_2 * PML_dux_dzl_old(i,j)
                else
                  rmemory_acoustic_dux_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + &
                                                             coef1_xz_2 * time_n * PML_dux_dzl(i,j) + &
                                                             coef2_xz_2 * time_nsub1 * PML_dux_dzl_old(i,j)
                endif
              endif

              if (stage_time_scheme == 6) then
                rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,1) = &
                       ALPHA_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,1) + &
                       deltat * (-bb_zx_1 * rmemory_acoustic_dux_dx(i,j,ispec_PML,1) + PML_dux_dxl(i,j))
                rmemory_acoustic_dux_dx(i,j,ispec_PML,1) = rmemory_acoustic_dux_dx(i,j,ispec_PML,1) + &
                       BETA_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,1)
                if (singularity_type_zx == 0) then
                  rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2) = &
                         ALPHA_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2) + &
                         deltat * (-bb_zx_2 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + PML_dux_dxl(i,j))
                  rmemory_acoustic_dux_dx(i,j,ispec_PML,2) = rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + &
                         BETA_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2)
                else
                  rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2) = &
                         ALPHA_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2) + &
                         deltat * (-bb_zx_2 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + PML_dux_dxl(i,j) * time_n)
                  rmemory_acoustic_dux_dx(i,j,ispec_PML,2) = rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + &
                         BETA_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2)
                endif

                rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,1) = &
                       ALPHA_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,1) + &
                       deltat * (-bb_xz_1 * rmemory_acoustic_dux_dz(i,j,ispec_PML,1) + PML_dux_dzl(i,j))
                rmemory_acoustic_dux_dz(i,j,ispec_PML,1) = rmemory_acoustic_dux_dz(i,j,ispec_PML,1) + &
                       BETA_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,1)
                if (singularity_type_xz == 0) then
                  rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2) = &
                         ALPHA_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2) + &
                         deltat * (-bb_xz_2 * rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + PML_dux_dzl(i,j))
                  rmemory_acoustic_dux_dz(i,j,ispec_PML,2) = rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + &
                         BETA_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2)
                else
                  rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2) = &
                         ALPHA_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2) + &
                         deltat * (-bb_xz_2 * rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + PML_dux_dzl(i,j) * time_n)
                  rmemory_acoustic_dux_dz(i,j,ispec_PML,2) = rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + &
                         BETA_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2)
                endif

              endif

              dux_dxl = A5 * PML_dux_dxl(i,j) + A6 * rmemory_acoustic_dux_dx(i,j,ispec_PML,1) + &
                                                A7 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2)
              dux_dzl = A8 * PML_dux_dzl(i,j) + A9 *  rmemory_acoustic_dux_dz(i,j,ispec_PML,1) + &
                                                A10 * rmemory_acoustic_dux_dz(i,j,ispec_PML,2)
            endif
          endif ! PML_BOUNDARY_CONDITIONS

          ! AXISYM case
          if (AXISYM) then
            if (is_on_the_axis(ispec) .and. i == 1) then ! dchi/dr=rho * u_r=0 on the axis
              dux_dxl = ZERO
            endif
          endif

          jacobianl = jacobian(i,j,ispec)

          ! if external density model
          if (assign_external_model) then
            rhol = rhoext(i,j,ispec)
          endif

          if (AXISYM) then
            if (is_on_the_axis(ispec) .and. i == 1) then
              xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
              r_xiplus1(i,j) = xxi
            else if (is_on_the_axis(ispec)) then
              r_xiplus1(i,j) = coord(1,ibool(i,j,ispec))/(xiglj(i)+ONE)
            endif
          endif

          ! for acoustic medium also add integration weights
          if (AXISYM) then
            if (is_on_the_axis(ispec)) then
              tempx1(i,j) = wzgll(j) * r_xiplus1(i,j) * jacobianl * (xixl * dux_dxl + xizl * dux_dzl) / rhol
              tempx2(i,j) = wxglj(i) * r_xiplus1(i,j) * jacobianl * (gammaxl * dux_dxl + gammazl * dux_dzl) / rhol
            else
              tempx1(i,j) = wzgll(j) * coord(1,ibool(i,j,ispec)) * jacobianl * (xixl * dux_dxl + xizl * dux_dzl) / rhol
              tempx2(i,j) = wxgll(i) * coord(1,ibool(i,j,ispec)) * jacobianl * (gammaxl * dux_dxl + gammazl * dux_dzl) / rhol
            endif
          else
            tempx1(i,j) = wzgll(j) * jacobianl * (xixl * dux_dxl + xizl * dux_dzl) / rhol
            tempx2(i,j) = wxgll(i) * jacobianl * (gammaxl * dux_dxl + gammazl * dux_dzl) / rhol
          endif
        enddo
      enddo

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          if (assign_external_model) then
            rhol = rhoext(i,j,ispec)
            cpl = vpext(i,j,ispec)
            !assuming that in fluid(acoustic) part input cpl is defined by sqrt(kappal/rhol), &
            !which is not the same as in cpl input in elastic part
            kappal = rhol * cpl * cpl ! CHECK Kappa
          else
            lambdal_relaxed = poroelastcoef(1,1,kmato(ispec))
            mul_relaxed = poroelastcoef(2,1,kmato(ispec))
            if (AXISYM) then ! CHECK kappa
              kappal  = lambdal_relaxed + TWO * mul_relaxed/3._CUSTOM_REAL
            else
              kappal  = lambdal_relaxed + mul_relaxed
            endif
            rhol = density(1,kmato(ispec))
          endif

          if (PML_BOUNDARY_CONDITIONS) then
            if (ispec_is_PML(ispec) .and. nspec_PML > 0) then
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
              ! the subroutine of l_parameter_computation is located at the end of compute_forces_viscoelastic.F90
              call l_parameter_computation(time_n,deltat,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                           CPML_region_local,A0,A1,A2,A3,A4,singularity_type,&
                                           bb_1,coef0_1,coef1_1,coef2_1,bb_2,coef0_2,coef1_2,coef2_2)

              if (stage_time_scheme == 1) then
                rmemory_potential_acoustic(1,i,j,ispec_PML) = coef0_1 * rmemory_potential_acoustic(1,i,j,ispec_PML) + &
                       coef1_1 * potential_acoustic(iglob) + coef2_1 * potential_acoustic_old(iglob)

                if (singularity_type == 0) then
                  rmemory_potential_acoustic(2,i,j,ispec_PML) = coef0_2 * rmemory_potential_acoustic(2,i,j,ispec_PML) + &
                         coef1_2 * potential_acoustic(iglob) + coef2_2 * potential_acoustic_old(iglob)
                else
                  rmemory_potential_acoustic(2,i,j,ispec_PML) = coef0_2 * rmemory_potential_acoustic(2,i,j,ispec_PML) + &
                         coef1_2 * time_n * potential_acoustic(iglob) + coef2_2 * time_nsub1 * potential_acoustic_old(iglob)
                endif
              endif

              if (stage_time_scheme == 6) then
                rmemory_potential_acoustic_LDDRK(1,i,j,ispec_PML) = &
                       ALPHA_LDDRK(i_stage) * rmemory_potential_acoustic_LDDRK(1,i,j,ispec_PML) + &
                       deltat * (-bb_1 * rmemory_potential_acoustic(1,i,j,ispec_PML) + potential_acoustic(iglob))
                rmemory_potential_acoustic(1,i,j,ispec_PML) = rmemory_potential_acoustic(1,i,j,ispec_PML) + &
                         BETA_LDDRK(i_stage) * rmemory_potential_acoustic_LDDRK(1,i,j,ispec_PML)
                if (singularity_type == 0) then
                  rmemory_potential_acoustic_LDDRK(2,i,j,ispec_PML) = &
                         ALPHA_LDDRK(i_stage) * rmemory_potential_acoustic_LDDRK(2,i,j,ispec_PML) + &
                         deltat * (-bb_2 * rmemory_potential_acoustic(2,i,j,ispec_PML) + potential_acoustic(iglob))
                  rmemory_potential_acoustic(2,i,j,ispec_PML) = rmemory_potential_acoustic(2,i,j,ispec_PML) + &
                         BETA_LDDRK(i_stage) * rmemory_potential_acoustic_LDDRK(2,i,j,ispec_PML)
                else
                  rmemory_potential_acoustic_LDDRK(2,i,j,ispec_PML) = &
                         ALPHA_LDDRK(i_stage) * rmemory_potential_acoustic_LDDRK(2,i,j,ispec_PML) + &
                         deltat * (-bb_2 * rmemory_potential_acoustic(2,i,j,ispec_PML) + potential_acoustic(iglob) * time_n)
                  rmemory_potential_acoustic(2,i,j,ispec_PML) = rmemory_potential_acoustic(2,i,j,ispec_PML) + &
                         BETA_LDDRK(i_stage) * rmemory_potential_acoustic_LDDRK(2,i,j,ispec_PML)
                endif
              endif

              if (AXISYM) then
                if (is_on_the_axis(ispec)) then
                  potential_dot_dot_acoustic_PML(i,j)= wxglj(i) * wzgll(j)/kappal * jacobian(i,j,ispec) * r_xiplus1(i,j) * &
                           (A1 * potential_dot_acoustic(iglob) + A2 * potential_acoustic(iglob) + &
                            A3 * rmemory_potential_acoustic(1,i,j,ispec_PML) + &
                            A4 * rmemory_potential_acoustic(2,i,j,ispec_PML))
                else
                  potential_dot_dot_acoustic_PML(i,j)= wxgll(i) * wzgll(j)/kappal * jacobian(i,j,ispec) &
                                                      * coord(1,ibool(i,j,ispec)) * &
                           (A1 * potential_dot_acoustic(iglob) + A2 * potential_acoustic(iglob) + &
                            A3 * rmemory_potential_acoustic(1,i,j,ispec_PML) + &
                            A4 * rmemory_potential_acoustic(2,i,j,ispec_PML))
                endif
              else
                  potential_dot_dot_acoustic_PML(i,j)= wxgll(i) * wzgll(j)/kappal * jacobian(i,j,ispec) * &
                           (A1 * potential_dot_acoustic(iglob) + A2 * potential_acoustic(iglob) + &
                            A3 * rmemory_potential_acoustic(1,i,j,ispec_PML) + A4 * rmemory_potential_acoustic(2,i,j,ispec_PML))
              endif
            endif
          endif ! PML_BOUNDARY_CONDITIONS

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
          if (AXISYM) then
            if (is_on_the_axis(ispec)) then
              do k = 1,NGLLX
                potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
                       (tempx1(k,j) * hprimeBarwglj_xx(k,i) + tempx2(i,k) * hprimewgll_zz(k,j))
              enddo
            else
              do k = 1,NGLLX
                potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
                         (tempx1(k,j) * hprimewgll_xx(k,i) + tempx2(i,k) * hprimewgll_zz(k,j))
              enddo
            endif
          else
            do k = 1,NGLLX
              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
                       (tempx1(k,j) * hprimewgll_xx(k,i) + tempx2(i,k) * hprimewgll_zz(k,j))
            enddo
          endif

          if (PML_BOUNDARY_CONDITIONS) then
            if (ispec_is_PML(ispec)) then
              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_dot_acoustic_PML(i,j)
            endif
          endif
        enddo ! second loop over the GLL points
      enddo

    endif ! end of test if acoustic element
  enddo ! end of loop over all spectral elements

  end subroutine compute_forces_acoustic
