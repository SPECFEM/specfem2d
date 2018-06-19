!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

! calculates contribution from each C-PML element to update acceleration to the global mesh

! second-order accurate convolution term calculation from equation (21) of
! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
! Anisotropic-medium PML for vector FETD with modified basis functions,
! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)


  subroutine pml_compute_accel_contribution_acoustic(ispec,nglob, &
                                                     potential_acoustic,potential_acoustic_old,potential_dot_acoustic, &
                                                     potential_dot_dot_acoustic_PML,r_xiplus1)

! for acoustic elements

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,CPML_X_ONLY,CPML_Z_ONLY,ALPHA_LDDRK,BETA_LDDRK,C_LDDRK, &
    NGLJ,TWO_THIRDS

  use specfem_par, only: time_stepping_scheme,i_stage,it,deltat, &
                         assign_external_model,rhoext,vpext,density,poroelastcoef,kmato, &
                         ibool,jacobian, &
                         wxgll,wzgll, &
                         rmemory_potential_acoustic,rmemory_potential_acoustic_LDDRK, &
                         AXISYM,is_on_the_axis,coord,wxglj

  ! PML arrays
  use specfem_par, only: nspec_PML,ispec_is_PML,spec_to_PML,region_CPML, &
                K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store

  implicit none

  integer,intent(in) :: ispec

  integer,intent(in) :: nglob
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: potential_acoustic,potential_acoustic_old,potential_dot_acoustic

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(out) :: potential_dot_dot_acoustic_PML
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ),intent(in) :: r_xiplus1

  ! local parameters
  integer :: i,j,iglob

  ! PML
  integer :: ispec_PML
  integer :: CPML_region_local
  double precision :: kappa_x,kappa_z,d_x,d_z,alpha_x,alpha_z,beta_x,beta_z
  double precision :: time_n,time_nsub1
  double precision :: A0,A1,A2,A3,A4,bb_1,coef0_1,coef1_1,coef2_1,bb_2,coef0_2,coef1_2,coef2_2

  ! material properties of the acoustic medium
  real(kind=CUSTOM_REAL) :: mul_relaxed,lambdal_relaxed,kappal,cpl,rhol
  real(kind=CUSTOM_REAL) :: fac

  ! checks if anything to do in this slice
  if (nspec_PML == 0) return
  if (.not. ispec_is_PML(ispec)) return

  ! timing
  select case (time_stepping_scheme)
  case (1)
    ! Newmark
    time_n = (it-1) * deltat
    time_nsub1 = (it-2) * deltat
  case (2)
    ! LDDRK
    time_n = (it-1) * deltat + C_LDDRK(i_stage) * deltat
  case default
    call stop_the_code('Sorry, time stepping scheme for PML accel contribution not implemented yet')
  end select

  ! local PML element
  ispec_PML = spec_to_PML(ispec)
  CPML_region_local = region_CPML(ispec)

  do j = 1,NGLLZ
    do i = 1,NGLLX

      kappa_x = K_x_store(i,j,ispec_PML)
      kappa_z = K_z_store(i,j,ispec_PML)

      d_x = d_x_store(i,j,ispec_PML)
      d_z = d_z_store(i,j,ispec_PML)

      alpha_x = alpha_x_store(i,j,ispec_PML)
      alpha_z = alpha_z_store(i,j,ispec_PML)

      beta_x = alpha_x + d_x / kappa_x
      beta_z = alpha_z + d_z / kappa_z

      ! the subroutine of l_parameter_computation is located at the end of compute_forces_viscoelastic.F90
      call l_parameter_computation(deltat,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                   CPML_region_local,A0,A1,A2,A3,A4, &
                                   bb_1,coef0_1,coef1_1,coef2_1,bb_2,coef0_2,coef1_2,coef2_2)

      iglob = ibool(i,j,ispec)

      select case (time_stepping_scheme)
      case (1)
        ! Newmark
        rmemory_potential_acoustic(1,i,j,ispec_PML) = coef0_1 * rmemory_potential_acoustic(1,i,j,ispec_PML) + &
               coef1_1 * potential_acoustic(iglob) + coef2_1 * potential_acoustic_old(iglob)

        rmemory_potential_acoustic(2,i,j,ispec_PML) = coef0_2 * rmemory_potential_acoustic(2,i,j,ispec_PML) + &
               coef1_2 * potential_acoustic(iglob) + coef2_2 * potential_acoustic_old(iglob)

      case (2)
        ! LDDRK
        rmemory_potential_acoustic_LDDRK(1,i,j,ispec_PML) = &
               ALPHA_LDDRK(i_stage) * rmemory_potential_acoustic_LDDRK(1,i,j,ispec_PML) + &
               deltat * (-bb_1 * rmemory_potential_acoustic(1,i,j,ispec_PML) + potential_acoustic(iglob))
        rmemory_potential_acoustic(1,i,j,ispec_PML) = rmemory_potential_acoustic(1,i,j,ispec_PML) + &
                 BETA_LDDRK(i_stage) * rmemory_potential_acoustic_LDDRK(1,i,j,ispec_PML)

        rmemory_potential_acoustic_LDDRK(2,i,j,ispec_PML) = &
               ALPHA_LDDRK(i_stage) * rmemory_potential_acoustic_LDDRK(2,i,j,ispec_PML) + &
               deltat * (-bb_2 * rmemory_potential_acoustic(2,i,j,ispec_PML) + potential_acoustic(iglob))
        rmemory_potential_acoustic(2,i,j,ispec_PML) = rmemory_potential_acoustic(2,i,j,ispec_PML) + &
               BETA_LDDRK(i_stage) * rmemory_potential_acoustic_LDDRK(2,i,j,ispec_PML)
      end select

      ! material properties
      if (assign_external_model) then
        rhol = rhoext(i,j,ispec)
        cpl = vpext(i,j,ispec)
        !assuming that in fluid(acoustic) part input cpl is defined by sqrt(kappal/rhol), &
        !which is not the same as in cpl input in elastic part
        kappal = rhol * cpl * cpl ! CHECK Kappa : it is ok here because we are in acoustic elements
      else
        rhol = density(1,kmato(ispec))
        lambdal_relaxed = poroelastcoef(1,1,kmato(ispec))
        mul_relaxed = poroelastcoef(2,1,kmato(ispec))
        if (AXISYM) then ! CHECK kappa
          kappal  = lambdal_relaxed + TWO_THIRDS * mul_relaxed
        else
          kappal  = lambdal_relaxed + mul_relaxed
        endif
      endif

      fac = jacobian(i,j,ispec) / kappal

      if (AXISYM) then
        if (is_on_the_axis(ispec)) then
          potential_dot_dot_acoustic_PML(i,j)= wxglj(i) * wzgll(j) * fac * r_xiplus1(i,j) * &
                   (A1 * potential_dot_acoustic(iglob) + A2 * potential_acoustic(iglob) + &
                    A3 * rmemory_potential_acoustic(1,i,j,ispec_PML) + &
                    A4 * rmemory_potential_acoustic(2,i,j,ispec_PML))
        else
          potential_dot_dot_acoustic_PML(i,j)= wxgll(i) * wzgll(j) * fac * coord(1,ibool(i,j,ispec)) * &
                   (A1 * potential_dot_acoustic(iglob) + A2 * potential_acoustic(iglob) + &
                    A3 * rmemory_potential_acoustic(1,i,j,ispec_PML) + &
                    A4 * rmemory_potential_acoustic(2,i,j,ispec_PML))
        endif
      else
        potential_dot_dot_acoustic_PML(i,j)= wxgll(i) * wzgll(j) * fac * &
                  (A1 * potential_dot_acoustic(iglob) + A2 * potential_acoustic(iglob) + &
                  A3 * rmemory_potential_acoustic(1,i,j,ispec_PML) + A4 * rmemory_potential_acoustic(2,i,j,ispec_PML))
      endif
    enddo
  enddo

  end subroutine pml_compute_accel_contribution_acoustic

!
!-------------------------------------------------------------------------------------
!


  subroutine pml_compute_accel_contribution_elastic(ispec,nglob, &
                                                    dummy_loc,displ_elastic_old,veloc_elastic, &
                                                    accel_elastic_PML,r_xiplus1)
! for elastic elements

  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLZ,NGLJ,TWO_THIRDS, &
                       CPML_X_ONLY,CPML_Z_ONLY,ALPHA_LDDRK,BETA_LDDRK,C_LDDRK

  use specfem_par, only: time_stepping_scheme,i_stage,deltat, &
                         assign_external_model,rhoext,density,kmato, &
                         ibool,jacobian,wxgll,wzgll, &
                         AXISYM,is_on_the_axis,coord,wxglj, &
                         rmemory_displ_elastic,rmemory_displ_elastic_LDDRK

  ! PML arrays
  use specfem_par, only: nspec_PML,ispec_is_PML,spec_to_PML,region_CPML, &
                K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store

  implicit none

  integer,intent(in) :: ispec,nglob

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ),intent(in) :: dummy_loc
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: veloc_elastic,displ_elastic_old

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ),intent(out) :: accel_elastic_PML

  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ),intent(in) :: r_xiplus1

  ! local parameters
  integer :: i,j,iglob

  ! PML
  integer :: ispec_PML,CPML_region_local
  double precision :: kappa_x,kappa_z,d_x,d_z,alpha_x,alpha_z,beta_x,beta_z
  double precision :: A0,A1,A2,A3,A4,bb_1,coef0_1,coef1_1,coef2_1,bb_2,coef0_2,coef1_2,coef2_2

  ! material properties of the acoustic medium
  real(kind=CUSTOM_REAL) :: rhol

  ! checks if anything to do in this slice
  if (nspec_PML == 0) return
  if (.not. ispec_is_PML(ispec)) return

  ! initializes
  accel_elastic_PML(:,:,:) = 0._CUSTOM_REAL

  ! local PML element
  ispec_PML = spec_to_PML(ispec)
  CPML_region_local = region_CPML(ispec)

  do j = 1,NGLLZ
    do i = 1,NGLLX

      iglob = ibool(i,j,ispec)

      kappa_x = K_x_store(i,j,ispec_PML)
      kappa_z = K_z_store(i,j,ispec_PML)

      d_x = d_x_store(i,j,ispec_PML)
      d_z = d_z_store(i,j,ispec_PML)

      alpha_x = alpha_x_store(i,j,ispec_PML)
      alpha_z = alpha_z_store(i,j,ispec_PML)

      beta_x = alpha_x + d_x / kappa_x
      beta_z = alpha_z + d_z / kappa_z

      call l_parameter_computation(deltat,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                   CPML_region_local,A0,A1,A2,A3,A4, &
                                   bb_1,coef0_1,coef1_1,coef2_1,bb_2,coef0_2,coef1_2,coef2_2)

      select case (time_stepping_scheme)
      case (1)
        ! Newmark
        rmemory_displ_elastic(1,1,i,j,ispec_PML) = coef0_1 * rmemory_displ_elastic(1,1,i,j,ispec_PML) + &
                                                   coef1_1 * dummy_loc(1,i,j) + coef2_1 * displ_elastic_old(1,iglob)
        rmemory_displ_elastic(1,2,i,j,ispec_PML) = coef0_1 * rmemory_displ_elastic(1,2,i,j,ispec_PML) + &
                                                   coef1_1 * dummy_loc(2,i,j) + coef2_1 * displ_elastic_old(2,iglob)

        rmemory_displ_elastic(2,1,i,j,ispec_PML) = coef0_2 * rmemory_displ_elastic(2,1,i,j,ispec_PML) + &
                                                   coef1_2 * dummy_loc(1,i,j) + coef2_2 * displ_elastic_old(1,iglob)
        rmemory_displ_elastic(2,2,i,j,ispec_PML) = coef0_2 * rmemory_displ_elastic(2,2,i,j,ispec_PML) + &
                                                   coef1_2 * dummy_loc(2,i,j) + coef2_2 * displ_elastic_old(2,iglob)

      case (2)
        ! LDDRK
        rmemory_displ_elastic_LDDRK(1,1,i,j,ispec_PML) = &
              ALPHA_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(1,1,i,j,ispec_PML) + &
              deltat * (-bb_1 * rmemory_displ_elastic(1,1,i,j,ispec_PML) + dummy_loc(1,i,j))
        rmemory_displ_elastic(1,1,i,j,ispec_PML) = rmemory_displ_elastic(1,1,i,j,ispec_PML) + &
              BETA_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(1,1,i,j,ispec_PML)

        rmemory_displ_elastic_LDDRK(1,2,i,j,ispec_PML) = &
              ALPHA_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(1,2,i,j,ispec_PML) + &
              deltat * (-bb_1 * rmemory_displ_elastic(1,2,i,j,ispec_PML) + dummy_loc(2,i,j))
        rmemory_displ_elastic(1,2,i,j,ispec_PML) = rmemory_displ_elastic(1,2,i,j,ispec_PML) + &
              BETA_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(1,2,i,j,ispec_PML)

        rmemory_displ_elastic_LDDRK(2,1,i,j,ispec_PML) = &
              ALPHA_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,1,i,j,ispec_PML) + &
              deltat * (-bb_2 * rmemory_displ_elastic(2,1,i,j,ispec_PML) + dummy_loc(1,i,j))
        rmemory_displ_elastic(2,1,i,j,ispec_PML) = rmemory_displ_elastic(2,1,i,j,ispec_PML) + &
              BETA_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,1,i,j,ispec_PML)

        rmemory_displ_elastic_LDDRK(2,2,i,j,ispec_PML) = &
              ALPHA_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,2,i,j,ispec_PML) + &
              deltat * (-bb_2 * rmemory_displ_elastic(2,2,i,j,ispec_PML) + dummy_loc(2,i,j))
        rmemory_displ_elastic(2,2,i,j,ispec_PML) = rmemory_displ_elastic(2,2,i,j,ispec_PML) + &
              BETA_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,2,i,j,ispec_PML)
      case default
        call stop_the_code('Time stepping scheme not implemented yet for PML accel contribution')
      end select

      if (assign_external_model) then
        rhol = rhoext(i,j,ispec)
      else
        rhol = density(1,kmato(ispec))
      endif

      if (AXISYM) then
        if (is_on_the_axis(ispec)) then
          accel_elastic_PML(1,i,j) = wxglj(i)*wzgll(j)*rhol*jacobian(i,j,ispec)*r_xiplus1(i,j) * &
               ( A1 * veloc_elastic(1,iglob) + A2 * dummy_loc(1,i,j) + &
                 A3 * rmemory_displ_elastic(1,1,i,j,ispec_PML) + A4 * rmemory_displ_elastic(2,1,i,j,ispec_PML))
          accel_elastic_PML(2,i,j) = wxglj(i)*wzgll(j)*rhol*jacobian(i,j,ispec)*r_xiplus1(i,j) * &
               ( A1 * veloc_elastic(2,iglob) + A2 * dummy_loc(2,i,j) + &
                 A3 * rmemory_displ_elastic(1,2,i,j,ispec_PML) + A4 * rmemory_displ_elastic(2,2,i,j,ispec_PML))
        else
          accel_elastic_PML(1,i,j) = wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec)*coord(1,ibool(i,j,ispec))* &
               ( A1 * veloc_elastic(1,iglob) + A2 * dummy_loc(1,i,j) + &
                 A3 * rmemory_displ_elastic(1,1,i,j,ispec_PML) + A4 * rmemory_displ_elastic(2,1,i,j,ispec_PML))
          accel_elastic_PML(2,i,j) = wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec)*coord(1,ibool(i,j,ispec))* &
               ( A1 * veloc_elastic(2,iglob) + A2 * dummy_loc(2,i,j) + &
                 A3 * rmemory_displ_elastic(1,2,i,j,ispec_PML) + A4 * rmemory_displ_elastic(2,2,i,j,ispec_PML))
        endif
      else
        ! default
        accel_elastic_PML(1,i,j) = wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * &
              ( A1 * veloc_elastic(1,iglob) + A2 * dummy_loc(1,i,j) + &
                A3 * rmemory_displ_elastic(1,1,i,j,ispec_PML) + A4 * rmemory_displ_elastic(2,1,i,j,ispec_PML))
        accel_elastic_PML(2,i,j) = wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * &
              ( A1 * veloc_elastic(2,iglob) + A2 * dummy_loc(2,i,j) + &
                A3 * rmemory_displ_elastic(1,2,i,j,ispec_PML) + A4 * rmemory_displ_elastic(2,2,i,j,ispec_PML))
      endif
    enddo
  enddo

  end subroutine pml_compute_accel_contribution_elastic

