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

  subroutine prepare_wavefields()

  use constants, only: IMAIN,USE_ENFORCE_FIELDS
  use specfem_par

  implicit none

  ! local parameters
  integer :: b_nglob_acoustic,b_nglob_elastic,b_nglob_poroelastic
  integer :: b_nspec_acoustic,b_nspec_elastic
  integer :: ier

  ! displacement, velocity, acceleration and inverse of the mass matrix for elastic elements
  if (myrank == 0) then
    write(IMAIN,*) 'Preparing array allocations'
    call flush_IMAIN()
  endif

! note: allocating fields does not yet map the arrays onto memory.
!       only when we first assign values to the arrays, the actually memory addresses get determined.
!
!       therefore, changing the order here could potentially change how near the next "neighbor" field address gets located,
!       which in turn could have an effect on data access/movement and thus code performance.

  !
  ! acoustic domains
  !
  if (ACOUSTIC_SIMULATION) then
    if (myrank == 0) then
      write(IMAIN,*) '  arrays for acoustic domains'
      call flush_IMAIN()
    endif
  endif

  ! potential, its first and second derivative, and inverse of the mass matrix for acoustic elements
  if (any_acoustic) then
    nglob_acoustic = nglob
  else
    ! dummy allocate unused arrays with fictitious size
    nglob_acoustic = 1
  endif

  allocate(potential_acoustic(nglob_acoustic), &
           potential_dot_acoustic(nglob_acoustic), &
           potential_dot_dot_acoustic(nglob_acoustic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating acoustic wavefield arrays')
  potential_acoustic(:) = 0._CUSTOM_REAL; potential_dot_acoustic(:) = 0._CUSTOM_REAL;
  potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL

  ! PML
  allocate(potential_acoustic_old(nglob_acoustic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating potential_acoustic_old array')
  potential_acoustic_old(:) = 0._CUSTOM_REAL

  allocate(rmass_inverse_acoustic(nglob_acoustic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating rmass_inverse_acoustic array')
  rmass_inverse_acoustic(:) = 0._CUSTOM_REAL

  ! intermediate fields
  ! LDDRK
  if (time_stepping_scheme == 2) then
    allocate(potential_acoustic_LDDRK(nglob_acoustic), &
             potential_dot_acoustic_LDDRK(nglob_acoustic), &
             potential_dot_acoustic_temp(nglob_acoustic),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating acoustic LDDRK arrays')
    potential_acoustic_LDDRK(:) = 0._CUSTOM_REAL; potential_dot_acoustic_LDDRK(:) = 0._CUSTOM_REAL
    potential_dot_acoustic_temp(:) = 0._CUSTOM_REAL
  endif

  ! RK4
  if (time_stepping_scheme == 3) then
    allocate(potential_acoustic_init_rk(nglob_acoustic), &
             potential_dot_acoustic_init_rk(nglob_acoustic), &
             potential_dot_dot_acoustic_rk(nglob_acoustic,NSTAGE_TIME_SCHEME), &
             potential_dot_acoustic_rk(nglob_acoustic,NSTAGE_TIME_SCHEME),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating acoustic RK arrays')
    potential_acoustic_init_rk(:) = 0._CUSTOM_REAL; potential_dot_acoustic_init_rk(:) = 0._CUSTOM_REAL
    potential_dot_dot_acoustic_rk(:,:) = 0._CUSTOM_REAL; potential_dot_acoustic_rk(:,:) = 0._CUSTOM_REAL
  endif

!! DK DK March 2018: this was missing in Quentin Brissaud's new variational formulation for viscoacoustic media; I added it
  if (ATTENUATION_VISCOACOUSTIC) then
    allocate(rmass_inverse_e1(nglob_acoustic,N_SLS),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating rmass_inverse_e1 arrays')
  else
    allocate(rmass_inverse_e1(1,1))
  endif
  rmass_inverse_e1(:,:) = 0.0_CUSTOM_REAL

  if (SIMULATION_TYPE == 3 .and. any_acoustic) then
    b_nglob_acoustic = nglob
    b_nspec_acoustic = nspec
  else
    ! dummy array allocations
    ! allocates unused arrays with fictitious size
    b_nglob_acoustic = 1
    b_nspec_acoustic = 1
  endif
  allocate(b_potential_acoustic(b_nglob_acoustic), &
           b_potential_dot_acoustic(b_nglob_acoustic), &
           b_potential_dot_dot_acoustic(b_nglob_acoustic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating acoustic backward wavefield arrays')
  b_potential_acoustic(:) = 0._CUSTOM_REAL; b_potential_dot_acoustic(:) = 0._CUSTOM_REAL
  b_potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL

  ! PML
  allocate(b_potential_acoustic_old(b_nglob_acoustic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating b_potential_acoustic_old arrays')
  b_potential_acoustic_old(:) = 0._CUSTOM_REAL

  ! coupling for SIMULATION_TYPE == 3
  ! not needed anymore, taking care of by re-ordering domain updates
  !allocate(potential_acoustic_adj_coupling(b_nglob_acoustic),stat=ier)
  !if (ier /= 0) call stop_the_code('Error allocating potential_acoustic_adj_coupling array')
  !potential_acoustic_adj_coupling(:) = 0.0_CUSTOM_REAL

  ! acoustic kernels
  allocate(rho_ac_kl(NGLLX,NGLLZ,b_nspec_acoustic), &
           kappa_ac_kl(NGLLX,NGLLZ,b_nspec_acoustic), &
           rhop_ac_kl(NGLLX,NGLLZ,b_nspec_acoustic), &
           alpha_ac_kl(NGLLX,NGLLZ,b_nspec_acoustic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating rho_ac_kl arrays')
  rho_ac_kl(:,:,:) = 0.0_CUSTOM_REAL; kappa_ac_kl(:,:,:) = 0.0_CUSTOM_REAL
  rhop_ac_kl(:,:,:) = 0.0_CUSTOM_REAL; alpha_ac_kl(:,:,:) = 0.0_CUSTOM_REAL

  ! for APPROXIMATE_HESS_KL
  allocate(rhorho_ac_Hessian_final2(NGLLX,NGLLZ,b_nspec_acoustic), &
           rhorho_ac_Hessian_final1(NGLLX,NGLLZ,b_nspec_acoustic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating rhorho_ac_Hessian_final arrays')
  rhorho_ac_Hessian_final1(:,:,:) = 0.0_CUSTOM_REAL; rhorho_ac_Hessian_final2(:,:,:) = 0.0_CUSTOM_REAL

  !
  ! elastic domains
  !
  ! user info
  if (ELASTIC_SIMULATION) then
    if (myrank == 0) then
      write(IMAIN,*) '  arrays for elastic domains'
      call flush_IMAIN()
    endif
  endif

  ! sets global points in this slice
  if (any_elastic) then
    nglob_elastic = nglob
  else
    ! dummy allocate unused arrays with fictitious size
    nglob_elastic = 1
  endif

  allocate(displ_elastic(NDIM,nglob_elastic), &
           veloc_elastic(NDIM,nglob_elastic), &
           accel_elastic(NDIM,nglob_elastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating elastic wavefield arrays')
  displ_elastic(:,:) = 0.0_CUSTOM_REAL; veloc_elastic(:,:) = 0.0_CUSTOM_REAL; accel_elastic(:,:) = 0.0_CUSTOM_REAL

  ! PML
  allocate(displ_elastic_old(NDIM,nglob_elastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating old elastic wavefield arrays')
  displ_elastic_old(:,:) = 0.0_CUSTOM_REAL

  allocate(rmass_inverse_elastic(NDIM,nglob_elastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating elastic mass matrix array')
  rmass_inverse_elastic(:,:) = 0.0_CUSTOM_REAL

  ! intermediate wavefields
  ! LDDRK
  if (time_stepping_scheme == 2) then
    allocate(displ_elastic_LDDRK(NDIM,nglob_elastic), &
             veloc_elastic_LDDRK(NDIM,nglob_elastic), &
             veloc_elastic_LDDRK_temp(NDIM,nglob_elastic),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating elastic LDDRK wavefield arrays')
    displ_elastic_LDDRK(:,:) = 0.0_CUSTOM_REAL; veloc_elastic_LDDRK(:,:) = 0.0_CUSTOM_REAL
    veloc_elastic_LDDRK_temp(:,:) = 0.0_CUSTOM_REAL
  endif

  ! RK4
  if (time_stepping_scheme == 3) then
    allocate(accel_elastic_rk(NDIM,nglob_elastic,NSTAGE_TIME_SCHEME), &
             veloc_elastic_rk(NDIM,nglob_elastic,NSTAGE_TIME_SCHEME), &
             veloc_elastic_initial_rk(NDIM,nglob_elastic), &
             displ_elastic_initial_rk(NDIM,nglob_elastic),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating elastic RK wavefield arrays')
    accel_elastic_rk(:,:,:) = 0.0_CUSTOM_REAL; veloc_elastic_rk(:,:,:) = 0.0_CUSTOM_REAL
    veloc_elastic_initial_rk(:,:) = 0.0_CUSTOM_REAL; displ_elastic_initial_rk(:,:) = 0.0_CUSTOM_REAL
  endif

  ! extra array if adjoint and kernels calculation
  if (SIMULATION_TYPE == 3 .and. any_elastic) then
    b_nglob_elastic = nglob
    b_nspec_elastic = nspec
  else
    ! dummy allocations
    b_nglob_elastic = 1
    b_nspec_elastic = 1
  endif

  allocate(b_displ_elastic(NDIM,b_nglob_elastic), &
           b_veloc_elastic(NDIM,b_nglob_elastic), &
           b_accel_elastic(NDIM,b_nglob_elastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating elastic backward wavefield arrays')
  b_displ_elastic(:,:) = 0.0_CUSTOM_REAL; b_veloc_elastic(:,:) = 0.0_CUSTOM_REAL; b_accel_elastic(:,:) = 0.0_CUSTOM_REAL

  ! PML
  allocate(b_displ_elastic_old(NDIM,b_nglob_elastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating b_displ_elastic_old array')
  b_displ_elastic_old(:,:) = 0.0_CUSTOM_REAL

  ! for SIMULATION_TYPE == 3 and coupled_acoustic_elastic
  ! not needed anymore, taking care of by re-ordering domain updates
  !allocate(accel_elastic_adj_coupling(NDIM,b_nglob_elastic),stat=ier)
  !if (ier /= 0) call stop_the_code('Error allocating accel_elastic_adj_coupling array')
  !accel_elastic_adj_coupling(:,:) = 0.0_CUSTOM_REAL

  ! kernels
  allocate(rho_kl(NGLLX,NGLLZ,b_nspec_elastic), &
           mu_kl(NGLLX,NGLLZ,b_nspec_elastic), &
           kappa_kl(NGLLX,NGLLZ,b_nspec_elastic), &
           rhop_kl(NGLLX,NGLLZ,b_nspec_elastic), &
           alpha_kl(NGLLX,NGLLZ,b_nspec_elastic), &
           beta_kl(NGLLX,NGLLZ,b_nspec_elastic), &
           bulk_c_kl(NGLLX,NGLLZ,b_nspec_elastic), &
           bulk_beta_kl(NGLLX,NGLLZ,b_nspec_elastic), &
           c11_kl(NGLLX,NGLLZ,b_nspec_elastic), &
           c13_kl(NGLLX,NGLLZ,b_nspec_elastic), &
           c15_kl(NGLLX,NGLLZ,b_nspec_elastic), &
           c33_kl(NGLLX,NGLLZ,b_nspec_elastic), &
           c35_kl(NGLLX,NGLLZ,b_nspec_elastic), &
           c55_kl(NGLLX,NGLLZ,b_nspec_elastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating elastic nspec kernel arrays')
  rho_kl(:,:,:) = 0.0_CUSTOM_REAL; mu_kl(:,:,:) = 0.0_CUSTOM_REAL; kappa_kl(:,:,:) = 0.0_CUSTOM_REAL
  rhop_kl(:,:,:) = 0.0_CUSTOM_REAL; alpha_kl(:,:,:) = 0.0_CUSTOM_REAL; beta_kl(:,:,:) = 0.0_CUSTOM_REAL
  bulk_c_kl(:,:,:) = 0.0_CUSTOM_REAL; bulk_beta_kl(:,:,:) = 0.0_CUSTOM_REAL
  c11_kl(:,:,:) = 0.0_CUSTOM_REAL; c13_kl(:,:,:) = 0.0_CUSTOM_REAL; c15_kl(:,:,:) = 0.0_CUSTOM_REAL
  c33_kl(:,:,:) = 0.0_CUSTOM_REAL; c35_kl(:,:,:) = 0.0_CUSTOM_REAL; c55_kl(:,:,:) = 0.0_CUSTOM_REAL

  ! for APPROXIMATE_HESS_KL
  allocate(rhorho_el_Hessian_final2(NGLLX,NGLLZ,b_nspec_elastic), &
           rhorho_el_Hessian_final1(NGLLX,NGLLZ,b_nspec_elastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating elastic Hessian kernel arrays')
  rhorho_el_Hessian_final1(:,:,:) = 0.0_CUSTOM_REAL; rhorho_el_Hessian_final2(:,:,:) = 0.0_CUSTOM_REAL

  !
  ! poro-elastic domains
  !
  if (POROELASTIC_SIMULATION) then
    if (myrank == 0) then
      write(IMAIN,*) '  arrays for poroelastic domains'
      call flush_IMAIN()
    endif
  endif

  ! sets number of points for this slice
  if (any_poroelastic) then
    nglob_poroelastic = nglob
  else
    ! dummy allocate unused arrays with fictitious size
    nglob_poroelastic = 1
  endif
  allocate(displs_poroelastic(NDIM,nglob_poroelastic), &
           velocs_poroelastic(NDIM,nglob_poroelastic), &
           accels_poroelastic(NDIM,nglob_poroelastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating poroelastic wavefield arrays')
  displs_poroelastic(:,:) = 0._CUSTOM_REAL; velocs_poroelastic(:,:) = 0._CUSTOM_REAL; accels_poroelastic(:,:) = 0._CUSTOM_REAL

  ! PML
  allocate(displs_poroelastic_old(NDIM,nglob_poroelastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating displs_poroelastic_old array')
  displs_poroelastic_old(:,:) = 0._CUSTOM_REAL

  allocate(displw_poroelastic(NDIM,nglob_poroelastic), &
           velocw_poroelastic(NDIM,nglob_poroelastic), &
           accelw_poroelastic(NDIM,nglob_poroelastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating poroelastic fluid wavefield arrays')
  displw_poroelastic(:,:) = 0._CUSTOM_REAL; velocw_poroelastic(:,:) = 0._CUSTOM_REAL; accelw_poroelastic(:,:) = 0._CUSTOM_REAL

  allocate(rmass_s_inverse_poroelastic(nglob_poroelastic), &
           rmass_w_inverse_poroelastic(nglob_poroelastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating poroelastic mass arrays')
  rmass_s_inverse_poroelastic(:) = 0._CUSTOM_REAL
  rmass_w_inverse_poroelastic(:) = 0._CUSTOM_REAL

  ! intermediate fields
  ! LDDRK
  if (time_stepping_scheme == 2) then
    allocate(displs_poroelastic_LDDRK(NDIM,nglob_poroelastic), &
             velocs_poroelastic_LDDRK(NDIM,nglob_poroelastic), &
             displw_poroelastic_LDDRK(NDIM,nglob_poroelastic), &
             velocw_poroelastic_LDDRK(NDIM,nglob_poroelastic),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating poroelastic LDDRK arrays')
    displs_poroelastic_LDDRK(:,:) = 0._CUSTOM_REAL; velocs_poroelastic_LDDRK(:,:) = 0._CUSTOM_REAL
    displw_poroelastic_LDDRK(:,:) = 0._CUSTOM_REAL; velocw_poroelastic_LDDRK(:,:) = 0._CUSTOM_REAL
  endif

  ! RK4
  if (time_stepping_scheme == 3) then
    allocate(accels_poroelastic_rk(NDIM,nglob_poroelastic,NSTAGE_TIME_SCHEME), &
             velocs_poroelastic_rk(NDIM,nglob_poroelastic,NSTAGE_TIME_SCHEME), &
             accelw_poroelastic_rk(NDIM,nglob_poroelastic,NSTAGE_TIME_SCHEME), &
             velocw_poroelastic_rk(NDIM,nglob_poroelastic,NSTAGE_TIME_SCHEME), &
             displs_poroelastic_initial_rk(NDIM,nglob_poroelastic), &
             velocs_poroelastic_initial_rk(NDIM,nglob_poroelastic), &
             displw_poroelastic_initial_rk(NDIM,nglob_poroelastic), &
             velocw_poroelastic_initial_rk(NDIM,nglob_poroelastic),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating poroelastic RK arrays')
    accels_poroelastic_rk(:,:,:) = 0._CUSTOM_REAL; velocs_poroelastic_rk(:,:,:) = 0._CUSTOM_REAL
    accelw_poroelastic_rk(:,:,:) = 0._CUSTOM_REAL; velocw_poroelastic_rk(:,:,:) = 0._CUSTOM_REAL
    velocs_poroelastic_initial_rk(:,:) = 0._CUSTOM_REAL; displs_poroelastic_initial_rk(:,:) = 0._CUSTOM_REAL
    velocw_poroelastic_initial_rk(:,:) = 0._CUSTOM_REAL; displw_poroelastic_initial_rk(:,:) = 0._CUSTOM_REAL
  endif

  ! extra array if adjoint and kernels calculation
  if (SIMULATION_TYPE == 3 .and. any_poroelastic) then
    b_nglob_poroelastic = nglob
    b_nspec_poroelastic = nspec
  else
    ! dummy allocations
    b_nglob_poroelastic = 1
    b_nspec_poroelastic = 1
  endif
  allocate(b_displs_poroelastic(NDIM,b_nglob_poroelastic), &
           b_velocs_poroelastic(NDIM,b_nglob_poroelastic), &
           b_accels_poroelastic(NDIM,b_nglob_poroelastic), &
           b_displw_poroelastic(NDIM,b_nglob_poroelastic), &
           b_velocw_poroelastic(NDIM,b_nglob_poroelastic), &
           b_accelw_poroelastic(NDIM,b_nglob_poroelastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating poroelastic backward wavefield arrays')
  b_displs_poroelastic(:,:) = 0.0_CUSTOM_REAL; b_velocs_poroelastic(:,:) = 0.0_CUSTOM_REAL
  b_accels_poroelastic(:,:) = 0.0_CUSTOM_REAL
  b_displw_poroelastic(:,:) = 0.0_CUSTOM_REAL; b_velocw_poroelastic(:,:) = 0.0_CUSTOM_REAL
  b_accelw_poroelastic(:,:) = 0.0_CUSTOM_REAL

  ! coupling when SIMULATION_TYPE == 3
  ! not needed anymore, taking care of by re-ordering domain updates
  !allocate(accels_poroelastic_adj_coupling(NDIM,b_nglob_poroelastic), &
  !         accelw_poroelastic_adj_coupling(NDIM,b_nglob_poroelastic),stat=ier)
  !if (ier /= 0) call stop_the_code('Error allocating accels_poroelastic_adj_coupling arrays')
  !accels_poroelastic_adj_coupling(:,:) = 0.0_CUSTOM_REAL; accelw_poroelastic_adj_coupling(:,:) = 0.0_CUSTOM_REAL

  ! strain for kernel simulations
  allocate(epsilondev_s(4,NGLLX,NGLLZ,b_nspec_poroelastic), &
           epsilondev_w(4,NGLLX,NGLLZ,b_nspec_poroelastic), &
           b_epsilondev_s(4,NGLLX,NGLLZ,b_nspec_poroelastic), &
           b_epsilondev_w(4,NGLLX,NGLLZ,b_nspec_poroelastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating epsilondev_s arrays')
  epsilondev_s(:,:,:,:) = 0.0_CUSTOM_REAL; epsilondev_w(:,:,:,:) = 0.0_CUSTOM_REAL
  b_epsilondev_s(:,:,:,:) = 0.0_CUSTOM_REAL; b_epsilondev_w(:,:,:,:) = 0.0_CUSTOM_REAL

  ! poroelastic kernels
  allocate(rhot_kl(NGLLX,NGLLZ,b_nspec_poroelastic), &
           rhof_kl(NGLLX,NGLLZ,b_nspec_poroelastic), &
           sm_kl(NGLLX,NGLLZ,b_nspec_poroelastic), &
           eta_kl(NGLLX,NGLLZ,b_nspec_poroelastic), &
           mufr_kl(NGLLX,NGLLZ,b_nspec_poroelastic), &
           B_kl(NGLLX,NGLLZ,b_nspec_poroelastic), &
           C_kl(NGLLX,NGLLZ,b_nspec_poroelastic), &
           M_kl(NGLLX,NGLLZ,b_nspec_poroelastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating rhot_kl kernel arrays')
  rhot_kl(:,:,:) = 0.0_CUSTOM_REAL; rhof_kl(:,:,:) = 0.0_CUSTOM_REAL
  sm_kl(:,:,:) = 0.0_CUSTOM_REAL; eta_kl(:,:,:) = 0.0_CUSTOM_REAL
  mufr_kl(:,:,:) = 0.0_CUSTOM_REAL; B_kl(:,:,:) = 0.0_CUSTOM_REAL
  C_kl(:,:,:) = 0.0_CUSTOM_REAL; M_kl(:,:,:) = 0.0_CUSTOM_REAL

  allocate(phi_kl(NGLLX,NGLLZ,b_nspec_poroelastic), &
           phib_kl(NGLLX,NGLLZ,b_nspec_poroelastic), &
           mufrb_kl(NGLLX,NGLLZ,b_nspec_poroelastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating phi_kl kernel arrays')
  phi_kl(:,:,:) = 0.0_CUSTOM_REAL; phib_kl(:,:,:) = 0.0_CUSTOM_REAL; mufrb_kl(:,:,:) = 0.0_CUSTOM_REAL

  allocate(rhob_kl(NGLLX,NGLLZ,b_nspec_poroelastic), &
           rhobb_kl(NGLLX,NGLLZ,b_nspec_poroelastic), &
           rhofb_kl(NGLLX,NGLLZ,b_nspec_poroelastic), &
           rhofbb_kl(NGLLX,NGLLZ,b_nspec_poroelastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating rhob_kl kernel arrays')
  rhob_kl(:,:,:) = 0.0_CUSTOM_REAL; rhobb_kl(:,:,:) = 0.0_CUSTOM_REAL
  rhofb_kl(:,:,:) = 0.0_CUSTOM_REAL; rhofbb_kl(:,:,:) = 0.0_CUSTOM_REAL

  allocate(cpI_kl(NGLLX,NGLLZ,b_nspec_poroelastic), &
           cpII_kl(NGLLX,NGLLZ,b_nspec_poroelastic), &
           cs_kl(NGLLX,NGLLZ,b_nspec_poroelastic), &
           ratio_kl(NGLLX,NGLLZ,b_nspec_poroelastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating cpI_kl kernel arrays')
  cpI_kl(:,:,:) = 0.0_CUSTOM_REAL; cpII_kl(:,:,:) = 0.0_CUSTOM_REAL
  cs_kl(:,:,:) = 0.0_CUSTOM_REAL; ratio_kl(:,:,:) = 0.0_CUSTOM_REAL

  if (COMPUTE_INTEGRATED_ENERGY_FIELD) then
    allocate(total_integrated_energy_field(nspec), &
             max_total_energy_field(nspec), &
             total_effective_duration_field(nspec), &
             integrated_kinetic_energy_field(nspec), &
             max_kinetic_energy_field(nspec), &
             integrated_potential_energy_field(nspec), &
             max_potential_energy_field(nspec), &
             kinetic_effective_duration_field(nspec), &
             potential_effective_duration_field(nspec),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating total_integrated_energy_field arrays')
    total_integrated_energy_field(:) = 0._CUSTOM_REAL; total_effective_duration_field(:) = 0._CUSTOM_REAL
    integrated_kinetic_energy_field(:) = 0._CUSTOM_REAL; integrated_potential_energy_field(:) = 0._CUSTOM_REAL
    kinetic_effective_duration_field(:) = 0._CUSTOM_REAL; potential_effective_duration_field(:) = 0._CUSTOM_REAL
    max_total_energy_field(:) = 0._CUSTOM_REAL
    max_kinetic_energy_field(:) = 0._CUSTOM_REAL
    max_potential_energy_field(:) = 0._CUSTOM_REAL
  endif

  ! iglob_is_forced array is used when USE_ENFORCE_FIELDS is .true. (it says if a GLL point is forced or not)
  ! needs full array since we have it checked in if-statements
  allocate(iglob_is_forced(nglob),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating iglob_is_forced array')
  iglob_is_forced(:) = .false.

  ! following arrays are only fully accessed when option USE_ENFORCE_FIELDS is .true.
  if (USE_ENFORCE_FIELDS) then
    allocate(acoustic_iglob_is_forced(nglob), &
             elastic_iglob_is_forced(nglob), &
             modeAmplitude(nglob_acoustic),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating forced arrays')
  else
    ! dummy
    allocate(acoustic_iglob_is_forced(1),elastic_iglob_is_forced(1),modeAmplitude(1))
  endif
  acoustic_iglob_is_forced(:) = .false.
  elastic_iglob_is_forced(:) = .false.
  modeAmplitude(:) = 0.0d0

  ! synchronizes all processes
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  done initialization'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_wavefields

