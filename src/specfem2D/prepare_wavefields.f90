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

  use constants, only: IMAIN,APPROXIMATE_HESS_KL
  use specfem_par

  implicit none

  ! local parameters
  integer :: nglob_acoustic_b,nglob_elastic_b,nglob_poroelastic_b
  integer :: ier

  ! displacement, velocity, acceleration and inverse of the mass matrix for elastic elements
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Preparing array allocations'
    call flush_IMAIN()
  endif

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

  ! PML
  allocate(displ_elastic_old(NDIM,nglob_elastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating old elastic wavefield arrays')

  if (SIMULATION_TYPE == 3) then
    if (coupled_acoustic_elastic) then
      allocate(accel_elastic_adj_coupling(NDIM,nglob_elastic))
    endif
  endif

  allocate(rmass_inverse_elastic(NDIM,nglob_elastic),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating elastic mass matrix array')

  if (time_stepping_scheme == 2) then
    allocate(displ_elastic_LDDRK(NDIM,nglob_elastic), &
             veloc_elastic_LDDRK(NDIM,nglob_elastic), &
             veloc_elastic_LDDRK_temp(NDIM,nglob_elastic),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating elastic LDDRK wavefield arrays')
  endif

  if (time_stepping_scheme == 3) then
    allocate(accel_elastic_rk(NDIM,nglob_elastic,stage_time_scheme), &
             veloc_elastic_rk(NDIM,nglob_elastic,stage_time_scheme), &
             veloc_elastic_initial_rk(NDIM,nglob_elastic), &
             displ_elastic_initial_rk(NDIM,nglob_elastic),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating elastic RK wavefield arrays')
  endif

  ! extra array if adjoint and kernels calculation
  if (SIMULATION_TYPE == 3 .and. any_elastic) then
    nglob_elastic_b = nglob
    nspec_elastic_b = nspec
  else
    ! dummy allocations
    nglob_elastic_b = 1
    nspec_elastic_b = 1
  endif

  allocate(b_displ_elastic(NDIM,nglob_elastic_b), &
           b_veloc_elastic(NDIM,nglob_elastic_b), &
           b_accel_elastic(NDIM,nglob_elastic_b),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating elastic backward wavefield arrays')

  allocate(b_displ_elastic_old(NDIM,nglob_elastic_b))

  ! kernels
  ! on global nodes
  allocate(rho_k(nglob_elastic_b))
  allocate(mu_k(nglob_elastic_b))
  allocate(kappa_k(nglob_elastic_b))
  allocate(c11_k(nglob_elastic_b))
  allocate(c13_k(nglob_elastic_b))
  allocate(c15_k(nglob_elastic_b))
  allocate(c33_k(nglob_elastic_b))
  allocate(c35_k(nglob_elastic_b))
  allocate(c55_k(nglob_elastic_b))
  ! on local nodes
  allocate(rho_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(mu_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(kappa_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(rhop_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(alpha_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(beta_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(bulk_c_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(bulk_beta_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(c11_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(c13_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(c15_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(c33_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(c35_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(c55_kl(NGLLX,NGLLZ,nspec_elastic_b))
  if (APPROXIMATE_HESS_KL) then
    allocate(rhorho_el_Hessian_final2(NGLLX,NGLLZ,nspec_elastic_b))
    allocate(rhorho_el_Hessian_final1(NGLLX,NGLLZ,nspec_elastic_b))
  endif

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
  allocate(displs_poroelastic(NDIM,nglob_poroelastic))
  allocate(displs_poroelastic_old(NDIM,nglob_poroelastic))
  allocate(velocs_poroelastic(NDIM,nglob_poroelastic))
  allocate(accels_poroelastic(NDIM,nglob_poroelastic))
  if (SIMULATION_TYPE == 3) then
    allocate(accels_poroelastic_adj_coupling(NDIM,nglob_poroelastic))
  endif
  allocate(rmass_s_inverse_poroelastic(nglob_poroelastic))

  allocate(displw_poroelastic(NDIM,nglob_poroelastic))
  allocate(velocw_poroelastic(NDIM,nglob_poroelastic))
  allocate(accelw_poroelastic(NDIM,nglob_poroelastic))
  if (SIMULATION_TYPE == 3) then
    allocate(accelw_poroelastic_adj_coupling(NDIM,nglob_poroelastic))
  endif
  allocate(rmass_w_inverse_poroelastic(nglob_poroelastic))

  if (time_stepping_scheme == 2) then
    allocate(displs_poroelastic_LDDRK(NDIM,nglob_poroelastic))
    allocate(velocs_poroelastic_LDDRK(NDIM,nglob_poroelastic))
    allocate(displw_poroelastic_LDDRK(NDIM,nglob_poroelastic))
    allocate(velocw_poroelastic_LDDRK(NDIM,nglob_poroelastic))
  endif

  if (time_stepping_scheme == 3) then
    allocate(accels_poroelastic_rk(NDIM,nglob_poroelastic,stage_time_scheme))
    allocate(velocs_poroelastic_rk(NDIM,nglob_poroelastic,stage_time_scheme))
    allocate(accelw_poroelastic_rk(NDIM,nglob_poroelastic,stage_time_scheme))
    allocate(velocw_poroelastic_rk(NDIM,nglob_poroelastic,stage_time_scheme))
    allocate(displs_poroelastic_initial_rk(NDIM,nglob_poroelastic))
    allocate(velocs_poroelastic_initial_rk(NDIM,nglob_poroelastic))
    allocate(displw_poroelastic_initial_rk(NDIM,nglob_poroelastic))
    allocate(velocw_poroelastic_initial_rk(NDIM,nglob_poroelastic))
  endif

  ! extra array if adjoint and kernels calculation
  if (SIMULATION_TYPE == 3 .and. any_poroelastic) then
    nglob_poroelastic_b = nglob
    nspec_poroelastic_b = nspec
  else
    ! dummy allocations
    nglob_poroelastic_b = 1
    nspec_poroelastic_b = 1
  endif
  allocate(b_displs_poroelastic(NDIM,nglob_poroelastic_b))
  allocate(b_velocs_poroelastic(NDIM,nglob_poroelastic_b))
  allocate(b_accels_poroelastic(NDIM,nglob_poroelastic_b))
  allocate(b_displw_poroelastic(NDIM,nglob_poroelastic_b))
  allocate(b_velocw_poroelastic(NDIM,nglob_poroelastic_b))
  allocate(b_accelw_poroelastic(NDIM,nglob_poroelastic_b))
  ! strain
  allocate(epsilondev_s(4,NGLLX,NGLLZ,nspec_poroelastic_b), &
           epsilondev_w(4,NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(b_epsilondev_s(4,NGLLX,NGLLZ,nspec_poroelastic_b), &
           b_epsilondev_w(4,NGLLX,NGLLZ,nspec_poroelastic_b))
  ! kernels
  ! on global nodes
  allocate(rhot_k(nglob_poroelastic_b))
  allocate(rhof_k(nglob_poroelastic_b))
  allocate(sm_k(nglob_poroelastic_b))
  allocate(eta_k(nglob_poroelastic_b))
  allocate(mufr_k(nglob_poroelastic_b))
  allocate(B_k(nglob_poroelastic_b))
  allocate(C_k(nglob_poroelastic_b))
  allocate(M_k(nglob_poroelastic_b))
  ! on local nodes
  allocate(rhot_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(rhof_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(sm_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(eta_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(mufr_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(B_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(C_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(M_kl(NGLLX,NGLLZ,nspec_poroelastic_b))

  allocate(phi_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(phib_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(mufrb_kl(NGLLX,NGLLZ,nspec_poroelastic_b))

  allocate(rhob_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(rhobb_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(rhofb_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(rhofbb_kl(NGLLX,NGLLZ,nspec_poroelastic_b))

  allocate(cpI_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(cpII_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(cs_kl(NGLLX,NGLLZ,nspec_poroelastic_b))

  allocate(ratio_kl(NGLLX,NGLLZ,nspec_poroelastic_b))

  if (COMPUTE_INTEGRATED_ENERGY_FIELD) then
    allocate(total_integrated_energy_field(nspec),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating total_integrated_energy_field array')
    total_integrated_energy_field(:) = 0._CUSTOM_REAL
    allocate(max_total_energy_field(nspec),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating max_total_energy_field array')
    max_total_energy_field(:) = 0._CUSTOM_REAL
    allocate(total_effective_duration_field(nspec),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating total_effective_duration_field array')
    total_effective_duration_field(:) = 0._CUSTOM_REAL
    allocate(integrated_kinetic_energy_field(nspec),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating integrated_kinetic_energy_field array')
    integrated_kinetic_energy_field(:) = 0._CUSTOM_REAL
    allocate(max_kinetic_energy_field(nspec),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating max_kinetic_energy_field array')
    max_kinetic_energy_field(:) = 0._CUSTOM_REAL
    allocate(integrated_potential_energy_field(nspec),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating integrated_potential_energy_field array')
    integrated_potential_energy_field(:) = 0._CUSTOM_REAL
    allocate(max_potential_energy_field(nspec),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating max_potential_energy_field array')
    max_potential_energy_field(:) = 0._CUSTOM_REAL
    allocate(kinetic_effective_duration_field(nspec),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating kinetic_effective_duration_field array')
    kinetic_effective_duration_field(:) = 0._CUSTOM_REAL
    allocate(potential_effective_duration_field(nspec),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating potential_effective_duration_field array')
    potential_effective_duration_field(:) = 0._CUSTOM_REAL
  endif
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

  allocate(potential_acoustic(nglob_acoustic))
  allocate(potential_acoustic_old(nglob_acoustic))
  if (SIMULATION_TYPE == 3) then
    allocate(potential_acoustic_adj_coupling(nglob_acoustic))
  endif

  allocate(potential_dot_acoustic(nglob_acoustic))
  allocate(potential_dot_dot_acoustic(nglob_acoustic))

  if (time_stepping_scheme == 2) then
    allocate(potential_acoustic_LDDRK(nglob_acoustic))
    allocate(potential_dot_acoustic_LDDRK(nglob_acoustic))
    allocate(potential_dot_acoustic_temp(nglob_acoustic))
  endif

  if (time_stepping_scheme == 3) then
    allocate(potential_acoustic_init_rk(nglob_acoustic))
    allocate(potential_dot_acoustic_init_rk(nglob_acoustic))
    allocate(potential_dot_dot_acoustic_rk(nglob_acoustic,stage_time_scheme))
    allocate(potential_dot_acoustic_rk(nglob_acoustic,stage_time_scheme))
  endif

  allocate(rmass_inverse_acoustic(nglob_acoustic))
!! DK DK March 2018: this was missing in Quentin Brissaud's new variational formulation for viscoacoustic media; I added it
  if (ATTENUATION_VISCOACOUSTIC) then
    allocate(rmass_inverse_e1(nglob_acoustic,N_SLS))
  else
    allocate(rmass_inverse_e1(1,1))
  endif

  if (SIMULATION_TYPE == 3 .and. any_acoustic) then
    nglob_acoustic_b = nglob
    nspec_acoustic_b = nspec
  else
    ! dummy array allocations
    ! allocates unused arrays with fictitious size
    nglob_acoustic_b = 1
    nspec_acoustic_b = 1
  endif
  allocate(b_potential_acoustic(nglob_acoustic_b))
  allocate(b_potential_acoustic_old(nglob_acoustic_b))
  allocate(b_potential_dot_acoustic(nglob_acoustic_b))
  allocate(b_potential_dot_dot_acoustic(nglob_acoustic_b))
  allocate(b_displ_ac(2,nglob_acoustic_b))
  allocate(b_accel_ac(2,nglob_acoustic_b))
  allocate(accel_ac(2,nglob_acoustic_b))
  ! kernels
  ! on global points
  allocate(rhol_ac_global(nglob_acoustic_b))
  allocate(kappal_ac_global(nglob_acoustic_b))
  ! on local points
  allocate(rho_ac_kl(NGLLX,NGLLZ,nspec_acoustic_b))
  allocate(kappa_ac_kl(NGLLX,NGLLZ,nspec_acoustic_b))
  allocate(rhop_ac_kl(NGLLX,NGLLZ,nspec_acoustic_b))
  allocate(alpha_ac_kl(NGLLX,NGLLZ,nspec_acoustic_b))
  if (APPROXIMATE_HESS_KL) then
    allocate(rhorho_ac_Hessian_final2(NGLLX,NGLLZ,nspec_acoustic_b))
    allocate(rhorho_ac_Hessian_final1(NGLLX,NGLLZ,nspec_acoustic_b))
  endif

  ! iglob_is_forced array is used when USE_ENFORCE_FIELDS is .true. (it says if a GLL point is forced or not)
  allocate(iglob_is_forced(nglob))
  allocate(acoustic_iglob_is_forced(nglob))
  allocate(elastic_iglob_is_forced(nglob))
  allocate(modeAmplitude(nglob_acoustic))
  iglob_is_forced(:) = .false.
  acoustic_iglob_is_forced(:) = .false.
  elastic_iglob_is_forced(:) = .false.
  modeAmplitude(:) = 0.0d0

  ! synchronizes all processes
  call synchronize_all()

  ! initializes wavefields
  if (myrank == 0) then
    write(IMAIN,*) '  wavefield initialization'
    call flush_IMAIN()
  endif

  ! initialize arrays to zero
  displ_elastic(:,:) = 0._CUSTOM_REAL
  displ_elastic_old(:,:) = 0._CUSTOM_REAL
  veloc_elastic(:,:) = 0._CUSTOM_REAL
  accel_elastic(:,:) = 0._CUSTOM_REAL

  if (SIMULATION_TYPE == 3 .and. any_elastic) then
    b_displ_elastic_old(:,:) = 0._CUSTOM_REAL
    b_displ_elastic(:,:) = 0._CUSTOM_REAL
    b_veloc_elastic(:,:) = 0._CUSTOM_REAL
    b_accel_elastic(:,:) = 0._CUSTOM_REAL
  endif

  if (time_stepping_scheme == 2) then
    displ_elastic_LDDRK(:,:) = 0._CUSTOM_REAL
    veloc_elastic_LDDRK(:,:) = 0._CUSTOM_REAL
    veloc_elastic_LDDRK_temp(:,:) = 0._CUSTOM_REAL
  endif

  if (time_stepping_scheme == 3) then
    accel_elastic_rk(:,:,:) = 0._CUSTOM_REAL
    veloc_elastic_rk(:,:,:) = 0._CUSTOM_REAL
    veloc_elastic_initial_rk(:,:) = 0._CUSTOM_REAL
    displ_elastic_initial_rk(:,:) = 0._CUSTOM_REAL
  endif

  displs_poroelastic(:,:) = 0._CUSTOM_REAL
  displs_poroelastic_old(:,:) = 0._CUSTOM_REAL
  velocs_poroelastic(:,:) = 0._CUSTOM_REAL
  accels_poroelastic(:,:) = 0._CUSTOM_REAL
  displw_poroelastic(:,:) = 0._CUSTOM_REAL
  velocw_poroelastic(:,:) = 0._CUSTOM_REAL
  accelw_poroelastic(:,:) = 0._CUSTOM_REAL

  if (time_stepping_scheme == 2) then
    displs_poroelastic_LDDRK(:,:) = 0._CUSTOM_REAL
    velocs_poroelastic_LDDRK(:,:) = 0._CUSTOM_REAL
    displw_poroelastic_LDDRK(:,:) = 0._CUSTOM_REAL
    velocw_poroelastic_LDDRK(:,:) = 0._CUSTOM_REAL
  endif

  if (time_stepping_scheme == 3) then
    accels_poroelastic_rk(:,:,:) = 0._CUSTOM_REAL
    velocs_poroelastic_rk(:,:,:) = 0._CUSTOM_REAL

    accelw_poroelastic_rk(:,:,:) = 0._CUSTOM_REAL
    velocw_poroelastic_rk(:,:,:) = 0._CUSTOM_REAL

    velocs_poroelastic_initial_rk(:,:) = 0._CUSTOM_REAL
    displs_poroelastic_initial_rk(:,:) = 0._CUSTOM_REAL

    velocw_poroelastic_initial_rk(:,:) = 0._CUSTOM_REAL
    displw_poroelastic_initial_rk(:,:) = 0._CUSTOM_REAL
  endif

  potential_acoustic(:) = 0._CUSTOM_REAL
  potential_acoustic_old(:) = 0._CUSTOM_REAL
  potential_dot_acoustic(:) = 0._CUSTOM_REAL
  potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL

  b_potential_acoustic(:) = 0._CUSTOM_REAL
  b_potential_acoustic_old(:) = 0._CUSTOM_REAL
  b_potential_dot_acoustic(:) = 0._CUSTOM_REAL
  b_potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL

  if (time_stepping_scheme == 2) then
    potential_acoustic_LDDRK(:) = 0._CUSTOM_REAL
    potential_dot_acoustic_LDDRK(:) = 0._CUSTOM_REAL
    potential_dot_acoustic_temp(:) = 0._CUSTOM_REAL
  endif

  if (time_stepping_scheme == 3) then
    potential_acoustic_init_rk(:) = 0._CUSTOM_REAL
    potential_dot_acoustic_init_rk(:) = 0._CUSTOM_REAL
    potential_dot_dot_acoustic_rk(:,:) = 0._CUSTOM_REAL
    potential_dot_acoustic_rk(:,:) = 0._CUSTOM_REAL
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  done initialization'
    call flush_IMAIN()
  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_wavefields

