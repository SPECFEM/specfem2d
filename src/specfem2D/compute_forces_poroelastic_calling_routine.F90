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

  subroutine compute_forces_poroelastic_main()

  use specfem_par

  implicit none

  ! local parameters
  ! for rk44
  double precision :: weight_rk

  ! checks if anything to do in this slice
  if (.not. any_poroelastic) return

  ! implement viscous attenuation for poroelastic media
  if (ATTENUATION_PORO_FLUID_PART) call compute_attenuation_poro_fluid_part()

  ! main solver for the poroelastic elements: first the solid (u_s) then the fluid (w)
  call compute_forces_poro_solid()
  call compute_forces_poro_fluid()

  ! viscous damping
  if (USE_PORO_VISCOUS_DAMPING) call compute_forces_poro_viscous_damping()

  ! Stacey absorbing boundary
  if (anyabs) then
    call compute_stacey_poro_fluid(f0_source(1))
    call compute_stacey_poro_solid(f0_source(1))
  endif

  ! add coupling with the acoustic side
  if (coupled_acoustic_poro) then
    call compute_coupling_poro_ac()
  endif

  ! add coupling with the elastic side
  if (coupled_elastic_poro) then
    call compute_coupling_poro_viscoelastic()
  endif

  ! add source
  if (.not. initialfield) then
    if (SIMULATION_TYPE == 1) then
      ! forward wavefield
      call compute_add_sources_poro(accels_poroelastic,accelw_poroelastic,it,i_stage)
    else if (SIMULATION_TYPE == 3) then
      ! backward wavefield
      call compute_add_sources_poro(b_accels_poroelastic,b_accelw_poroelastic,NSTEP-it+1,stage_time_scheme-i_stage+1)
      ! adjoint wavefield
      call compute_add_sources_poro_adjoint()
    endif
  endif

  ! assembling accels_proelastic & accelw_poroelastic for poroelastic elements
#ifdef USE_MPI
  if (NPROC > 1 .and. ninterface_poroelastic > 0) then
    call assemble_MPI_vector_po(accels_poroelastic,accelw_poroelastic)
    if (SIMULATION_TYPE == 3) then
      call assemble_MPI_vector_po(b_accels_poroelastic,b_accelw_poroelastic)
    endif
  endif
#endif

  ! multiply by the inverse of the mass matrix and update velocity
  if (time_stepping_scheme == 1) then
    accels_poroelastic(1,:) = accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
    accels_poroelastic(2,:) = accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)
    velocs_poroelastic = velocs_poroelastic + deltatover2*accels_poroelastic

    accelw_poroelastic(1,:) = accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
    accelw_poroelastic(2,:) = accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)
    velocw_poroelastic = velocw_poroelastic + deltatover2*accelw_poroelastic

    if (SIMULATION_TYPE == 3) then
      b_accels_poroelastic(1,:) = b_accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
      b_accels_poroelastic(2,:) = b_accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)
      b_velocs_poroelastic = b_velocs_poroelastic + b_deltatover2*b_accels_poroelastic

      b_accelw_poroelastic(1,:) = b_accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
      b_accelw_poroelastic(2,:) = b_accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)
      b_velocw_poroelastic = b_velocw_poroelastic + b_deltatover2*b_accelw_poroelastic
    endif
  endif

  if (time_stepping_scheme == 2) then
    accels_poroelastic(1,:) = accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
    accels_poroelastic(2,:) = accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)

    velocs_poroelastic_LDDRK = ALPHA_LDDRK(i_stage) * velocs_poroelastic_LDDRK + deltat * accels_poroelastic
    displs_poroelastic_LDDRK = ALPHA_LDDRK(i_stage) * displs_poroelastic_LDDRK + deltat * velocs_poroelastic
    velocs_poroelastic = velocs_poroelastic + BETA_LDDRK(i_stage) * velocs_poroelastic_LDDRK
    displs_poroelastic = displs_poroelastic + BETA_LDDRK(i_stage) * displs_poroelastic_LDDRK

    accelw_poroelastic(1,:) = accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
    accelw_poroelastic(2,:) = accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)

    velocw_poroelastic_LDDRK = ALPHA_LDDRK(i_stage) * velocw_poroelastic_LDDRK + deltat * accelw_poroelastic
    displw_poroelastic_LDDRK = ALPHA_LDDRK(i_stage) * displw_poroelastic_LDDRK + deltat * velocw_poroelastic
    velocw_poroelastic = velocw_poroelastic + BETA_LDDRK(i_stage) * velocw_poroelastic_LDDRK
    displw_poroelastic = displw_poroelastic + BETA_LDDRK(i_stage) * displw_poroelastic_LDDRK
  endif

  if (time_stepping_scheme == 3) then
    accels_poroelastic(1,:) = accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
    accels_poroelastic(2,:) = accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)

    accels_poroelastic_rk(1,:,i_stage) = deltat * accels_poroelastic(1,:)
    accels_poroelastic_rk(2,:,i_stage) = deltat * accels_poroelastic(2,:)
    velocs_poroelastic_rk(1,:,i_stage) = deltat * velocs_poroelastic(1,:)
    velocs_poroelastic_rk(2,:,i_stage) = deltat * velocs_poroelastic(2,:)

    accelw_poroelastic(1,:) = accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
    accelw_poroelastic(2,:) = accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)

    accelw_poroelastic_rk(1,:,i_stage) = deltat * accelw_poroelastic(1,:)
    accelw_poroelastic_rk(2,:,i_stage) = deltat * accelw_poroelastic(2,:)
    velocw_poroelastic_rk(1,:,i_stage) = deltat * velocw_poroelastic(1,:)
    velocw_poroelastic_rk(2,:,i_stage) = deltat * velocw_poroelastic(2,:)

    if (i_stage==1 .or. i_stage==2 .or. i_stage==3) then
      if (i_stage == 1) weight_rk = 0.5d0
      if (i_stage == 2) weight_rk = 0.5d0
      if (i_stage == 3) weight_rk = 1.0d0

      if (i_stage==1) then
        velocs_poroelastic_initial_rk(1,:) = velocs_poroelastic(1,:)
        velocs_poroelastic_initial_rk(2,:) = velocs_poroelastic(2,:)
        displs_poroelastic_initial_rk(1,:) = displs_poroelastic(1,:)
        displs_poroelastic_initial_rk(2,:) = displs_poroelastic(2,:)

        velocw_poroelastic_initial_rk(1,:) = velocw_poroelastic(1,:)
        velocw_poroelastic_initial_rk(2,:) = velocw_poroelastic(2,:)
        displw_poroelastic_initial_rk(1,:) = displw_poroelastic(1,:)
        displw_poroelastic_initial_rk(2,:) = displw_poroelastic(2,:)
      endif

      velocs_poroelastic(1,:) = velocs_poroelastic_initial_rk(1,:) + weight_rk * accels_poroelastic_rk(1,:,i_stage)
      velocs_poroelastic(2,:) = velocs_poroelastic_initial_rk(2,:) + weight_rk * accels_poroelastic_rk(2,:,i_stage)
      displs_poroelastic(1,:) = displs_poroelastic_initial_rk(1,:) + weight_rk * velocs_poroelastic_rk(1,:,i_stage)
      displs_poroelastic(2,:) = displs_poroelastic_initial_rk(2,:) + weight_rk * velocs_poroelastic_rk(2,:,i_stage)

      velocw_poroelastic(1,:) = velocw_poroelastic_initial_rk(1,:) + weight_rk * accelw_poroelastic_rk(1,:,i_stage)
      velocw_poroelastic(2,:) = velocw_poroelastic_initial_rk(2,:) + weight_rk * accelw_poroelastic_rk(2,:,i_stage)
      displw_poroelastic(1,:) = displw_poroelastic_initial_rk(1,:) + weight_rk * velocw_poroelastic_rk(1,:,i_stage)
      displw_poroelastic(2,:) = displw_poroelastic_initial_rk(2,:) + weight_rk * velocw_poroelastic_rk(2,:,i_stage)

    else if (i_stage==4) then

      velocs_poroelastic(1,:) = velocs_poroelastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
      (accels_poroelastic_rk(1,:,1) + 2.0d0 * accels_poroelastic_rk(1,:,2) + &
      2.0d0 * accels_poroelastic_rk(1,:,3) + accels_poroelastic_rk(1,:,4))

      velocs_poroelastic(2,:) = velocs_poroelastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
      (accels_poroelastic_rk(2,:,1) + 2.0d0 * accels_poroelastic_rk(2,:,2) + &
       2.0d0 * accels_poroelastic_rk(2,:,3) + accels_poroelastic_rk(2,:,4))

      displs_poroelastic(1,:) = displs_poroelastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
      (velocs_poroelastic_rk(1,:,1) + 2.0d0 * velocs_poroelastic_rk(1,:,2) + &
       2.0d0 * velocs_poroelastic_rk(1,:,3) + velocs_poroelastic_rk(1,:,4))

      displs_poroelastic(2,:) = displs_poroelastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
      (velocs_poroelastic_rk(2,:,1) + 2.0d0 * velocs_poroelastic_rk(2,:,2) + &
       2.0d0 * velocs_poroelastic_rk(2,:,3) + velocs_poroelastic_rk(2,:,4))

      velocw_poroelastic(1,:) = velocw_poroelastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
      (accelw_poroelastic_rk(1,:,1) + 2.0d0 * accelw_poroelastic_rk(1,:,2) + &
       2.0d0 * accelw_poroelastic_rk(1,:,3) + accelw_poroelastic_rk(1,:,4))

      velocw_poroelastic(2,:) = velocw_poroelastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
      (accelw_poroelastic_rk(2,:,1) + 2.0d0 * accelw_poroelastic_rk(2,:,2) + &
       2.0d0 * accelw_poroelastic_rk(2,:,3) + accelw_poroelastic_rk(2,:,4))

      displw_poroelastic(1,:) = displw_poroelastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
      (velocw_poroelastic_rk(1,:,1) + 2.0d0 * velocw_poroelastic_rk(1,:,2) + &
       2.0d0 * velocw_poroelastic_rk(1,:,3) + velocw_poroelastic_rk(1,:,4))

      displw_poroelastic(2,:) = displw_poroelastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
      (velocw_poroelastic_rk(2,:,1) + 2.0d0 * velocw_poroelastic_rk(2,:,2) + &
       2.0d0 * velocw_poroelastic_rk(2,:,3) + velocw_poroelastic_rk(2,:,4))
    endif
  endif

  ! Explanation of the code below, from Christina Morency and Yang Luo, January 2012:
  !
  ! Coupled elastic-poroelastic simulations imply continuity of traction and
  ! displacement at the interface.
  ! For the traction we pass on both sides n*(T + Te)/2 , that is, the average
  ! between the total stress (from the poroelastic part) and the elastic stress.
  ! For the displacement, we enforce its continuity in the assembling stage,
  ! realizing that continuity of displacement correspond to the continuity of
  ! the acceleration we have:
  !
  ! accel_elastic = rmass_inverse_elastic * force_elastic
  ! accels_poroelastic = rmass_s_inverse_poroelastic * force_poroelastic
  !
  ! Therefore, continuity of acceleration gives
  !
  ! accel = (force_elastic + force_poroelastic)/
  !     (1/rmass_inverse_elastic + 1/rmass_inverse_poroelastic)
  !
  ! Then
  !
  ! accel_elastic = accel
  ! accels_poroelastic = accel
  ! accelw_poroelastic = 0
  !
  ! From there, the velocity and displacement are updated.
  ! Note that force_elastic and force_poroelastic are the right hand sides of
  ! the equations we solve, that is, the acceleration terms before the
  ! division by the inverse of the mass matrices. This is why in the code below
  ! we first need to recover the accelerations (which are then
  ! the right hand sides forces) and the velocities before the update.
  !
  ! This implementation highly helped stability especially with unstructured meshes.
  !
  if (coupled_elastic_poro) then
    call compute_coupling_poro_viscoelastic_for_stabilization()
  endif

  end subroutine compute_forces_poroelastic_main

