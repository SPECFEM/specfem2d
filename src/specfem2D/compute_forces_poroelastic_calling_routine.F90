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

  use constants,only: USE_PORO_VISCOUS_DAMPING,ALPHA_LDDRK,BETA_LDDRK
  use specfem_par

  implicit none

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

  ! multiply by the inverse of the mass matrix
  ! solid
  accels_poroelastic(1,:) = accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
  accels_poroelastic(2,:) = accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)
  ! fluid
  accelw_poroelastic(1,:) = accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
  accelw_poroelastic(2,:) = accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)

  if (SIMULATION_TYPE == 3) then
    b_accels_poroelastic(1,:) = b_accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
    b_accels_poroelastic(2,:) = b_accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)

    b_accelw_poroelastic(1,:) = b_accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
    b_accelw_poroelastic(2,:) = b_accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)
  endif

  ! update velocity
  select case (time_stepping_scheme)
  case (1)
    ! Newmark
    call update_veloc_poroelastic_Newmark()
    if (SIMULATION_TYPE == 3) call update_veloc_poroelastic_Newmark_backward()
  case (2)
    ! LDDRK
    call update_veloc_poroelastic_LDDRK()
    if (SIMULATION_TYPE == 3) stop 'LDDRK scheme for poroelastic kernel simulation not implemented yet'
  case (3)
    ! Runge-Kutta
    call update_veloc_poroelastic_RK()
    if (SIMULATION_TYPE == 3) stop 'RK scheme for poroelastic kernel simulation not implemented yet'
  case default
    stop 'Time stepping scheme not implemented yet for poroelastic case'
  end select

  ! imposes continuity for stabilization
  if (coupled_elastic_poro) then
    call compute_coupling_poro_viscoelastic_for_stabilization()
  endif

  end subroutine compute_forces_poroelastic_main

