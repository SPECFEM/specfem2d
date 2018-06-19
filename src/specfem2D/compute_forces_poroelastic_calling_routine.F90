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

  subroutine compute_forces_poroelastic_main()

  use constants, only: USE_PORO_VISCOUS_DAMPING,ALPHA_LDDRK,BETA_LDDRK
  use specfem_par

  implicit none

  ! non-blocking MPI
  ! iphase: iphase = 1 is for computing outer elements (on MPI interface),
  !         iphase = 2 is for computing inner elements
  integer :: iphase

  ! checks if anything to do in this slice
  if (.not. any_poroelastic) return

  ! implement viscous attenuation for poroelastic media
  if (ATTENUATION_PORO_FLUID_PART) call compute_attenuation_poro_fluid_part()

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase = 1,2

    ! main solver for the poroelastic elements: first the solid (u_s) then the fluid (w)
    call compute_forces_poro_solid(displs_poroelastic,displw_poroelastic,accels_poroelastic,epsilondev_s,iphase)
    call compute_forces_poro_fluid(displs_poroelastic,displw_poroelastic,accelw_poroelastic,epsilondev_w,iphase)

    ! compute additional contributions
    if (iphase == 1) then
      ! viscous damping
      if (USE_PORO_VISCOUS_DAMPING) call compute_forces_poro_viscous_damping()

      ! Stacey absorbing boundary
      if (anyabs) then
        call compute_stacey_poro_fluid(f0_source(1))
        call compute_stacey_poro_solid(f0_source(1))
      endif

      ! add coupling with the acoustic side
      if (coupled_acoustic_poro) then
        call compute_coupling_poro_ac(potential_dot_dot_acoustic,accels_poroelastic,accelw_poroelastic,1)
      endif

      ! add coupling with the elastic side
      if (coupled_elastic_poro) then
        call compute_coupling_poro_viscoelastic(displ_elastic,displs_poroelastic,displw_poroelastic,accels_poroelastic)
      endif

      ! add source
      if (.not. initialfield) then
        if (SIMULATION_TYPE == 1) then
          ! forward wavefield
          call compute_add_sources_poro(accels_poroelastic,accelw_poroelastic,it,i_stage)
        else if (SIMULATION_TYPE == 3) then
          ! adjoint wavefield
          call compute_add_sources_poro_adjoint()
        endif
      endif
    endif

#ifdef USE_MPI
    ! assembling accels_proelastic & accelw_poroelastic for poroelastic elements
    if (NPROC > 1 .and. ninterface_poroelastic > 0) then
      if (iphase == 1) then
        call assemble_MPI_vector_po_s(accels_poroelastic,accelw_poroelastic)
      else
        call assemble_MPI_vector_po_w(accels_poroelastic,accelw_poroelastic)
      endif
    endif
#endif

  enddo ! iphase

  ! multiply by the inverse of the mass matrix
  ! solid
  accels_poroelastic(1,:) = accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
  accels_poroelastic(2,:) = accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)
  ! fluid
  accelw_poroelastic(1,:) = accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
  accelw_poroelastic(2,:) = accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)

  ! update velocity
  select case (time_stepping_scheme)
  case (1)
    ! Newmark
    call update_veloc_poroelastic_Newmark()
  case (2)
    ! LDDRK
    call update_veloc_poroelastic_LDDRK()
  case (3)
    ! Runge-Kutta
    call update_veloc_poroelastic_RK()
  case default
    call stop_the_code('Time stepping scheme not implemented yet for poroelastic case')
  end select

  ! imposes continuity for stabilization
  if (coupled_elastic_poro) then
    call compute_coupling_poro_viscoelastic_for_stabilization(veloc_elastic,accel_elastic, &
                                                              velocs_poroelastic,accels_poroelastic, &
                                                              velocw_poroelastic,accelw_poroelastic,deltatover2)
  endif

  end subroutine compute_forces_poroelastic_main


!
!-------------------------------------------------------------------------------------
!

  subroutine compute_forces_poroelastic_main_backward()

  use constants, only: USE_PORO_VISCOUS_DAMPING,ALPHA_LDDRK,BETA_LDDRK
  use specfem_par

  implicit none

  ! non-blocking MPI
  ! iphase: iphase = 1 is for computing outer elements (on MPI interface),
  !         iphase = 2 is for computing inner elements
  integer :: iphase

  ! safety check
  if (SIMULATION_TYPE /= 3) return

  ! checks if anything to do in this slice
  if (.not. any_poroelastic) return

  ! implement viscous attenuation for poroelastic media
  if (ATTENUATION_PORO_FLUID_PART) call stop_the_code( &
'ATTENUATION_PORO_FLUID_PART not implemented yet for backward/kernel simulations')

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase = 1,2

    ! main solver for the poroelastic elements: first the solid (u_s) then the fluid (w)
    call compute_forces_poro_solid(b_displs_poroelastic,b_displw_poroelastic,b_accels_poroelastic,b_epsilondev_s,iphase)
    call compute_forces_poro_fluid(b_displs_poroelastic,b_displw_poroelastic,b_accelw_poroelastic,b_epsilondev_w,iphase)

    ! compute additional contributions
    if (iphase == 1) then
      ! viscous damping
      if (USE_PORO_VISCOUS_DAMPING) call compute_forces_poro_viscous_damping_backward()

      ! Stacey absorbing boundary
      if (anyabs) then
        call compute_stacey_poro_fluid_backward()
        call compute_stacey_poro_solid_backward()
      endif

      ! add coupling with the acoustic side
      if (coupled_acoustic_poro) then
        call compute_coupling_poro_ac(b_potential_dot_dot_acoustic,b_accels_poroelastic,b_accelw_poroelastic,3)
      endif

      ! add coupling with the elastic side
      if (coupled_elastic_poro) then
        call compute_coupling_poro_viscoelastic(b_displ_elastic,b_displs_poroelastic,b_displw_poroelastic,b_accels_poroelastic)
      endif

      ! add source
      if (.not. initialfield) then
        ! backward wavefield
        call compute_add_sources_poro(b_accels_poroelastic,b_accelw_poroelastic,NSTEP-it+1,stage_time_scheme-i_stage+1)
      endif
    endif

#ifdef USE_MPI
    ! assembling accels_proelastic & accelw_poroelastic for poroelastic elements
    if (NPROC > 1 .and. ninterface_poroelastic > 0) then
      if (iphase == 1) then
        call assemble_MPI_vector_po_s(b_accels_poroelastic,b_accelw_poroelastic)
      else
        call assemble_MPI_vector_po_w(b_accels_poroelastic,b_accelw_poroelastic)
      endif
    endif
#endif

  enddo ! iphase

  ! multiply by the inverse of the mass matrix
  ! solid
  b_accels_poroelastic(1,:) = b_accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
  b_accels_poroelastic(2,:) = b_accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)
  ! fluid
  b_accelw_poroelastic(1,:) = b_accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
  b_accelw_poroelastic(2,:) = b_accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)

  ! update velocity
  select case (time_stepping_scheme)
  case (1)
    ! Newmark
    call update_veloc_poroelastic_Newmark_backward()
  case (2)
    ! LDDRK
    if (SIMULATION_TYPE == 3) call stop_the_code('LDDRK scheme for poroelastic kernel simulation not implemented yet')
  case (3)
    ! Runge-Kutta
    if (SIMULATION_TYPE == 3) call stop_the_code('RK scheme for poroelastic kernel simulation not implemented yet')
  case default
    call stop_the_code('Time stepping scheme not implemented yet for poroelastic case')
  end select

  ! imposes continuity for stabilization
  if (coupled_elastic_poro) then
    call compute_coupling_poro_viscoelastic_for_stabilization(b_veloc_elastic,b_accel_elastic, &
                                                              b_velocs_poroelastic,b_accels_poroelastic, &
                                                              b_velocw_poroelastic,b_accelw_poroelastic,b_deltatover2)
  endif

  end subroutine compute_forces_poroelastic_main_backward

