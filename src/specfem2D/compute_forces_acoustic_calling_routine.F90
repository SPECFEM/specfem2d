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

  subroutine compute_forces_acoustic_main()

  use specfem_par

  implicit none

  ! local parameters
  integer :: i
  ! for rk44
  double precision :: weight_rk

  ! main solver for the acoustic elements
  if (any_acoustic) then
    call compute_forces_acoustic(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                 potential_acoustic,potential_acoustic_old,PML_BOUNDARY_CONDITIONS)
  endif

  ! add acoustic forcing at a rigid boundary
  if (any_acoustic) then
    if (ACOUSTIC_FORCING) then
      call add_acoustic_forcing_at_rigid_boundary(potential_dot_dot_acoustic)
    endif
  endif

  ! add coupling with the elastic side
  if (coupled_acoustic_elastic) then
    if (SIMULATION_TYPE == 1) then
      call compute_coupling_acoustic_el(displ_elastic,displ_elastic_old,potential_dot_dot_acoustic)
    endif

    ! coupling for adjoint wavefields
    if (SIMULATION_TYPE == 3) then
      ! note: handles adjoint runs coupling between adjoint potential and adjoint elastic wavefield
      !       adjoint definition: \partial_t^2 \bfs^\dagger = - \frac{1}{\rho} \bfnabla \phi^\dagger
      !
      ! coupling with adjoint wavefields
      call compute_coupling_acoustic_el(accel_elastic_adj_coupling,displ_elastic_old,potential_dot_dot_acoustic)
    endif
  endif

  ! add coupling with the poroelastic side
  if (coupled_acoustic_poro) then
    call compute_coupling_acoustic_po()
  endif

  ! add force source
  if (any_acoustic) then
    if (.not. initialfield) then
      if (SIMULATION_TYPE == 1) then
        call compute_add_sources_acoustic(potential_dot_dot_acoustic,it,i_stage)
      else if (SIMULATION_TYPE == 3) then
        ! adjoint sources
        call compute_add_sources_acoustic_adjoint()
      endif
    endif
  endif

  ! assembling potential_dot_dot or b_potential_dot_dot for acoustic elements
#ifdef USE_MPI
  if (NPROC > 1 .and. any_acoustic .and. ninterface_acoustic > 0) then
    call assemble_MPI_vector_ac(potential_dot_dot_acoustic)

    if (time_stepping_scheme == 2) then
      if (i_stage==1 .and. it == 1 .and. (.not. initialfield)) then
        potential_dot_acoustic_temp(:) = potential_dot_acoustic(:)
        call assemble_MPI_vector_ac(potential_dot_acoustic)
      endif
    endif
  endif
#endif

  if (PML_BOUNDARY_CONDITIONS) then
    if (any_acoustic .and. nglob_interface > 0) then
      if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
        do i = 1, nglob_interface
          write(72)potential_dot_dot_acoustic(point_interface(i)),&
                   potential_dot_acoustic(point_interface(i)),&
                   potential_acoustic(point_interface(i))
        enddo
      endif
    endif
  endif

  ! multiply by the inverse of the mass matrix and update velocity
  if (any_acoustic) then
    ! free surface for an acoustic medium
    if (nelem_acoustic_surface > 0) then
      call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                         potential_acoustic)
    endif

    if (time_stepping_scheme == 1) then

      !! DK DK this should be vectorized
      potential_dot_dot_acoustic(:) = potential_dot_dot_acoustic(:) * rmass_inverse_acoustic(:)

      potential_dot_acoustic(:) = potential_dot_acoustic(:) + deltatover2 * potential_dot_dot_acoustic(:)

      ! update the potential field (use a new array here) for coupling terms
      potential_acoustic_adj_coupling(:) = potential_acoustic(:) + deltat * potential_dot_acoustic(:) + &
                                        deltatsquareover2 * potential_dot_dot_acoustic(:)
    endif

    if (time_stepping_scheme == 2) then
      !! DK DK this should be vectorized
      potential_dot_dot_acoustic(:) = potential_dot_dot_acoustic(:) * rmass_inverse_acoustic(:)

      potential_dot_acoustic_LDDRK(:) = alpha_LDDRK(i_stage) * potential_dot_acoustic_LDDRK(:) + &
                                        deltat * potential_dot_dot_acoustic(:)
      potential_acoustic_LDDRK(:) = alpha_LDDRK(i_stage) * potential_acoustic_LDDRK(:) + &
                                    deltat * potential_dot_acoustic(:)

      if (i_stage==1 .and. it == 1 .and. (.not. initialfield)) then
        !! DK DK this should be vectorized
        potential_dot_acoustic_temp(:) = potential_dot_acoustic_temp(:) + &
                                         beta_LDDRK(i_stage) * potential_dot_acoustic_LDDRK(:)
        potential_dot_acoustic(:) = potential_dot_acoustic_temp(:)
      else
        potential_dot_acoustic(:) = potential_dot_acoustic(:) + beta_LDDRK(i_stage) * potential_dot_acoustic_LDDRK(:)
      endif

      !! DK DK this should be vectorized
      potential_acoustic(:) = potential_acoustic(:) + beta_LDDRK(i_stage) * potential_acoustic_LDDRK(:)
    endif

    if (time_stepping_scheme == 3) then
      !! DK DK this should be vectorized
      potential_dot_dot_acoustic(:) = potential_dot_dot_acoustic(:) * rmass_inverse_acoustic(:)

      potential_dot_dot_acoustic_rk(:,i_stage) = deltat * potential_dot_dot_acoustic(:)
      potential_dot_acoustic_rk(:,i_stage) = deltat * potential_dot_acoustic(:)

      if (i_stage == 1 .or. i_stage == 2 .or. i_stage == 3) then
        if (i_stage == 1) weight_rk = 0.5d0
        if (i_stage == 2) weight_rk = 0.5d0
        if (i_stage == 3) weight_rk = 1.0d0

        if (i_stage == 1) then
!! DK DK this should be vectorized
          potential_dot_acoustic_init_rk(:) = potential_dot_acoustic(:)
          potential_acoustic_init_rk(:) = potential_acoustic(:)
        endif
!! DK DK this should be vectorized
        potential_dot_acoustic(:) = potential_dot_acoustic_init_rk(:) + &
                                    weight_rk * potential_dot_dot_acoustic_rk(:,i_stage)
        potential_acoustic(:) = potential_acoustic_init_rk(:) + weight_rk * potential_dot_acoustic_rk(:,i_stage)
      else if (i_stage == 4) then
!! DK DK this should be vectorized
        potential_dot_acoustic(:) = potential_dot_acoustic_init_rk(:) + &
                                    1.0d0 / 6.0d0 * ( potential_dot_dot_acoustic_rk(:,1) + &
                                                      2.0d0 * potential_dot_dot_acoustic_rk(:,2) + &
                                                      2.0d0 * potential_dot_dot_acoustic_rk(:,3) + &
                                                      potential_dot_dot_acoustic_rk(:,4) )

!! DK DK this should be vectorized
        potential_acoustic(:) = potential_acoustic_init_rk(:) + &
                                1.0d0 / 6.0d0 * ( potential_dot_acoustic_rk(:,1) + &
                                                  2.0d0 * potential_dot_acoustic_rk(:,2) + &
                                                  2.0d0 * potential_dot_acoustic_rk(:,3) + &
                                                  potential_dot_acoustic_rk(:,4) )
      endif
    endif
  endif ! of if (any_acoustic)

  end subroutine compute_forces_acoustic_main

!
!-------------------------------------------------------------------------------------
!

  subroutine compute_forces_acoustic_main_backward()

  use specfem_par

  implicit none

  ! local parameters
  integer :: i,it_temp,istage_temp

  ! checks
  if (SIMULATION_TYPE /= 3 ) return

  ! timing
  if (UNDO_ATTENUATION) then
    it_temp = NSTEP - (iteration_on_subset * NT_DUMP_ATTENUATION - it_of_this_subset + 1)
    istage_temp = i_stage
  else
    it_temp = NSTEP - it + 1
    istage_temp = stage_time_scheme - i_stage + 1
  endif

  ! main solver for the acoustic elements
  if (any_acoustic) then
    if (PML_BOUNDARY_CONDITIONS) then
      call rebuild_value_on_PML_interface_acoustic(it_temp)
    endif

    call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                       b_potential_acoustic)

    if (UNDO_ATTENUATION) then
      call compute_forces_acoustic(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                   b_potential_acoustic,b_potential_acoustic_old,.false.)

    else
      call compute_forces_acoustic_backward(b_potential_dot_dot_acoustic,b_potential_acoustic)
    endif

    if (PML_BOUNDARY_CONDITIONS) then
      call rebuild_value_on_PML_interface_acoustic(it_temp)
    endif
  endif ! end of test if any acoustic element

  ! add acoustic forcing at a rigid boundary
  if (any_acoustic) then
    if (ACOUSTIC_FORCING) then
      call add_acoustic_forcing_at_rigid_boundary(b_potential_dot_dot_acoustic)
    endif
  endif

  ! add coupling with the elastic side
  if (coupled_acoustic_elastic) then
    call compute_coupling_acoustic_el_backward(b_displ_elastic,b_potential_dot_dot_acoustic)
  endif

  ! add coupling with the poroelastic side
  if (coupled_acoustic_poro) then
    call compute_coupling_acoustic_po_backward()
  endif

  if (PML_BOUNDARY_CONDITIONS) then
    call rebuild_value_on_PML_interface_acoustic(it_temp)
  endif

  ! add force source
  if (any_acoustic) then
    if (.not. initialfield) then
      ! backward wavefield
      call compute_add_sources_acoustic(b_potential_dot_dot_acoustic,it_temp,istage_temp)
    endif
  endif

  ! assembling potential_dot_dot or b_potential_dot_dot for acoustic elements
#ifdef USE_MPI
  if (NPROC > 1 .and. any_acoustic .and. ninterface_acoustic > 0) then
    call assemble_MPI_vector_ac(b_potential_dot_dot_acoustic)
  endif
#endif

  if (PML_BOUNDARY_CONDITIONS) then
    if (any_acoustic .and. nglob_interface > 0) then
      do i = 1, nglob_interface
        b_potential_dot_dot_acoustic(point_interface(i)) = pml_interface_history_potential_dot_dot(i,it_temp)
      enddo
    endif
  endif

  ! multiply by the inverse of the mass matrix and update velocity
  if (any_acoustic) then
    ! free surface for an acoustic medium
    if (nelem_acoustic_surface > 0) then
      call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                         b_potential_acoustic)
    endif

    if (time_stepping_scheme == 1) then
      !! DK DK this should be vectorized
      b_potential_dot_dot_acoustic(:) = b_potential_dot_dot_acoustic(:) * rmass_inverse_acoustic(:)

      b_potential_dot_acoustic(:) = b_potential_dot_acoustic(:) + b_deltatover2 * b_potential_dot_dot_acoustic(:)
    endif

  endif ! of if (any_acoustic)

  end subroutine compute_forces_acoustic_main_backward



