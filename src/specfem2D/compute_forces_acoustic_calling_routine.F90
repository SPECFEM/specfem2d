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

  ! checks if anything to do in this slice
  if (.not. any_acoustic) return

  ! free surface for an acoustic medium
  call enforce_acoustic_free_surface(minus_pressure_acoustic,minus_int_pressure_acoustic,minus_int_int_pressure_acoustic)

  ! main solver for the acoustic elements
  call compute_forces_acoustic(minus_pressure_acoustic,minus_int_pressure_acoustic,minus_int_int_pressure_acoustic, &
                               PML_BOUNDARY_CONDITIONS,minus_int_int_pressure_acoustic_old)

  ! Stacey boundary conditions
  if (STACEY_BOUNDARY_CONDITIONS) then
    call compute_stacey_acoustic(minus_pressure_acoustic,minus_int_pressure_acoustic)
  endif

  ! PML boundary conditions
  if (PML_BOUNDARY_CONDITIONS) then
    call pml_boundary_acoustic(minus_pressure_acoustic,minus_int_pressure_acoustic, &
                               minus_int_int_pressure_acoustic,minus_int_int_pressure_acoustic_old)
  endif

  ! add acoustic forcing at a rigid boundary
  if (ACOUSTIC_FORCING) then
    call add_acoustic_forcing_at_rigid_boundary(minus_pressure_acoustic)
  endif

  ! applies to coupling in case of MPI partitioning:
  !   coupling interfaces might not be properly detected if one material domain is only in an another slice.
  !   in such a case, the common nodes would not be detected as belonging to a coupling interface.
  !   something to do in future...

  ! add coupling with the elastic side
  if (coupled_acoustic_elastic) then
    if (SIMULATION_TYPE == 1) then
      call compute_coupling_acoustic_el(displ_elastic,displ_elastic_old,minus_pressure_acoustic)
    endif

    ! coupling for adjoint wavefields
    if (SIMULATION_TYPE == 3) then
      ! note: handles adjoint runs coupling between adjoint pressure and adjoint elastic wavefield
      !       adjoint definition: \partial_t^2 \bfs^\dagger = - \frac{1}{\rho} \bfnabla \phi^\dagger
      !
      ! coupling with adjoint wavefields
      call compute_coupling_acoustic_el(accel_elastic_adj_coupling,displ_elastic_old,minus_pressure_acoustic)
    endif
  endif

  ! add coupling with the poroelastic side
  if (coupled_acoustic_poro) then
    call compute_coupling_acoustic_po()
  endif

  ! add force source
  if (.not. initialfield) then
    if (SIMULATION_TYPE == 1) then
      call compute_add_sources_acoustic(minus_pressure_acoustic,it,i_stage)
    else if (SIMULATION_TYPE == 3) then
      ! adjoint sources
      call compute_add_sources_acoustic_adjoint()
    endif
  endif

  ! assembling minus_pressure or b_minus_pressure for acoustic elements
#ifdef USE_MPI
  if (NPROC > 1 .and. ninterface_acoustic > 0) then
    call assemble_MPI_vector_ac(minus_pressure_acoustic)

    if (time_stepping_scheme == 2) then
      if (i_stage==1 .and. it == 1 .and. (.not. initialfield)) then
        minus_int_pressure_acoustic_temp(:) = minus_int_pressure_acoustic(:)
        call assemble_MPI_vector_ac(minus_int_pressure_acoustic)
      endif
    endif
  endif
#endif

  if (PML_BOUNDARY_CONDITIONS) then
    if (nglob_interface > 0) then
      if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
        do i = 1, nglob_interface
          write(72) minus_pressure_acoustic(point_interface(i)), minus_int_pressure_acoustic(point_interface(i)), &
                    minus_int_int_pressure_acoustic(point_interface(i))
        enddo
      endif
    endif
  endif

  ! multiply by the inverse of the mass matrix and update velocity
  if (time_stepping_scheme == 1) then
    ! Newmark scheme
    !! DK DK this should be vectorized
    minus_pressure_acoustic(:) = minus_pressure_acoustic(:) * rmass_inverse_acoustic(:)

    minus_int_pressure_acoustic(:) = minus_int_pressure_acoustic(:) + deltatover2 * minus_pressure_acoustic(:)

    ! update the pressure field (use a new array here) for coupling terms
    if (SIMULATION_TYPE == 3) then
      minus_int_int_pressure_acoustic_adj_coupling(:) = minus_int_int_pressure_acoustic(:) + &
                                           deltat * minus_int_pressure_acoustic(:) + &
                                           deltatsquareover2 * minus_pressure_acoustic(:)
    endif
  endif

  if (time_stepping_scheme == 2) then
    ! LDDRK scheme
    !! DK DK this should be vectorized
    minus_pressure_acoustic(:) = minus_pressure_acoustic(:) * rmass_inverse_acoustic(:)

    minus_int_pressure_acoustic_LDDRK(:) = ALPHA_LDDRK(i_stage) * minus_int_pressure_acoustic_LDDRK(:) + &
                                      deltat * minus_pressure_acoustic(:)
    minus_int_int_pressure_acoustic_LDDRK(:) = ALPHA_LDDRK(i_stage) * minus_int_int_pressure_acoustic_LDDRK(:) + &
                                  deltat * minus_int_pressure_acoustic(:)

    if (i_stage==1 .and. it == 1 .and. (.not. initialfield)) then
      !! DK DK this should be vectorized
      minus_int_pressure_acoustic_temp(:) = minus_int_pressure_acoustic_temp(:) + &
                                       BETA_LDDRK(i_stage) * minus_int_pressure_acoustic_LDDRK(:)
      minus_int_pressure_acoustic(:) = minus_int_pressure_acoustic_temp(:)
    else
      minus_int_pressure_acoustic(:) = minus_int_pressure_acoustic(:) + BETA_LDDRK(i_stage) * minus_int_pressure_acoustic_LDDRK(:)
    endif

    !! DK DK this should be vectorized
    minus_int_int_pressure_acoustic(:) = minus_int_int_pressure_acoustic(:) + &
        BETA_LDDRK(i_stage) * minus_int_int_pressure_acoustic_LDDRK(:)
  endif

  if (time_stepping_scheme == 3) then
    ! RK scheme
    !! DK DK this should be vectorized
    minus_pressure_acoustic(:) = minus_pressure_acoustic(:) * rmass_inverse_acoustic(:)

    minus_pressure_acoustic_rk(:,i_stage) = deltat * minus_pressure_acoustic(:)
    minus_int_pressure_acoustic_rk(:,i_stage) = deltat * minus_int_pressure_acoustic(:)

    if (i_stage == 1 .or. i_stage == 2 .or. i_stage == 3) then
      if (i_stage == 1) weight_rk = 0.5d0
      if (i_stage == 2) weight_rk = 0.5d0
      if (i_stage == 3) weight_rk = 1.0d0

      if (i_stage == 1) then
!! DK DK this should be vectorized
        minus_int_pressure_acoustic_init_rk(:) = minus_int_pressure_acoustic(:)
        minus_int_int_pressure_acoustic_init_rk(:) = minus_int_int_pressure_acoustic(:)
      endif
!! DK DK this should be vectorized
      minus_int_pressure_acoustic(:) = minus_int_pressure_acoustic_init_rk(:) + &
                                  weight_rk * minus_pressure_acoustic_rk(:,i_stage)
      minus_int_int_pressure_acoustic(:) = minus_int_int_pressure_acoustic_init_rk(:) + &
                                  weight_rk * minus_int_pressure_acoustic_rk(:,i_stage)
    else if (i_stage == 4) then
!! DK DK this should be vectorized
      minus_int_pressure_acoustic(:) = minus_int_pressure_acoustic_init_rk(:) + &
                                  1.0d0 / 6.0d0 * ( minus_pressure_acoustic_rk(:,1) + &
                                                    2.0d0 * minus_pressure_acoustic_rk(:,2) + &
                                                    2.0d0 * minus_pressure_acoustic_rk(:,3) + &
                                                    minus_pressure_acoustic_rk(:,4) )

!! DK DK this should be vectorized
      minus_int_int_pressure_acoustic(:) = minus_int_int_pressure_acoustic_init_rk(:) + &
                              1.0d0 / 6.0d0 * ( minus_int_pressure_acoustic_rk(:,1) + &
                                                2.0d0 * minus_int_pressure_acoustic_rk(:,2) + &
                                                2.0d0 * minus_int_pressure_acoustic_rk(:,3) + &
                                                minus_int_pressure_acoustic_rk(:,4) )
    endif
  endif

  ! free surface for an acoustic medium
  call enforce_acoustic_free_surface(minus_pressure_acoustic,minus_int_pressure_acoustic,minus_int_int_pressure_acoustic)

  end subroutine compute_forces_acoustic_main

!
!-------------------------------------------------------------------------------------
!

  subroutine compute_forces_acoustic_main_backward()

  use specfem_par

  implicit none

  ! local parameters
  integer :: it_temp,istage_temp

  ! checks
  if (SIMULATION_TYPE /= 3 ) return

  ! checks if anything to do in this slice
  if (.not. any_acoustic) return

  ! timing
  if (UNDO_ATTENUATION) then
    ! time increment
    ! example: NSTEP = 800, NT_DUMP_ATTENUATION = 500 -> 1. subset: it_temp = (2-1)*500 + 1 = 501,502,..,800
    !                                                 -> 2. subset: it_temp = (2-2)*500 + 1 = 1,2,..,500
    it_temp = (NSUBSET_ITERATIONS - iteration_on_subset)*NT_DUMP_ATTENUATION + it_of_this_subset
    ! time scheme
    istage_temp = i_stage
  else
    ! time increment
    it_temp = NSTEP - it + 1
    ! time scheme
    istage_temp = stage_time_scheme - i_stage + 1
  endif

  ! PML
  if (PML_BOUNDARY_CONDITIONS) then
    call rebuild_value_on_PML_interface_acoustic(it_temp)
  endif

  ! free surface for an acoustic medium
  call enforce_acoustic_free_surface(b_minus_pressure_acoustic,b_minus_int_pressure_acoustic,b_minus_int_int_pressure_acoustic)

  ! main solver for the acoustic elements
  if (UNDO_ATTENUATION) then
    call compute_forces_acoustic(b_minus_pressure_acoustic,b_minus_int_pressure_acoustic,b_minus_int_int_pressure_acoustic, &
                                 .false.,b_minus_int_int_pressure_acoustic_old)
  else
    call compute_forces_acoustic_backward(b_minus_pressure_acoustic,b_minus_int_int_pressure_acoustic)
  endif

  ! Stacey boundary conditions
  if (STACEY_BOUNDARY_CONDITIONS) then
    if (UNDO_ATTENUATION) then
      call compute_stacey_acoustic(b_minus_pressure_acoustic,b_minus_int_pressure_acoustic)
    else
      call compute_stacey_acoustic_backward(b_minus_pressure_acoustic)
    endif
  endif

  ! PML boundary conditions
  if (PML_BOUNDARY_CONDITIONS) then
    call pml_boundary_acoustic(b_minus_pressure_acoustic,b_minus_int_pressure_acoustic, &
                               b_minus_int_int_pressure_acoustic,b_minus_int_int_pressure_acoustic_old)
  endif

  ! PML
  if (PML_BOUNDARY_CONDITIONS) then
    call rebuild_value_on_PML_interface_acoustic(it_temp)
  endif

  ! add acoustic forcing at a rigid boundary
  if (ACOUSTIC_FORCING) then
    call add_acoustic_forcing_at_rigid_boundary(b_minus_pressure_acoustic)
  endif

  ! add coupling with the elastic side
  if (coupled_acoustic_elastic) then
    call compute_coupling_acoustic_el_backward(b_displ_elastic,b_minus_pressure_acoustic)
  endif

  ! add coupling with the poroelastic side
  if (coupled_acoustic_poro) then
    call compute_coupling_acoustic_po_backward()
  endif

  if (PML_BOUNDARY_CONDITIONS) then
    call rebuild_value_on_PML_interface_acoustic(it_temp)
  endif

  ! add force source
  if (.not. initialfield) then
    ! backward wavefield
    call compute_add_sources_acoustic(b_minus_pressure_acoustic,it_temp,istage_temp)
  endif

  ! assembling minus_pressure or b_minus_pressure for acoustic elements
#ifdef USE_MPI
  if (NPROC > 1 .and. ninterface_acoustic > 0) then
    call assemble_MPI_vector_ac(b_minus_pressure_acoustic)
  endif
#endif

  if (PML_BOUNDARY_CONDITIONS) then
    call rebuild_value_on_PML_interface_acoustic_accel(it_temp)
  endif

  ! multiply by the inverse of the mass matrix and update velocity
  if (time_stepping_scheme == 1) then
    !! DK DK this should be vectorized
    b_minus_pressure_acoustic(:) = b_minus_pressure_acoustic(:) * rmass_inverse_acoustic(:)

    b_minus_int_pressure_acoustic(:) = b_minus_int_pressure_acoustic(:) + b_deltatover2 * b_minus_pressure_acoustic(:)
  endif

  ! free surface for an acoustic medium
  call enforce_acoustic_free_surface(b_minus_pressure_acoustic,b_minus_int_pressure_acoustic,b_minus_int_int_pressure_acoustic)

  end subroutine compute_forces_acoustic_main_backward



