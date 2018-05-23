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

  subroutine compute_forces_viscoacoustic_main()

  use constants, only: SOURCE_IS_MOVING,USE_ENFORCE_FIELDS,ALPHA_LDDRK,BETA_LDDRK,ZERO,USE_A_STRONG_FORMULATION_FOR_E1

  use specfem_par

  implicit none

  ! local parameters
  integer :: i,iglob,i_sls
  ! non-blocking MPI
  ! iphase: iphase = 1 is for computing outer elements (on MPI interface),
  !         iphase = 2 is for computing inner elements
  integer :: iphase

  ! checks if anything to do in this slice
  if ((.not. any_acoustic) .and. (.not. SOURCE_IS_MOVING)) return

  ! free surface for an acoustic medium
  call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic)

!! DK DK QUENTIN visco begin
  ! viscoacoustic attenuation for fluid media
  if (ATTENUATION_VISCOACOUSTIC .and. (.not. USE_A_STRONG_FORMULATION_FOR_E1) ) then

    ! If Newmark scheme we compute the memory variables with a Delta_t offset from the displacement
    if (time_stepping_scheme == 1) then

      do iphase = 1,2

        call compute_attenuation_acoustic_integration(potential_acoustic,ispec_is_acoustic,PML_BOUNDARY_CONDITIONS,iphase,dot_e1)

#ifdef USE_MPI
        ! assembling potential_dot_dot or b_potential_dot_dot for acoustic elements
        if (NPROC > 1 .and. ninterface_acoustic > 0) then
        ! loop over relaxation mechanisms
        do i_sls = 1,N_SLS
          if (iphase == 1) then
            call assemble_MPI_scalar_ac_s_e1(dot_e1(:,i_sls),dot_e1,0)
          else
            call assemble_MPI_scalar_ac_w_e1(dot_e1(:,i_sls),dot_e1,0)
          endif
        enddo
        endif
#endif

      enddo

      ! multiply by the inverse of the mass matrix
      !! DK DK this should be vectorized
      if (USE_ENFORCE_FIELDS) then
        do iglob = 1,nglob_acoustic
          if (.not. iglob_is_forced(iglob)) then
            do i_sls = 1,N_SLS
              dot_e1(iglob,i_sls) = dot_e1(iglob,i_sls) * rmass_inverse_e1(iglob,i_sls)
            enddo
          endif
        enddo
      else
        dot_e1(:,:) = rmass_inverse_e1(:,:)*dot_e1(:,:)
      endif

    endif ! of if test on time_stepping_scheme == 1

    call update_memory_var_acous_weak_form(dot_e1)

  endif

    ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
    do iphase = 1,2

    ! main solver for the acoustic elements
    call compute_forces_viscoacoustic(potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic, &
                                 PML_BOUNDARY_CONDITIONS,potential_acoustic_old,iphase,e1_acous_sf,sum_forces_old)

    ! PML boundary conditions enforces zero potentials on boundary
    if (PML_BOUNDARY_CONDITIONS) then
      call pml_boundary_acoustic(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                 potential_acoustic,potential_acoustic_old)
    endif

    ! computes additional contributions
    if (iphase == 1) then
      ! Stacey boundary conditions
      if (STACEY_ABSORBING_CONDITIONS) then
        call compute_stacey_acoustic(potential_dot_dot_acoustic,potential_dot_acoustic)
      endif

      ! add acoustic forcing at a rigid boundary
      if (ACOUSTIC_FORCING .and. (.not. USE_ENFORCE_FIELDS)) then
        call add_acoustic_forcing_at_rigid_boundary(potential_dot_dot_acoustic,dot_e1)
      endif

      ! applies to coupling in case of MPI partitioning:
      !   coupling interfaces might not be properly detected if one material domain is only in an another slice.
      !   in such a case, the common nodes would not be detected as belonging to a coupling interface.
      !   something to do in future...

      ! add coupling with the elastic side
      if (coupled_acoustic_elastic) then
        if (SIMULATION_TYPE == 1) then
          call compute_coupling_acoustic_el(displ_elastic,displ_elastic_old,potential_dot_dot_acoustic,dot_e1)
        endif

        ! coupling for adjoint wavefields
        if (SIMULATION_TYPE == 3) then
          ! note: handles adjoint runs coupling between adjoint potential and adjoint elastic wavefield
          !       adjoint definition: \partial_t^2 \bfs^\dagger = - \frac{1}{\rho} \bfnabla \phi^\dagger
          !
          ! coupling with adjoint wavefields
          call compute_coupling_acoustic_el(accel_elastic_adj_coupling,displ_elastic_old,potential_dot_dot_acoustic,dot_e1)
        endif
      endif

      ! add coupling with the poroelastic side
      if (coupled_acoustic_poro) then
        call compute_coupling_acoustic_po(dot_e1)
      endif

      ! add force source
      if (.not. initialfield) then
        if (SIMULATION_TYPE == 1) then
          if (SOURCE_IS_MOVING) then
            call compute_add_sources_acoustic_moving_source(potential_dot_dot_acoustic,it,i_stage)
          else
            call compute_add_sources_acoustic(potential_dot_dot_acoustic,it,i_stage)
          endif
        else if (SIMULATION_TYPE == 3) then
          ! adjoint sources
          call compute_add_sources_acoustic_adjoint()
        endif
      endif
    endif ! iphase == 1

#ifdef USE_MPI
    ! assembling potential_dot_dot or b_potential_dot_dot for acoustic elements
    if (NPROC > 1 .and. ninterface_acoustic > 0) then
      if (iphase == 1) then
        call assemble_MPI_scalar_ac_s_e1(potential_dot_dot_acoustic,dot_e1,N_SLS)
      else
        call assemble_MPI_scalar_ac_w_e1(potential_dot_dot_acoustic,dot_e1,N_SLS)
      endif
      if (time_stepping_scheme == 2) then
        ! LDDRK
        if (i_stage == 1 .and. it == 1 .and. iphase == 2 .and. (.not. initialfield)) then
          potential_dot_acoustic_temp(:) = potential_dot_acoustic(:)
          call assemble_MPI_scalar_ac_blocking(potential_dot_acoustic)
        endif
      endif
    endif
#endif

  enddo ! iphase

  ! PML saves interface values
  if (PML_BOUNDARY_CONDITIONS) then
    if (nglob_interface > 0) then
      if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
        do i = 1, nglob_interface
          write(72) potential_dot_dot_acoustic(point_interface(i)), potential_dot_acoustic(point_interface(i)), &
                    potential_acoustic(point_interface(i))
        enddo
      endif
    endif
  endif

  ! multiply by the inverse of the mass matrix
  !! DK DK this should be vectorized
  if (USE_ENFORCE_FIELDS) then
    do iglob = 1,nglob_acoustic
      if (.not. iglob_is_forced(iglob)) then
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) * rmass_inverse_acoustic(iglob)
        if (ATTENUATION_VISCOACOUSTIC .and. (.not. USE_A_STRONG_FORMULATION_FOR_E1) .and. time_stepping_scheme > 1) &
                dot_e1(iglob,:) = dot_e1(iglob,:) * rmass_inverse_e1(iglob,:)
      endif
    enddo
  else
    potential_dot_dot_acoustic(:) = potential_dot_dot_acoustic(:) * rmass_inverse_acoustic(:)
    if (ATTENUATION_VISCOACOUSTIC .and. (.not. USE_A_STRONG_FORMULATION_FOR_E1) .and. time_stepping_scheme > 1) &
        dot_e1(:,:) = rmass_inverse_e1(:,:)*dot_e1(:,:)
  endif

!! DK DK QUENTIN visco end

  ! update velocity
  select case (time_stepping_scheme)
  case (1)
    ! Newmark scheme
    call update_veloc_acoustic_Newmark()
  case (2)
    ! LDDRK scheme
    call update_veloc_acoustic_LDDRK()
  case (3)
    ! RK scheme
    call update_veloc_acoustic_RK()
  case default
    call stop_the_code('Invalid time stepping scheme for compute forces routine!')
  end select

  ! free surface for an acoustic medium
  call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic)

  end subroutine compute_forces_viscoacoustic_main

!
!-------------------------------------------------------------------------------------
!

  subroutine compute_forces_viscoacoustic_main_backward()

  use specfem_par

  implicit none

  ! local parameters
  integer :: it_temp,istage_temp
  ! non-blocking MPI
  ! iphase: iphase = 1 is for computing outer elements (on MPI interface),
  !         iphase = 2 is for computing inner elements
  integer :: iphase

  ! checks
  if (SIMULATION_TYPE /= 3 ) return

  ! checks if anything to do in this slice
  if (.not. any_acoustic) return

  ! timing
  if (UNDO_ATTENUATION_AND_OR_PML) then
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

  ! PML restores interface values
  if (PML_BOUNDARY_CONDITIONS) then
    call rebuild_value_on_PML_interface_acoustic(it_temp,b_potential_acoustic,b_potential_dot_acoustic)
  endif

  ! free surface for an acoustic medium
  call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic)

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase = 1,2

    ! main solver for the acoustic elements
    call compute_forces_viscoacoustic(b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic, &
                                     .false.,b_potential_acoustic_old,iphase,b_e1_acous_sf,b_sum_forces_old)

    ! PML boundary conditions
    if (PML_BOUNDARY_CONDITIONS) then
      ! enforces zero potentials on boundary
      call pml_boundary_acoustic(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                 b_potential_acoustic,b_potential_acoustic_old)

      ! restores potentials on interface
      call rebuild_value_on_PML_interface_acoustic(it_temp,b_potential_acoustic,b_potential_dot_acoustic)
    endif

    ! computes additional contributions
    if (iphase == 1) then
      ! Stacey boundary conditions
      if (STACEY_ABSORBING_CONDITIONS) then
        if (UNDO_ATTENUATION_AND_OR_PML) then
!! DK DK March 2018: following Quentin Brissaud's new variational implementation of viscoacousticity,
!! DK DK March 2018: I am not sure if the last argument, dot_e1, that I added here is right, or if
!! DK DK March 2018: a backward one should be created (b_dot_e1); far more likely it should be b_dot_e1 to create (not done yet)
          call compute_stacey_acoustic(b_potential_dot_dot_acoustic,b_potential_dot_acoustic)
        else
          call compute_stacey_acoustic_backward(b_potential_dot_dot_acoustic)
        endif
      endif

      ! add acoustic forcing at a rigid boundary
      if (ACOUSTIC_FORCING) then
!! DK DK March 2018: following Quentin Brissaud's new variational implementation of viscoacousticity,
!! DK DK March 2018: I am not sure if the last argument, dot_e1, that I added here is right, or if
!! DK DK March 2018: a backward one should be created (b_dot_e1); far more likely it should be b_dot_e1 to create (not done yet)
        call add_acoustic_forcing_at_rigid_boundary(b_potential_dot_dot_acoustic,dot_e1)
      endif

      ! add coupling with the elastic side
      if (coupled_acoustic_elastic) then
        call compute_coupling_acoustic_el_backward(b_displ_elastic,b_potential_dot_dot_acoustic)
      endif

      ! add coupling with the poroelastic side
      if (coupled_acoustic_poro) then
        call compute_coupling_acoustic_po_backward()
      endif

      ! PML restores interface values
      if (PML_BOUNDARY_CONDITIONS) then
        call rebuild_value_on_PML_interface_acoustic(it_temp,b_potential_acoustic,b_potential_dot_acoustic)
      endif

      ! add force source
      if (.not. initialfield) then
        ! backward wavefield
        call compute_add_sources_acoustic(b_potential_dot_dot_acoustic,it_temp,istage_temp)
      endif

    endif ! iphase == 1

#ifdef USE_MPI
    ! assembling potential_dot_dot or b_potential_dot_dot for acoustic elements
    if (NPROC > 1 .and. ninterface_acoustic > 0) then
      if (iphase == 1) then
        call assemble_MPI_scalar_ac_s(b_potential_dot_dot_acoustic)
      else
        call assemble_MPI_scalar_ac_w(b_potential_dot_dot_acoustic)
      endif
    endif
#endif

  enddo ! iphase

  ! PML restores interface values
  if (PML_BOUNDARY_CONDITIONS) then
    call rebuild_value_on_PML_interface_acoustic_accel(it_temp,b_potential_dot_dot_acoustic)
  endif

  ! multiply by the inverse of the mass matrix
  b_potential_dot_dot_acoustic(:) = b_potential_dot_dot_acoustic(:) * rmass_inverse_acoustic(:)

  ! update velocity
  select case (time_stepping_scheme)
  case (1)
    ! Newmark
    call update_veloc_acoustic_Newmark_backward()
  case default
    call stop_the_code('Sorry, time stepping scheme not implemented yet for backward computations')
  end select

  ! free surface for an acoustic medium
  call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic)

  end subroutine compute_forces_viscoacoustic_main_backward



