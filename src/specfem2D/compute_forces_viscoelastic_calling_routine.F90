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

  subroutine compute_forces_viscoelastic_main()

  use constants, only: SOURCE_IS_MOVING,USE_ENFORCE_FIELDS,ALPHA_LDDRK,BETA_LDDRK
  use specfem_par
  use specfem_par_noise

  implicit none

  ! local parameters
  integer :: i,iglob
  ! non-blocking MPI
  ! iphase: iphase = 1 is for computing outer elements (on MPI interface),
  !         iphase = 2 is for computing inner elements
  integer :: iphase

  ! checks if anything to do in this slice
  if ((.not. any_elastic) .and. (.not. SOURCE_IS_MOVING)) return

  ! enforces vanishing wavefields on axis
  if (AXISYM) then
    call enforce_zero_radial_displacements_on_the_axis()
  endif

  ! distinguishes two runs: for elements on MPI interfaces (outer), and elements within the partitions (inner)
  do iphase = 1,2

    ! main solver for the elastic elements
    ! visco-elastic term
    call compute_forces_viscoelastic(accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old,dux_dxl_old,duz_dzl_old, &
                                     dux_dzl_plus_duz_dxl_old,PML_BOUNDARY_CONDITIONS,e1,e11,e13,iphase)

    ! computes additional contributions to acceleration field
    if (iphase == 1) then

      ! Stacey boundary conditions
      if (STACEY_ABSORBING_CONDITIONS) then
        call compute_stacey_elastic(accel_elastic,veloc_elastic)
      endif

      ! PML boundary
      if (PML_BOUNDARY_CONDITIONS) then
        call pml_boundary_elastic(accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old)
      endif

      ! add coupling with the acoustic side
      if (coupled_acoustic_elastic) then
        call compute_coupling_viscoelastic_ac()
      endif

      ! add coupling with the poroelastic side
      if (coupled_elastic_poro) then
        call compute_coupling_viscoelastic_po()
      endif

      ! add force source
      if (.not. initialfield) then

        select case(NOISE_TOMOGRAPHY)
        case (0)
          ! earthquake/force source
          if (SIMULATION_TYPE == 1) then
            if (SOURCE_IS_MOVING) then
              call compute_add_sources_viscoelastic_moving_source(accel_elastic,it,i_stage)
            else
              call compute_add_sources_viscoelastic(accel_elastic,it,i_stage)
            endif
          endif

        case (1)
          ! noise source at master station
          call add_point_source_noise()

        case (2)
          ! inject generating wavefield for noise simulations
          call add_surface_movie_noise(accel_elastic)
        end select

        ! adjoint wavefield source
        if (SIMULATION_TYPE == 3) then
          ! adjoint sources
          call compute_add_sources_viscoelastic_adjoint()
        endif
      endif

    endif

    ! enforces vanishing wavefields on axis
    if (AXISYM) then
      call enforce_zero_radial_displacements_on_the_axis()
    endif

#ifdef USE_MPI
    ! LDDRK
    ! daniel: when is this needed? veloc_elastic at it == 1 and i_stage == 1 is zero for non-initialfield simulations.
    !         todo - please check...
    if (time_stepping_scheme == 2) then
      if (i_stage == 1 .and. it == 1 .and. iphase == 2 .and. (.not. initialfield)) then
        ! debug
        !print *,'debug veloc min/max = ',minval(veloc_elastic(:,:)),maxval(veloc_elastic(:,:))
        veloc_elastic_LDDRK_temp(:,:) = veloc_elastic(:,:)

        ! assembles velocity
        call assemble_MPI_vector_el_blocking(veloc_elastic)
      endif
    endif

    ! assemble all the contributions between slices using MPI
    if (NPROC > 1 .and. ninterface_elastic > 0) then
      ! collects all contributions on shared degrees of freedom
      !call assemble_MPI_vector_el_blocking(accel_elastic)

      ! sends out MPI interface data
      if (iphase == 1) then
        ! sends accel values to corresponding MPI interface neighbors
        call assemble_MPI_vector_el_s(accel_elastic)
      else
        ! waits for send/receive requests to be completed and assembles values
        call assemble_MPI_vector_el_w(accel_elastic)
      endif
    endif
#endif

  enddo ! iphase

  ! saves boundary condition for reconstruction
  if (PML_BOUNDARY_CONDITIONS) then
    if (nglob_interface > 0) then
      if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
        do i = 1, nglob_interface
          write(71) accel_elastic(1,point_interface(i)),accel_elastic(2,point_interface(i)), &
                    veloc_elastic(1,point_interface(i)),veloc_elastic(2,point_interface(i)), &
                    displ_elastic(1,point_interface(i)),displ_elastic(2,point_interface(i))
        enddo
      endif
    endif
  endif

  ! multiply by the inverse of the mass matrix and update velocity
  !! DK DK this should be vectorized
  if (USE_ENFORCE_FIELDS) then
    do iglob = 1,nglob_elastic
      if (.not. iglob_is_forced(iglob)) then
        accel_elastic(:,iglob) = accel_elastic(:,iglob) * rmass_inverse_elastic(:,iglob)
      endif
    enddo
  else
    accel_elastic(:,:) = accel_elastic(:,:) * rmass_inverse_elastic(:,:)
  endif

  ! time stepping
  select case (time_stepping_scheme)
  case (1)
    ! Newmark
    call update_veloc_elastic_Newmark()
  case (2)
    ! LDDRK
    call update_veloc_elastic_LDDRK()
  case (3)
    ! RK
    call update_veloc_elastic_RK()
  end select

  end subroutine compute_forces_viscoelastic_main

!
!-------------------------------------------------------------------------------------
!

  subroutine compute_forces_viscoelastic_main_backward()

  use constants, only: NOISE_SAVE_EVERYWHERE
  use specfem_par
  use specfem_par_noise

  implicit none

  ! local parameters
  integer :: it_temp,istage_temp
  ! non blocking MPI
  ! iphase: iphase = 1 is for computing outer elements (on MPI interface),
  !         iphase = 2 is for computing inner elements
  integer :: iphase

  ! checks if anything to do
  if (SIMULATION_TYPE /= 3 ) return

  ! checks if anything to do in this slice
  if (.not. any_elastic) return

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
    ! example: NSTEP = 800 -> 800,799,..,1
    it_temp = NSTEP - it + 1
    ! time scheme
    istage_temp = stage_time_scheme - i_stage + 1
  endif

  ! main solver for the elastic elements

  !ZN currently we do not support plane wave source in adjoint inversion
  if (PML_BOUNDARY_CONDITIONS) then
    call rebuild_value_on_PML_interface_viscoelastic(it_temp)
  endif

  ! distinguishes two runs: for elements on MPI interfaces (outer), and elements within the partitions (inner)
  do iphase = 1,2

    call compute_forces_viscoelastic(b_accel_elastic,b_veloc_elastic,b_displ_elastic,b_displ_elastic_old, &
                                     b_dux_dxl_old,b_duz_dzl_old, &
                                     b_dux_dzl_plus_duz_dxl_old,.false.,b_e1,b_e11,b_e13,iphase)

    ! computes additional contributions
    if (iphase == 1) then
      ! Stacey boundary conditions
      if (STACEY_ABSORBING_CONDITIONS) then
        if (UNDO_ATTENUATION_AND_OR_PML) then
          call compute_stacey_elastic(b_accel_elastic,b_veloc_elastic)
        else
          call compute_stacey_elastic_backward(b_accel_elastic)
        endif
      endif

      ! PML boundary
      if (PML_BOUNDARY_CONDITIONS) then
        call pml_boundary_elastic(b_accel_elastic,b_veloc_elastic,b_displ_elastic,b_displ_elastic_old)
      endif

      if (PML_BOUNDARY_CONDITIONS) then
        call rebuild_value_on_PML_interface_viscoelastic(it_temp)
      endif

      ! add coupling with the acoustic side
      if (coupled_acoustic_elastic) then
        call compute_coupling_viscoelastic_ac_backward()
      endif

      ! add coupling with the poroelastic side
      if (coupled_elastic_poro) then
        call compute_coupling_viscoelastic_po_backward()
      endif

      ! only on forward arrays so far implemented...
      !if (AXISYM) then
      !  call enforce_zero_radial_displacements_on_the_axis()
      !endif

      ! add force source
      if (.not. initialfield) then

        select case (NOISE_TOMOGRAPHY)
        case (0)
          ! earthquake/force source
          ! for backward wavefield
          call compute_add_sources_viscoelastic(b_accel_elastic,it_temp,istage_temp)

        case (3)
          ! noise simulation
          ! reconstruction/backward wavefield
          ! injects generating wavefield sources
          if (.not. NOISE_SAVE_EVERYWHERE) call add_surface_movie_noise(b_accel_elastic)
        end select
      endif ! if not using an initial field

    endif ! iphase

#ifdef USE_MPI
    ! assembling accel_elastic for elastic elements
    if (NPROC > 1 .and. ninterface_elastic > 0) then
      if (iphase == 1) then
        call assemble_MPI_vector_el_s(b_accel_elastic)
      else
        call assemble_MPI_vector_el_w(b_accel_elastic)
      endif
    endif
#endif

  enddo ! iphase

  if (PML_BOUNDARY_CONDITIONS) then
    call rebuild_value_on_PML_interface_viscoelastic_accel(it_temp)
  endif

  ! multiply by the inverse of the mass matrix and update velocity
  !! DK DK this should be vectorized
  b_accel_elastic(:,:) = b_accel_elastic(:,:) * rmass_inverse_elastic(:,:)

  ! time stepping
  select case (time_stepping_scheme)
  case (1)
    ! Newmark
    call update_veloc_elastic_Newmark_backward()
  case default
    call stop_the_code('Time stepping scheme not implemented yet in viscoelastic backward routine')
  end select

  end subroutine compute_forces_viscoelastic_main_backward

