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

  subroutine compute_forces_viscoelastic_main()

  use specfem_par
  use specfem_par_noise

  implicit none

  ! local parameters
  integer :: i
  ! for rk44
  double precision :: weight_rk

  ! main solver for the elastic elements

  ! checks if anything to do in this slice
  if (.not. any_elastic) return

  ! enforces vanishing wavefields on axis
  if (AXISYM) then
    call enforce_zero_radial_displacements_on_the_axis()
  endif

  ! visco-elastic term
  call compute_forces_viscoelastic(accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old, &
                                   PML_BOUNDARY_CONDITIONS,e1,e11,e13)

  ! Stacey boundary conditions
  if (STACEY_BOUNDARY_CONDITIONS) then
    call compute_stacey_elastic(accel_elastic,veloc_elastic,displ_elastic)
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

  ! enforces vanishing wavefields on axis
  if (AXISYM) then
    call enforce_zero_radial_displacements_on_the_axis()
  endif

  ! add force source
  if (.not. initialfield) then
    if (SIMULATION_TYPE == 1) then
      call compute_add_sources_viscoelastic(accel_elastic,it,i_stage)
    else if (SIMULATION_TYPE == 3) then
      ! adjoint sources
      call compute_add_sources_viscoelastic_adjoint()
    endif

    !<NOISE_TOMOGRAPHY
    ! inject wavefield sources for noise simulations
    if (NOISE_TOMOGRAPHY == 1) then
      call add_point_source_noise()
    else if (NOISE_TOMOGRAPHY == 2) then
      call add_surface_movie_noise(accel_elastic)
    endif
    !>NOISE_TOMOGRAPHY
  endif

  ! enforces vanishing wavefields on axis
  if (AXISYM) then
    call enforce_zero_radial_displacements_on_the_axis()
  endif

  ! assembling accel_elastic for elastic elements
#ifdef USE_MPI
  if (NPROC > 1 .and. ninterface_elastic > 0) then
    if (time_stepping_scheme == 2) then
      if (i_stage==1 .and. it == 1 .and. (.not. initialfield)) then
        veloc_elastic_LDDRK_temp = veloc_elastic
        call assemble_MPI_vector_el(veloc_elastic)
      endif
    endif

    call assemble_MPI_vector_el(accel_elastic)
  endif
#endif

  if (PML_BOUNDARY_CONDITIONS) then
    if (nglob_interface > 0) then
      if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
        do i = 1, nglob_interface
          write(71) accel_elastic(1,point_interface(i)),accel_elastic(2,point_interface(i)),&
                    accel_elastic(3,point_interface(i)),&
                    veloc_elastic(1,point_interface(i)),veloc_elastic(2,point_interface(i)),&
                    veloc_elastic(3,point_interface(i)),&
                    displ_elastic(1,point_interface(i)),displ_elastic(2,point_interface(i)),&
                    displ_elastic(3,point_interface(i))
        enddo
      endif
    endif
  endif

  ! multiply by the inverse of the mass matrix and update velocity
  !! DK DK this should be vectorized
  accel_elastic(1,:) = accel_elastic(1,:) * rmass_inverse_elastic_one(:)
  accel_elastic(2,:) = accel_elastic(2,:) * rmass_inverse_elastic_one(:)
  accel_elastic(3,:) = accel_elastic(3,:) * rmass_inverse_elastic_three(:)

  if (time_stepping_scheme == 1) then
    !! DK DK this should be vectorized
    veloc_elastic(:,:) = veloc_elastic(:,:) + deltatover2 * accel_elastic(:,:)
  endif

  if (time_stepping_scheme == 2) then
    !! DK DK this should be vectorized
    veloc_elastic_LDDRK(:,:) = ALPHA_LDDRK(i_stage) * veloc_elastic_LDDRK(:,:) + deltat * accel_elastic(:,:)
    displ_elastic_LDDRK(:,:) = ALPHA_LDDRK(i_stage) * displ_elastic_LDDRK(:,:) + deltat * veloc_elastic(:,:)
    if (i_stage==1 .and. it == 1 .and. (.not. initialfield)) then
      veloc_elastic_LDDRK_temp(:,:) = veloc_elastic_LDDRK_temp(:,:) + BETA_LDDRK(i_stage) * veloc_elastic_LDDRK(:,:)
      veloc_elastic(:,:) = veloc_elastic_LDDRK_temp(:,:)
    else
      veloc_elastic(:,:) = veloc_elastic(:,:) + BETA_LDDRK(i_stage) * veloc_elastic_LDDRK(:,:)
    endif
    displ_elastic(:,:) = displ_elastic(:,:) + BETA_LDDRK(i_stage) * displ_elastic_LDDRK(:,:)
  endif

  if (time_stepping_scheme == 3) then
    !! DK DK this should be vectorized
    accel_elastic_rk(1,:,i_stage) = deltat * accel_elastic(1,:)
    accel_elastic_rk(2,:,i_stage) = deltat * accel_elastic(2,:)
    accel_elastic_rk(3,:,i_stage) = deltat * accel_elastic(3,:)

    veloc_elastic_rk(1,:,i_stage) = deltat * veloc_elastic(1,:)
    veloc_elastic_rk(2,:,i_stage) = deltat * veloc_elastic(2,:)
    veloc_elastic_rk(3,:,i_stage) = deltat * veloc_elastic(3,:)

    if (i_stage == 1 .or. i_stage == 2 .or. i_stage == 3) then

      if (i_stage == 1) weight_rk = 0.5d0
      if (i_stage == 2) weight_rk = 0.5d0
      if (i_stage == 3) weight_rk = 1.0d0

      if (i_stage == 1) then
        !! DK DK this should be vectorized
        veloc_elastic_initial_rk(1,:) = veloc_elastic(1,:)
        veloc_elastic_initial_rk(2,:) = veloc_elastic(2,:)
        veloc_elastic_initial_rk(3,:) = veloc_elastic(3,:)

        displ_elastic_initial_rk(1,:) = displ_elastic(1,:)
        displ_elastic_initial_rk(2,:) = displ_elastic(2,:)
        displ_elastic_initial_rk(3,:) = displ_elastic(3,:)
      endif

      !! DK DK this should be vectorized
      veloc_elastic(1,:) = veloc_elastic_initial_rk(1,:) + weight_rk * accel_elastic_rk(1,:,i_stage)
      veloc_elastic(2,:) = veloc_elastic_initial_rk(2,:) + weight_rk * accel_elastic_rk(2,:,i_stage)
      veloc_elastic(3,:) = veloc_elastic_initial_rk(3,:) + weight_rk * accel_elastic_rk(3,:,i_stage)

      displ_elastic(1,:) = displ_elastic_initial_rk(1,:) + weight_rk * veloc_elastic_rk(1,:,i_stage)
      displ_elastic(2,:) = displ_elastic_initial_rk(2,:) + weight_rk * veloc_elastic_rk(2,:,i_stage)
      displ_elastic(3,:) = displ_elastic_initial_rk(3,:) + weight_rk * veloc_elastic_rk(3,:,i_stage)

    else if (i_stage == 4) then
      !! DK DK this should be vectorized
      veloc_elastic(1,:) = veloc_elastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
                           ( accel_elastic_rk(1,:,1) + 2.0d0 * accel_elastic_rk(1,:,2) + &
                             2.0d0 * accel_elastic_rk(1,:,3) + accel_elastic_rk(1,:,4) )
      veloc_elastic(2,:) = veloc_elastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
                           ( accel_elastic_rk(2,:,1) + 2.0d0 * accel_elastic_rk(2,:,2) + &
                             2.0d0 * accel_elastic_rk(2,:,3) + accel_elastic_rk(2,:,4) )
      veloc_elastic(3,:) = veloc_elastic_initial_rk(3,:) + 1.0d0 / 6.0d0 * &
                           ( accel_elastic_rk(3,:,1) + 2.0d0 * accel_elastic_rk(3,:,2) + &
                             2.0d0 * accel_elastic_rk(3,:,3) + accel_elastic_rk(3,:,4) )

      displ_elastic(1,:) = displ_elastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
                           ( veloc_elastic_rk(1,:,1) + 2.0d0 * veloc_elastic_rk(1,:,2) + &
                             2.0d0 * veloc_elastic_rk(1,:,3) + veloc_elastic_rk(1,:,4) )
      displ_elastic(2,:) = displ_elastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
                           ( veloc_elastic_rk(2,:,1) + 2.0d0 * veloc_elastic_rk(2,:,2) + &
                             2.0d0 * veloc_elastic_rk(2,:,3) + veloc_elastic_rk(2,:,4))
      displ_elastic(3,:) = displ_elastic_initial_rk(3,:) + 1.0d0 / 6.0d0 * &
                           ( veloc_elastic_rk(3,:,1) + 2.0d0 * veloc_elastic_rk(3,:,2) + &
                             2.0d0 * veloc_elastic_rk(3,:,3) + veloc_elastic_rk(3,:,4))
    endif
  endif

  end subroutine compute_forces_viscoelastic_main

!
!-------------------------------------------------------------------------------------
!

  subroutine compute_forces_viscoelastic_main_backward()

  use specfem_par
  use specfem_par_noise

  implicit none

  ! local parameters
  integer :: it_temp,istage_temp

  ! checks if anything to do
  if (SIMULATION_TYPE /= 3 ) return

  ! checks if anything to do in this slice
  if (.not. any_elastic) return

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

  if (UNDO_ATTENUATION) then
    call compute_forces_viscoelastic(b_accel_elastic,b_veloc_elastic,b_displ_elastic,b_displ_elastic_old, &
                                     .false.,b_e1,b_e11,b_e13)
  else
    ! todo: maybe should be b_e1,b_e11,.. here, please check...
    call compute_forces_viscoelastic_backward(b_accel_elastic,b_displ_elastic,b_displ_elastic_old, &
                                              e1,e11,e13)
  endif

  ! Stacey boundary conditions
  if (STACEY_BOUNDARY_CONDITIONS) then
    if (UNDO_ATTENUATION) then
      call compute_stacey_elastic(b_accel_elastic,b_veloc_elastic,b_displ_elastic)
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
    ! backward wavefield
    call compute_add_sources_viscoelastic(b_accel_elastic,it_temp,istage_temp)

    !<NOISE_TOMOGRAPHY
    ! inject wavefield sources for noise simulations
    if (NOISE_TOMOGRAPHY == 3) then
      if (.not. save_everywhere) then
        call add_surface_movie_noise(b_accel_elastic)
      endif
    endif
    !>NOISE_TOMOGRAPHY
  endif ! if not using an initial field

  ! only on forward arrays so far implemented...
  !if (AXISYM) then
  !  call enforce_zero_radial_displacements_on_the_axis()
  !endif

  ! assembling accel_elastic for elastic elements
#ifdef USE_MPI
  if (NPROC > 1 .and. ninterface_elastic > 0) then
    call assemble_MPI_vector_el(b_accel_elastic)
  endif
#endif

  if (PML_BOUNDARY_CONDITIONS) then
    call rebuild_value_on_PML_interface_viscoelastic_accel(it_temp)
  endif

  ! multiply by the inverse of the mass matrix and update velocity
  !! DK DK this should be vectorized
  b_accel_elastic(1,:) = b_accel_elastic(1,:) * rmass_inverse_elastic_one(:)
  b_accel_elastic(2,:) = b_accel_elastic(2,:) * rmass_inverse_elastic_one(:)
  b_accel_elastic(3,:) = b_accel_elastic(3,:) * rmass_inverse_elastic_three(:)

  b_veloc_elastic = b_veloc_elastic + b_deltatover2*b_accel_elastic

  end subroutine compute_forces_viscoelastic_main_backward

