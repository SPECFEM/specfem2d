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

  subroutine iterate_time_undoatt()

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par
  use specfem_par_noise,only: NOISE_TOMOGRAPHY

  implicit none

#ifdef USE_MPI
  include "precision.h"
#endif

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_potential_acoustic_buffer,b_potential_dot_dot_acoustic_buffer
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_displ_elastic_buffer,b_accel_elastic_buffer

  integer :: i,j,iteration_on_subset,it_of_this_subset,it_backward
  integer :: it_temp,seismo_current_temp
  integer :: ier
  integer, parameter :: it_begin = 1

  ! time
  character(len=8) :: datein
  character(len=10) :: timein
  character(len=5) :: zone
  integer, dimension(8) :: time_values
  integer :: year,mon,day,hr,minutes,timestamp


  if (myrank == 0) write(IMAIN,400)
!
!----          s t a r t   t i m e   i t e r a t i o n s
!

  ! count elapsed wall-clock time
  call date_and_time(datein,timein,zone,time_values)
  ! time_values(1): year
  ! time_values(2): month of the year
  ! time_values(3): day of the month
  ! time_values(5): hour of the day
  ! time_values(6): minutes of the hour
  ! time_values(7): seconds of the minute
  ! time_values(8): milliseconds of the second

  ! get timestamp in minutes of current date and time
  year = time_values(1)
  mon = time_values(2)
  day = time_values(3)
  hr = time_values(5)
  minutes = time_values(6)
  call convtime(timestamp,year,mon,day,hr,minutes)

  ! convert to seconds instead of minutes, to be more precise for 2D runs, which can be fast
  timestamp_seconds_start = timestamp*60.d0 + time_values(7) + time_values(8)/1000.d0

! *********************************************************
! ************* MAIN LOOP OVER THE TIME STEPS *************
! *********************************************************

  ! safety checks
  if (GPU_MODE ) call exit_MPI(myrank,'for undo_attenuation, GPU_MODE is not supported')
  if (time_stepping_scheme /= 1 ) call exit_MPI(myrank,'for undo_attenuation, only Newmark scheme has implemented ')
  if (any_gravitoacoustic ) call exit_MPI(myrank,'undo_attenuation has not implemented for gravitoacoustic yet')
  if (any_poroelastic ) call exit_MPI(myrank,'undo_attenuation has not implemented for poroelastic simulation yet')
  if (NOISE_TOMOGRAPHY /= 0 ) call exit_MPI(myrank,'for undo_attenuation, NOISE_TOMOGRAPHY is not supported')
  if (AXISYM ) call exit_MPI(myrank,'Just axisymmetric FORWARD simulations are possible so far')

  if (SIMULATION_TYPE == 3) then
    if (any_acoustic) then
      allocate(b_potential_acoustic_buffer(nglob,NT_DUMP_ATTENUATION),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'error allocating b_potential_acoustic')
      allocate(b_potential_dot_dot_acoustic_buffer(nglob,NT_DUMP_ATTENUATION),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'error allocating b_potential_dot_dot_acoustic')
    endif

    if (any_elastic) then
      allocate(b_displ_elastic_buffer(3,nglob,NT_DUMP_ATTENUATION),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'error allocating b_displ_elastic')
      allocate(b_accel_elastic_buffer(3,nglob,NT_DUMP_ATTENUATION),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'error allocating b_accel_elastic')

      if (ATTENUATION_VISCOELASTIC_SOLID) then
        allocate(b_e1(NGLLX,NGLLZ,nspec_allocate,N_SLS))
        allocate(b_e11(NGLLX,NGLLZ,nspec_allocate,N_SLS))
        allocate(b_e13(NGLLX,NGLLZ,nspec_allocate,N_SLS))

        b_e1(:,:,:,:) = 0._CUSTOM_REAL
        b_e11(:,:,:,:) = 0._CUSTOM_REAL
        b_e13(:,:,:,:) = 0._CUSTOM_REAL
      endif
    endif
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! synchronize all processes to make sure everybody is ready to start time loop
  call synchronize_all()

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Starting time iteration loop in undoing attenuation...'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! initialize variables for writing seismograms
  seismo_offset = it_begin-1
  seismo_current = 0

  ! *********************************************************
  ! ************* MAIN LOOP OVER THE TIME STEPS *************
  ! *********************************************************

  it = 0
  do iteration_on_subset = 1, NSTEP / NT_DUMP_ATTENUATION

    ! wavefield storage
    if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
      ! saves forward wavefields
      call save_forward_arrays_undoatt(iteration_on_subset)

    else if (SIMULATION_TYPE == 3) then
      ! reads in last stored forward wavefield
      call read_forward_arrays_undoatt(iteration_on_subset)

      ! intermediate storage of it and seismo_current positions
      it_temp = it
      seismo_current_temp = seismo_current
    endif

    ! time loop within this iteration subset
    select case( SIMULATION_TYPE )
    case( 1 )
      ! forward and adjoint simulations
      ! subset loop
      do it_of_this_subset = 1, NT_DUMP_ATTENUATION

        it = it + 1

        do i_stage = 1, stage_time_scheme ! is equal to 1 if Newmark because only one stage then

          ! *********************************************************
          ! ************* update_displacement_newmark
          ! *********************************************************
          if (any_acoustic) then
            if (nelem_acoustic_surface > 0) then
              call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                                 potential_acoustic)
            endif
            call update_displacement_newmark_acoustic(deltat,deltatover2,deltatsquareover2,&
                                                      potential_dot_dot_acoustic,potential_dot_acoustic,&
                                                      potential_acoustic,potential_acoustic_old, &
                                                      PML_BOUNDARY_CONDITIONS)
          endif

          if (any_elastic) then
            call update_displacement_newmark_elastic(deltat,deltatover2,deltatsquareover2,&
                                                     accel_elastic,veloc_elastic,&
                                                     displ_elastic,displ_elastic_old,&
                                                     PML_BOUNDARY_CONDITIONS)
          endif
! *********************************************************
! ************* main solver for the acoustic elements
! *********************************************************
          if (any_acoustic) then
            call compute_forces_acoustic(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                         potential_acoustic,potential_acoustic_old,PML_BOUNDARY_CONDITIONS)
            ! *********************************************************
            ! ************* add acoustic forcing at a rigid boundary
            ! *********************************************************
            if (ACOUSTIC_FORCING) then
              call add_acoustic_forcing_at_rigid_boundary(potential_dot_dot_acoustic)
            endif
          endif
          ! *********************************************************
          ! ************* add coupling with the elastic side
          ! *********************************************************
          if (coupled_acoustic_elastic) then
            call compute_coupling_acoustic_el(displ_elastic,displ_elastic_old,potential_dot_dot_acoustic)
          endif
          ! ************************************************************************************
          ! ************************************ add force source
          ! ************************************************************************************
          if (any_acoustic) then
            if (.not. initialfield) then
              call compute_add_sources_acoustic(potential_dot_dot_acoustic,it,i_stage)
            endif
          endif
          ! ************************************************************************************
          ! ********** assembling potential_dot_dot or b_potential_dot_dot for acoustic elements
          ! ************************************************************************************
#ifdef USE_MPI
          if (NPROC > 1 .and. any_acoustic .and. ninterface_acoustic > 0) then
            call assemble_MPI_vector_ac(potential_dot_dot_acoustic)
          endif
#endif
          ! ************************************************************************************
          ! ********** save or read value on PML interface: potential_*_acoustic
          ! ************************************************************************************
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
          ! ************************************************************************************
          ! ************* multiply by the inverse of the mass matrix and update velocity
          ! ************************************************************************************
          if (any_acoustic) then
            ! free surface for an acoustic medium
            if (nelem_acoustic_surface > 0) then
              call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                                 potential_acoustic)
            endif

            !! DK DK this should be vectorized
            potential_dot_dot_acoustic = potential_dot_dot_acoustic * rmass_inverse_acoustic
            potential_dot_acoustic = potential_dot_acoustic + deltatover2*potential_dot_dot_acoustic
          endif
! *********************************************************
! ************* main solver for the elastic elements
! *********************************************************
          if (any_elastic) then
            call compute_forces_viscoelastic(accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old, &
                                             PML_BOUNDARY_CONDITIONS,e1,e11,e13)
          endif
          ! *********************************************************
          ! ************* add coupling with the acoustic side
          ! *********************************************************
          if (coupled_acoustic_elastic ) call compute_coupling_viscoelastic_ac()

          ! ************************************************************************************
          ! ************* add force source
          ! ************************************************************************************
          if (any_elastic) then
            if (.not. initialfield) then
              if (SIMULATION_TYPE == 1) then
                call compute_add_sources_viscoelastic(accel_elastic,it,i_stage)
              endif
            endif
          endif
          ! ************************************************************************************
          ! ************* assembling accel_elastic for elastic elements
          ! ************************************************************************************
#ifdef USE_MPI
          if (NPROC > 1 .and. any_elastic .and. ninterface_elastic > 0) then
            call assemble_MPI_vector_el(accel_elastic)
          endif
#endif
          ! ************************************************************************************
          ! ********** save or read value on PML interface: *_elastic
          ! ************************************************************************************
          if (PML_BOUNDARY_CONDITIONS) then
            if (any_elastic .and. nglob_interface > 0) then
              if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
                do i = 1, nglob_interface
                  write(71)accel_elastic(1,point_interface(i)),accel_elastic(2,point_interface(i)),&
                           accel_elastic(3,point_interface(i)),&
                           veloc_elastic(1,point_interface(i)),veloc_elastic(2,point_interface(i)),&
                           veloc_elastic(3,point_interface(i)),&
                           displ_elastic(1,point_interface(i)),displ_elastic(2,point_interface(i)),&
                           displ_elastic(3,point_interface(i))
                enddo
              endif
            endif
          endif
          ! ************************************************************************************
          ! ************* multiply by the inverse of the mass matrix and update velocity
          ! ************************************************************************************
          if (any_elastic) then
            !! DK DK this should be vectorized
            accel_elastic(1,:) = accel_elastic(1,:) * rmass_inverse_elastic_one
            accel_elastic(2,:) = accel_elastic(2,:) * rmass_inverse_elastic_one
            accel_elastic(3,:) = accel_elastic(3,:) * rmass_inverse_elastic_three

            veloc_elastic = veloc_elastic + deltatover2 * accel_elastic
          endif

        enddo ! i_stage

        ! computes kinetic and potential energy
        if (output_energy) then
          call it_compute_and_output_energy()
        endif

        ! display time step and max of norm of displacement
        if (mod(it,NSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then
          call check_stability()
        endif

        ! loop on all the receivers to compute and store the seismograms
        call write_seismograms()

        ! display results at given time steps
        call write_movie_output()

      enddo ! subset loop

    case( 3 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! kernel simulations

      ! reconstructs forward wavefield based on last stored wavefield data
      !
      ! note: we step forward in time here, starting from last snapshot.
      !       the newly computed, reconstructed forward wavefields (b_displ_..) get stored in buffers.

      ! subset loop
!****************************************************************************************************backward
!****************************************************************************************************
      do it_of_this_subset = 1, NT_DUMP_ATTENUATION

        it = it + 1
        do i_stage = 1, stage_time_scheme
          !****************************************************************************************************backward_inner_loop
          !****************************************************************************************************
          ! *********************************************************
          ! ************* update_displacement_newmark
          ! *********************************************************
          if (any_acoustic) then
            call update_displacement_newmark_acoustic(deltat,deltatover2,deltatsquareover2,&
                                                      b_potential_dot_dot_acoustic,b_potential_dot_acoustic,&
                                                      b_potential_acoustic,b_potential_acoustic_old, &
                                                      .false.)
          endif

          if (any_elastic) then
            call update_displacement_newmark_elastic(deltat,deltatover2,deltatsquareover2,&
                                                     b_accel_elastic,b_veloc_elastic,&
                                                     b_displ_elastic,b_displ_elastic_old,&
                                                     .false.)
          endif

! *********************************************************
! ************* main solver for the acoustic elements
! *********************************************************
          if (any_acoustic) then
            if (nelem_acoustic_surface > 0) then
              call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                                 b_potential_acoustic)
            endif

            if (PML_BOUNDARY_CONDITIONS) then
              it_backward = NSTEP-(iteration_on_subset*NT_DUMP_ATTENUATION-it_of_this_subset+1)
              call rebuild_value_on_PML_interface_acoustic(it_backward)  !ZN need to modify the it value inside this code
            endif
            call compute_forces_acoustic(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                         b_potential_acoustic,b_potential_acoustic_old,.false.)
            if (PML_BOUNDARY_CONDITIONS) then
              it_backward = NSTEP-(iteration_on_subset*NT_DUMP_ATTENUATION-it_of_this_subset+1)
              call rebuild_value_on_PML_interface_acoustic(it_backward)  !ZN need to modify the it value inside this code
            endif
            ! *********************************************************
            ! ************* add acoustic forcing at a rigid boundary
            ! *********************************************************
            if (ACOUSTIC_FORCING) then
              call add_acoustic_forcing_at_rigid_boundary(b_potential_dot_dot_acoustic)
            endif
          endif
          ! *********************************************************
          ! ************* add coupling with the elastic side
          ! *********************************************************
          if (coupled_acoustic_elastic) then
            call compute_coupling_acoustic_el_backward(b_displ_elastic,b_potential_dot_dot_acoustic)
          endif

          if (PML_BOUNDARY_CONDITIONS) then
            it_backward = NSTEP-(iteration_on_subset*NT_DUMP_ATTENUATION-it_of_this_subset+1)
            call rebuild_value_on_PML_interface_acoustic(it_backward)  !ZN need to modify the it value inside this code
          endif

          ! ************************************************************************************
          ! ************************************ add force source
          ! ************************************************************************************
          if (any_acoustic) then
            if (.not. initialfield) then  !ZN maybe need to change it
              it_backward = NSTEP-(iteration_on_subset*NT_DUMP_ATTENUATION-it_of_this_subset+1)
              call compute_add_sources_acoustic(b_potential_dot_dot_acoustic,it_backward,i_stage)
            endif
          endif
          ! ************************************************************************************
          ! ********** assembling b_potential_dot_dot for acoustic elements
          ! ************************************************************************************
#ifdef USE_MPI
          if (NPROC > 1 .and. any_acoustic .and. ninterface_acoustic > 0) then
            call assemble_MPI_vector_ac(b_potential_dot_dot_acoustic)
          endif
#endif
          ! ************************************************************************************
          ! ********** rebuild b_potential_dot_dot_acoustic on PML interface
          ! ************************************************************************************
          if (PML_BOUNDARY_CONDITIONS) then
            if (any_acoustic .and. nglob_interface > 0) then !ZN may need to change it
                do i = 1, nglob_interface
                  it_backward = NSTEP-(iteration_on_subset*NT_DUMP_ATTENUATION-it_of_this_subset+1)
                  b_potential_dot_dot_acoustic(point_interface(i)) = pml_interface_history_potential_dot_dot(i,it_backward)
                enddo
            endif
          endif
          ! ************************************************************************************
          ! ************* multiply by the inverse of the mass matrix and update velocity
          ! ************************************************************************************
          if (any_acoustic) then
            ! free surface for an acoustic medium
            if (nelem_acoustic_surface > 0) then
              if (SIMULATION_TYPE == 3) then
                call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                                   b_potential_acoustic)
              endif
            endif

            if (time_stepping_scheme == 1) then
              !! DK DK this should be vectorized
              b_potential_dot_dot_acoustic = b_potential_dot_dot_acoustic * rmass_inverse_acoustic
              b_potential_dot_acoustic = b_potential_dot_acoustic + deltatover2*b_potential_dot_dot_acoustic
            endif
          endif ! any_acoustic
! *********************************************************
! ************* main solver for the elastic elements
! *********************************************************
          if (any_elastic) then
            if (PML_BOUNDARY_CONDITIONS) then
              it_backward = NSTEP-(iteration_on_subset*NT_DUMP_ATTENUATION-it_of_this_subset+1)
              call rebuild_value_on_PML_interface_viscoelastic(it_backward)
            endif

            call compute_forces_viscoelastic(b_accel_elastic,b_veloc_elastic,b_displ_elastic,b_displ_elastic_old, &
                                             .false.,b_e1,b_e11,b_e13)

            if (PML_BOUNDARY_CONDITIONS) then
              it_backward = NSTEP-(iteration_on_subset*NT_DUMP_ATTENUATION-it_of_this_subset+1)
              call rebuild_value_on_PML_interface_viscoelastic(it_backward)
            endif
          endif ! any_elastic
          ! *********************************************************
          ! ************* add coupling with the acoustic side
          ! *********************************************************
          if (coupled_acoustic_elastic ) call compute_coupling_viscoelastic_ac_backward()

          if (PML_BOUNDARY_CONDITIONS) then
            it_backward = NSTEP-(iteration_on_subset*NT_DUMP_ATTENUATION-it_of_this_subset+1)
            call rebuild_value_on_PML_interface_viscoelastic(it_backward)
          endif

          ! ************************************************************************************
          ! ************* add force source
          ! ************************************************************************************
          if (any_elastic) then
            if (.not. initialfield) then
              it_backward = NSTEP-(iteration_on_subset*NT_DUMP_ATTENUATION-it_of_this_subset+1)
              call compute_add_sources_viscoelastic(b_accel_elastic,it_backward,i_stage)
            endif ! if not using an initial field
          endif !if (any_elastic)
          ! ************************************************************************************
          ! ************* assembling accel_elastic for elastic elements
          ! ************************************************************************************
#ifdef USE_MPI
          if (NPROC > 1 .and. any_elastic .and. ninterface_elastic > 0) then
            call assemble_MPI_vector_el(b_accel_elastic)
          endif
#endif
          ! ************************************************************************************
          ! ********** rebuild b_accel_elastic on PML interface
          ! ************************************************************************************
          if (PML_BOUNDARY_CONDITIONS) then
            if (any_elastic .and. nglob_interface > 0) then  !the it may have to change
              if (SIMULATION_TYPE == 3) then
                do i = 1, nglob_interface
                  it_backward = NSTEP-(iteration_on_subset*NT_DUMP_ATTENUATION-it_of_this_subset+1)
                  b_accel_elastic(1,point_interface(i)) = pml_interface_history_accel(1,i,it_backward)
                  b_accel_elastic(2,point_interface(i)) = pml_interface_history_accel(2,i,it_backward)
                  b_accel_elastic(3,point_interface(i)) = pml_interface_history_accel(3,i,it_backward)
                enddo
              endif
            endif
          endif

          ! ************************************************************************************
          ! ************* multiply by the inverse of the mass matrix and update velocity
          ! ************************************************************************************
          if (any_elastic) then
            !! DK DK this should be vectorized
            b_accel_elastic(1,:) = b_accel_elastic(1,:) * rmass_inverse_elastic_one(:)
            b_accel_elastic(2,:) = b_accel_elastic(2,:) * rmass_inverse_elastic_one(:)
            b_accel_elastic(3,:) = b_accel_elastic(3,:) * rmass_inverse_elastic_three(:)

            b_veloc_elastic = b_veloc_elastic + deltatover2*b_accel_elastic
          endif

          !****************************************************************************************************backward_inner_loop
          !****************************************************************************************************
        enddo

        ! stores wavefield in buffers
        if (any_acoustic) then
          b_potential_acoustic_buffer(:,it_of_this_subset) = b_potential_acoustic(:)
          b_potential_dot_dot_acoustic_buffer(:,it_of_this_subset) = b_potential_dot_dot_acoustic(:)
        endif

        if (any_elastic) then
          b_displ_elastic_buffer(:,:,it_of_this_subset) = b_displ_elastic(:,:)
          b_accel_elastic_buffer(:,:,it_of_this_subset) = b_accel_elastic(:,:)
        endif

      enddo ! subset loop
!****************************************************************************************************backward
!****************************************************************************************************
      it = it_temp
      seismo_current = seismo_current_temp
      ! adjoint wavefield simulation
      do it_of_this_subset = 1, NT_DUMP_ATTENUATION
        ! reads backward/reconstructed wavefield from buffers
        ! note: uses wavefield at corresponding time (NSTEP - it + 1 ), i.e. we have now time-reversed wavefields
        if (any_acoustic) then
          do j = 1,NGLOB
            b_potential_acoustic(j) = b_potential_acoustic_buffer(j,NT_DUMP_ATTENUATION-it_of_this_subset+1)
            b_potential_dot_dot_acoustic(j) = b_potential_dot_dot_acoustic_buffer(j,NT_DUMP_ATTENUATION-it_of_this_subset+1)
          enddo
        endif

        if (any_elastic) then
          do i = 1,3
            do j = 1,NGLOB
              b_displ_elastic(i,j) = b_displ_elastic_buffer(i,j,NT_DUMP_ATTENUATION-it_of_this_subset+1)
              b_accel_elastic(i,j) = b_accel_elastic_buffer(i,j,NT_DUMP_ATTENUATION-it_of_this_subset+1)
            enddo
          enddo
        endif

        it = it + 1

        do i_stage = 1, stage_time_scheme
!****************************************************************************************************adjoint
!****************************************************************************************************
          ! *********************************************************
          ! ************* update_displacement_newmark
          ! *********************************************************
          if (any_acoustic) then
            ! free surface for an acoustic medium
            if (nelem_acoustic_surface > 0) then
              call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                                 potential_acoustic)
            endif

            call update_displacement_newmark_acoustic(deltat,deltatover2,deltatsquareover2,&
                                                      potential_dot_dot_acoustic,potential_dot_acoustic,&
                                                      potential_acoustic,potential_acoustic_old, &
                                                      PML_BOUNDARY_CONDITIONS)

            !ZN here we remove the trick introduced by Luoyang to stabilized the adjoint simulation
            !ZN However in order to keep the current code be consistent, we still keep potential_acoustic_adj_coupling
            !ZN the final goal should remove the *adj_coupling
            potential_acoustic_adj_coupling = potential_acoustic

          endif

          if (any_elastic) then
            call update_displacement_newmark_elastic(deltat,deltatover2,deltatsquareover2,&
                                                     accel_elastic,veloc_elastic,&
                                                     displ_elastic,displ_elastic_old,&
                                                     PML_BOUNDARY_CONDITIONS)
          endif
! *********************************************************
! ************* main solver for the elastic elements
! *********************************************************
          if (any_elastic) then
            !ZN currently we do not support plane wave source in adjoint inversion
            call compute_forces_viscoelastic(accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old, &
                                             PML_BOUNDARY_CONDITIONS,e1,e11,e13)

          endif !if (any_elastic)
          ! *********************************************************
          ! ************* add coupling with the acoustic side
          ! *********************************************************
          if (coupled_acoustic_elastic) then
            call compute_coupling_viscoelastic_ac()
          endif

          ! ************************************************************************************
          ! ************************************ add force source
          ! ************************************************************************************
          if (any_elastic) then
            if (.not. initialfield) then
              call compute_add_sources_viscoelastic_adjoint()
            endif
          endif

          ! ************************************************************************************
          ! ************* assembling accel_elastic for elastic elements
          ! ************************************************************************************
#ifdef USE_MPI
          if (NPROC > 1 .and. any_elastic .and. ninterface_elastic > 0) then
            call assemble_MPI_vector_el(accel_elastic)
          endif
#endif
          ! ************************************************************************************
          ! ************* multiply by the inverse of the mass matrix and update velocity
          ! ************************************************************************************
          if (any_elastic) then
            !! DK DK this should be vectorized
            accel_elastic(1,:) = accel_elastic(1,:) * rmass_inverse_elastic_one(:)
            accel_elastic(2,:) = accel_elastic(2,:) * rmass_inverse_elastic_one(:)
            accel_elastic(3,:) = accel_elastic(3,:) * rmass_inverse_elastic_three(:)

            veloc_elastic = veloc_elastic + deltatover2*accel_elastic
          endif

          if (coupled_acoustic_elastic) then
#ifdef FORCE_VECTORIZATION
              do i = 1,3*nglob_elastic
                accel_elastic_adj_coupling(i,1) = accel_elastic(i,1)
              enddo
#else
              accel_elastic_adj_coupling = accel_elastic
#endif
          endif
! *********************************************************
! ************* main solver for the acoustic elements
! *********************************************************
          if (any_acoustic) then
            call compute_forces_acoustic(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                         potential_acoustic,potential_acoustic_old,PML_BOUNDARY_CONDITIONS)
          endif ! end of test if any acoustic element

          ! *********************************************************
          ! ************* add acoustic forcing at a rigid boundary
          ! *********************************************************
          if (any_acoustic) then
            if (ACOUSTIC_FORCING) then
              call add_acoustic_forcing_at_rigid_boundary(potential_dot_dot_acoustic)
            endif
          endif ! end of test if any acoustic element

          ! *********************************************************
          ! ************* add coupling with the elastic side
          ! *********************************************************
          if (coupled_acoustic_elastic) then
            if (SIMULATION_TYPE == 3) then
              accel_elastic_adj_coupling2 = - accel_elastic_adj_coupling
              call compute_coupling_acoustic_el(accel_elastic_adj_coupling2,displ_elastic_old,potential_dot_dot_acoustic)
            endif
          endif
          ! ************************************************************************************
          ! ************************************ add force source
          ! ************************************************************************************
          if (any_acoustic) then
            if (.not. initialfield) then
              call compute_add_sources_acoustic_adjoint()
            endif
          endif
          ! ************************************************************************************
          ! ********** assembling potential_dot_dot for acoustic elements
          ! ************************************************************************************
#ifdef USE_MPI
          if (NPROC > 1 .and. any_acoustic .and. ninterface_acoustic > 0) then
            call assemble_MPI_vector_ac(potential_dot_dot_acoustic)
          endif
#endif
          ! ************************************************************************************
          ! ************* multiply by the inverse of the mass matrix and update velocity
          ! ************************************************************************************

          if (any_acoustic) then
            ! free surface for an acoustic medium
            if (nelem_acoustic_surface > 0) then
              call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                                 potential_acoustic)
            endif

            if (time_stepping_scheme == 1) then
              !! DK DK this should be vectorized
              potential_dot_dot_acoustic = potential_dot_dot_acoustic * rmass_inverse_acoustic
              potential_dot_acoustic = potential_dot_acoustic + deltatover2*potential_dot_dot_acoustic
            endif
          endif ! of if (any_acoustic)

!****************************************************************************************************adjoint
!****************************************************************************************************
        enddo ! i_stage

        ! computes kinetic and potential energy
        if (output_energy) then
          call it_compute_and_output_energy()
        endif

        ! display time step and max of norm of displacement
        if (mod(it,NSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then
          call check_stability()
        endif

        ! loop on all the receivers to compute and store the seismograms
        call write_seismograms()

        ! kernels calculation
        if (SIMULATION_TYPE == 3) then
          call compute_kernels()
        endif

        ! display results at given time steps
        call write_movie_output()

      enddo ! subset loop
!****************************************************************************************************
!****************************************************************************************************
    end select ! SIMULATION_TYPE

  enddo   ! end of main time loop

  !
  !---- end of time iteration loop
  !

  ! frees undo_attenuation buffers
  if (SIMULATION_TYPE == 3) then
    if (any_acoustic) then
      deallocate(b_potential_acoustic,b_potential_dot_dot_acoustic)
    endif

    if (any_elastic) then
      deallocate(b_displ_elastic,b_accel_elastic)
    endif

  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !----  formats
  400 format(/1x,41('=')/,' =  T i m e  e v o l u t i o n  l o o p  ='/1x,41('=')/)

  end subroutine iterate_time_undoatt
