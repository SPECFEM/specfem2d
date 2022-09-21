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


  subroutine iterate_time_undoatt()

  use constants, only: IMAIN, USE_A_STRONG_FORMULATION_FOR_E1
  use specfem_par
  use specfem_par_gpu, only: Mesh_pointer
  use specfem_par_movie, only: MOVIE_SIMULATION

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_potential_acoustic_buffer
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_displ_elastic_buffer,b_accel_elastic_buffer
  double precision :: sizeval

  integer :: it_temp,iframe_kernel,nframes_kernel,size_buffer
  integer :: ier

  ! time
  character(len=8) :: datein
  character(len=10) :: timein
  character(len=5) :: zone
  integer, dimension(8) :: time_values
  integer :: year,mon,day,hr,minutes,timestamp

  ! timing
  double precision, external :: wtime
  double precision :: start_time_of_time_loop,duration_of_time_loop_in_seconds
  logical :: compute_b_wavefield

  ! checks if anything to do
  if (.not. UNDO_ATTENUATION_AND_OR_PML) return

  ! user output
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
  if (GPU_MODE .and. ATTENUATION_VISCOELASTIC) &
    call exit_MPI(myrank,'For UNDO_ATTENUATION, GPU_MODE does not support ATTENUATION_VISCOELASTIC')
  if (any_poroelastic) &
    call exit_MPI(myrank,'UNDO_ATTENUATION is not implemented for poroelastic simulation yet')
  if (NOISE_TOMOGRAPHY /= 0) &
    call exit_MPI(myrank,'For UNDO_ATTENUATION, NOISE_TOMOGRAPHY is not supported')
  if (ATTENUATION_VISCOACOUSTIC .and. .not. USE_A_STRONG_FORMULATION_FOR_E1) &
    call exit_MPI(myrank,'For UNDO_ATTENUATION with viscoacousticity, &
                 & USE_A_STRONG_FORMULATION_FOR_E1 has to be set to true (in setup/constants.h)')
  if (NOISE_TOMOGRAPHY == 3) &
    call stop_the_code('UNDO_ATTENUATION for noise kernels not implemented yet')

  ! number of time subsets for time loop
  NSUBSET_ITERATIONS = ceiling( dble(NSTEP)/dble(NT_DUMP_ATTENUATION) )

  ! get the maximum number of frames to save
  if (NTSTEP_BETWEEN_COMPUTE_KERNELS == 1) then
    size_buffer = NT_DUMP_ATTENUATION
  else
    size_buffer = NT_DUMP_ATTENUATION / NTSTEP_BETWEEN_COMPUTE_KERNELS + 1
  endif

  ! checks
  if (NSUBSET_ITERATIONS <= 0) call exit_MPI(myrank,'Error invalid number of time subsets for undoing attenuation')

  ! user output
  if (SAVE_FORWARD .or. SIMULATION_TYPE == 3) then
    if (myrank == 0) then
      write(IMAIN,*) 'undoing attenuation:'
      write(IMAIN,*) '  total number of time subsets                     = ',NSUBSET_ITERATIONS
      write(IMAIN,*) '  wavefield snapshots at every NT_DUMP_ATTENUATION = ',NT_DUMP_ATTENUATION
      if (any_acoustic) then
        ! buffer(nglob,NT_DUMP_ATTENUATION) in MB
        sizeval = dble(nglob) * dble(size_buffer) * dble(CUSTOM_REAL) / 1024.d0 / 1024.d0
        write(IMAIN,*) '  size of acoustic wavefield buffer per slice      = ', sngl(sizeval),'MB'
      endif
      if (any_elastic) then
        ! buffer(3,nglob,NT_DUMP_ATTENUATION) in MB
        sizeval = dble(NDIM) * dble(nglob) * dble(size_buffer) * dble(CUSTOM_REAL) / 1024.d0 / 1024.d0
        if (APPROXIMATE_HESS_KL) sizeval = 2 * sizeval
        write(IMAIN,*) '  size of elastic wavefield buffer per slice       = ', sngl(sizeval),'MB'
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! allocates buffers
  if (SIMULATION_TYPE == 3) then
    if (any_acoustic) then
      allocate(b_potential_acoustic_buffer(nglob,size_buffer),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'error allocating b_potential_acoustic')
    endif
    if (any_elastic) then
      allocate(b_displ_elastic_buffer(NDIM,nglob,size_buffer),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'error allocating b_displ_elastic')
      if (APPROXIMATE_HESS_KL) allocate(b_accel_elastic_buffer(NDIM,nglob,size_buffer),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'error allocating b_accel_elastic')
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

  ! time loop increments begin/end
  it_begin = 1
  it_end = NSTEP

  ! initialize variables for writing seismograms
  seismo_offset = it_begin - 1
  seismo_current = 0

  ! *********************************************************
  ! ************* MAIN LOOP OVER THE TIME STEPS *************
  ! *********************************************************

  ! initializes time increments
  it = it_begin - 1

  ! timing
  start_time_of_time_loop = wtime()

  ! loops over time subsets
  do iteration_on_subset = 1, NSUBSET_ITERATIONS

    ! wavefield storage
    if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
      ! saves forward wavefields
      call save_forward_arrays_undoatt()

    else if (SIMULATION_TYPE == 3) then
      ! reads in last stored forward wavefield
      call read_forward_arrays_undoatt()
    endif

    ! time loop within this iteration subset
    select case (SIMULATION_TYPE)
    case (1)
      ! forward and adjoint simulations

      ! increment end of this subset
      if (iteration_on_subset < NSUBSET_ITERATIONS) then
        ! takes full length of subset
        it_subset_end = NT_DUMP_ATTENUATION
      else
        ! loops over remaining steps in last subset
        it_subset_end = NSTEP - (iteration_on_subset-1)*NT_DUMP_ATTENUATION
      endif
      ! checks end index
      if (it_subset_end > NT_DUMP_ATTENUATION) &
        call exit_MPI(myrank,'Error invalid buffer index for undoing attenuation')

      ! subset loop
      ! We don't compute the forward reconstructed wavefield in the loop below
      compute_b_wavefield = .false.

      do it_of_this_subset = 1, it_subset_end

        it = it + 1
        ! compute current time
        current_timeval = (it-1) * DT

        ! display time step and max of norm of displacement
        if (mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) call check_stability()

        do i_stage = 1, NSTAGE_TIME_SCHEME ! is equal to 1 if Newmark because only one stage then

          ! update displacement
          select case(time_stepping_scheme)
          case (1)
            ! Newmark
            if (any_acoustic) then
              if (GPU_MODE) then
                call update_displacement_newmark_GPU_acoustic(compute_b_wavefield)
              else
                call update_displ_acoustic_forward()
              endif
            endif
            if (any_elastic) then
              call update_displ_elastic_forward()
            endif
            if (any_poroelastic) then
              call update_displ_poroelastic_forward()
            endif

          case (2,3)
            ! LDDRK, RK4
            if (any_acoustic) then
              if (.not. GPU_MODE) potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL
            endif
            if (any_elastic) then
              if (.not. GPU_MODE) accel_elastic(:,:) = 0._CUSTOM_REAL
            endif
            if (any_poroelastic) then
              if (.not. GPU_MODE) then
                accels_poroelastic(:,:) = 0._CUSTOM_REAL
                accelw_poroelastic(:,:) = 0._CUSTOM_REAL
              endif
            endif

          case (4)
            ! symplectic PEFRL
            if (any_acoustic) then
              call update_displ_symplectic_acoustic_forward()
            endif
            if (any_elastic) then
              call update_displ_symplectic_elastic_forward()
            endif
            if (any_poroelastic) then
              call update_displ_symplectic_poroelastic_forward()
            endif

          case default
            call stop_the_code('Error time scheme for update displacement not implemented yet in iterate_time_undoatt()')
          end select

          ! acoustic domains
          if (any_acoustic) then
            if (.not. GPU_MODE) then
              call compute_forces_viscoacoustic_main()
            else
              ! on GPU
              call compute_forces_viscoacoustic_GPU(compute_b_wavefield)
            endif
          endif

          ! elastic domains
          if (any_elastic) then
            if (.not. GPU_MODE) then
              call compute_forces_viscoelastic_main()
            else
              ! on GPU
              call compute_forces_viscoelastic_GPU(compute_b_wavefield)
            endif
          endif

          if (any_poroelastic) call exit_MPI(myrank,'UNDO_ATTENUATION is not implemented for poroelastic simulation yet')
        enddo ! i_stage

        ! computes kinetic and potential energy
        if (OUTPUT_ENERGY) then
          call compute_and_output_energy()
        endif

        ! loop on all the receivers to compute and store the seismograms
        call write_seismograms()

        ! display results at given time steps
        if (MOVIE_SIMULATION) then
          call write_movie_output(compute_b_wavefield)
        endif

        ! first step of noise tomography, i.e., save a surface movie at every time step
        if (NOISE_TOMOGRAPHY == 1) then
          call noise_save_surface_movie()
        endif

      enddo ! subset loop

    case (3)
      ! kernel simulations
      ! intermediate storage of it position
      it_temp = it
      it = NT_DUMP_ATTENUATION * (NSUBSET_ITERATIONS - iteration_on_subset)
      ! increment end of this subset
      if (iteration_on_subset == 1) then
        ! loops over remaining steps in last forward subset
        it_subset_end = NSTEP - (NSUBSET_ITERATIONS-1)*NT_DUMP_ATTENUATION
      else
        ! takes full length of subset
        it_subset_end = NT_DUMP_ATTENUATION
      endif
      ! checks end index
      if (it_subset_end > NT_DUMP_ATTENUATION) &
        call exit_MPI(myrank,'Error invalid buffer index for undoing attenuation')

      ! reconstructs forward wavefield based on last stored wavefield data
      !
      ! note: we step forward in time here, starting from last snapshot.
      !       the newly computed, reconstructed forward wavefields (b_displ_..) get stored in buffers.

      iframe_kernel = 0

      ! subset loop
      do it_of_this_subset = 1, it_subset_end

        it = it + 1
        ! We compute the forward reconstructed wavefield first
        compute_b_wavefield = .true.

        do i_stage = 1, NSTAGE_TIME_SCHEME
          ! backward_inner_loop
          ! update displacement
          select case(time_stepping_scheme)
          case (1)
            ! Newmark
            if (GPU_MODE) then
              if (any_acoustic) call update_displacement_newmark_GPU_acoustic(compute_b_wavefield)
            else
              call update_displ_Newmark_backward()
            endif

          case (2,3)
            ! LDDRK, RK4
            if (any_acoustic) then
              if (.not. GPU_MODE) b_potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL
            endif
            if (any_elastic) then
              if (.not. GPU_MODE) b_accel_elastic(:,:) = 0._CUSTOM_REAL
            endif
            if (any_poroelastic) then
              if (.not. GPU_MODE) then
                b_accels_poroelastic(:,:) = 0._CUSTOM_REAL
                b_accelw_poroelastic(:,:) = 0._CUSTOM_REAL
              endif
            endif

          case (4)
            ! symplectic PEFRL
            if (any_acoustic) then
              call update_displ_symplectic_acoustic_backward()
            endif
            if (any_elastic) then
              call update_displ_symplectic_elastic_backward()
            endif
            if (any_poroelastic) then
              call update_displ_symplectic_poroelastic_backward()
            endif

          case default
            call stop_the_code('Error time scheme for update displacement not implemented yet in iterate_time_undoatt()')
          end select

          ! acoustic domains
          if (any_acoustic) then
            if (.not. GPU_MODE) then
              call compute_forces_viscoacoustic_main_backward()
            else
              ! on GPU
              call compute_forces_viscoacoustic_GPU(compute_b_wavefield)
            endif
          endif

          ! elastic domains
          if (any_elastic) then
            if (.not. GPU_MODE) then
              call compute_forces_viscoelastic_main_backward()
            else
              ! on GPU
              call compute_forces_viscoelastic_GPU(compute_b_wavefield)
            endif
          endif

          if (any_poroelastic) call exit_MPI(myrank,'UNDO_ATTENUATION is not implemented for poroelastic simulation yet')
        enddo

        ! stores wavefield in buffers
        if (mod(NSTEP-it+1,NTSTEP_BETWEEN_COMPUTE_KERNELS) == 0) then

          iframe_kernel = iframe_kernel + 1

          if (any_acoustic) then
            if (GPU_MODE) then
              call transfer_b_potential_ac_from_device(nglob,b_potential_acoustic_buffer(:,iframe_kernel),Mesh_pointer)
            else
              b_potential_acoustic_buffer(:,iframe_kernel) = b_potential_acoustic(:)
            endif
          endif

          if (any_elastic) then
            b_displ_elastic_buffer(:,:,iframe_kernel) = b_displ_elastic(:,:)
            if (APPROXIMATE_HESS_KL) b_accel_elastic_buffer(:,:,iframe_kernel) = b_accel_elastic(:,:)
          endif

        endif

        ! display results at given time steps
        if (MOVIE_SIMULATION) then
          call write_movie_output(compute_b_wavefield)
        endif

      enddo ! subset loop

      ! resets current it position
      it = it_temp
      nframes_kernel = iframe_kernel
      iframe_kernel = 0

      ! adjoint wavefield simulation
      do it_of_this_subset = 1, it_subset_end

        it = it + 1
        if (mod(it,NTSTEP_BETWEEN_COMPUTE_KERNELS) == 0) then

          iframe_kernel = iframe_kernel + 1

          ! reads backward/reconstructed wavefield from buffers
          ! note: uses wavefield at corresponding time (NSTEP - it + 1 ), i.e. we have now time-reversed wavefields
          if (any_acoustic) then
            if (GPU_MODE) then
              call transfer_b_potential_ac_to_device(nglob, &
                                        b_potential_acoustic_buffer(:,nframes_kernel-iframe_kernel+1),Mesh_pointer)
            else
              b_potential_acoustic(:) = b_potential_acoustic_buffer(:,nframes_kernel-iframe_kernel+1)
            endif
          endif

          ! copy the reconstructed wavefield for kernel integration
          if (any_elastic) then
             b_displ_elastic(:,:) = b_displ_elastic_buffer(:,:,nframes_kernel-iframe_kernel+1)
            if (APPROXIMATE_HESS_KL) then
               b_accel_elastic(:,:) = b_accel_elastic_buffer(:,:,nframes_kernel-iframe_kernel+1)
            endif
          endif

        endif ! mod(it,NTSTEP_BETWEEN_COMPUTE_KERNELS) == 0

        ! compute current time
        current_timeval = (it-1) * DT

        ! display time step and max of norm of displacement
        if ((.not. GPU_MODE) .and. mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then
          call check_stability()
        endif

        ! we only compute the adjoint wavefield on the next loop
        compute_b_wavefield = .false.

        do i_stage = 1, NSTAGE_TIME_SCHEME
          ! adjoint
          ! update displacement
          select case(time_stepping_scheme)
          case (1)
            ! Newmark
            if (any_acoustic) then
              if (GPU_MODE) then
                call update_displacement_newmark_GPU_acoustic(compute_b_wavefield)
              else
                call update_displ_acoustic_forward()
              endif
            endif
            if (any_elastic) then
              call update_displ_elastic_forward()
            endif
            if (any_poroelastic) then
              call update_displ_poroelastic_forward()
            endif
          case (2,3)
            ! LDDRK, RK4
            if (any_acoustic) then
              if (.not. GPU_MODE) potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL
            endif
            if (any_elastic) then
              if (.not. GPU_MODE) accel_elastic(:,:) = 0._CUSTOM_REAL
            endif
            if (any_poroelastic) then
              if (.not. GPU_MODE) then
                accels_poroelastic(:,:) = 0._CUSTOM_REAL
                accelw_poroelastic(:,:) = 0._CUSTOM_REAL
              endif
            endif

          case (4)
            ! symplectic PEFRL
            if (any_acoustic) then
              call update_displ_symplectic_acoustic_forward()
            endif
            if (any_elastic) then
              call update_displ_symplectic_elastic_forward()
            endif
            if (any_poroelastic) then
              call update_displ_symplectic_poroelastic_forward()
            endif

          case default
            call stop_the_code('Error time scheme for update displacement not implemented yet in iterate_time_undoatt()')
          end select

          ! not needed anymore, taking care of by re-ordering domain updates
          !if (any_acoustic) then
          !  !ZN here we remove the trick introduced by Luoyang to stabilized the adjoint simulation
          !  !ZN However in order to keep the current code be consistent, we still keep potential_acoustic_adj_coupling
          !  !ZN the final goal should remove the *adj_coupling
          !  potential_acoustic_adj_coupling(:) = potential_acoustic(:)
          !endif

          ! note: for adjoint wavefields, the coupling terms change and the order of the domain updates should be:
          !       1. elastic domain (and/or poroelastic domain)
          !       2. acoustic domain

          if (any_poroelastic) call exit_MPI(myrank,'UNDO_ATTENUATION is not implemented for poroelastic simulation yet')

          ! main solver for the elastic elements
          if (any_elastic) then
            if (.not. GPU_MODE) then
              call compute_forces_viscoelastic_main()
            else
              ! on GPU
              call compute_forces_viscoelastic_GPU(compute_b_wavefield)
            endif
          endif

          ! for coupling with adjoint wavefields, stores temporary old wavefield
          ! not needed anymore, taking care of by re-ordering domain updates
          !if (coupled_acoustic_elastic) then
          !  ! handles adjoint runs coupling between adjoint potential and adjoint elastic wavefield
          !  ! adjoint definition: \partial_t^2 \bfs^\dagger = - \frac{1}{\rho} \bfnabla \phi^\dagger
          !  accel_elastic_adj_coupling(:,:) = - accel_elastic(:,:)
          !endif

          ! acoustic domains
          if (any_acoustic) then
            if (.not. GPU_MODE) then
              call compute_forces_viscoacoustic_main()
            else
              ! on GPU
              call compute_forces_viscoacoustic_GPU(compute_b_wavefield)
            endif
          endif

        enddo !i_stage

        ! computes kinetic and potential energy
        if (OUTPUT_ENERGY) then
          call compute_and_output_energy()
        endif

        ! loop on all the receivers to compute and store the seismograms
        !call write_seismograms()

        ! kernels calculation
        call compute_kernels()

        ! display results at given time steps
        if (MOVIE_SIMULATION) then
          call write_movie_output(compute_b_wavefield)
        endif

      enddo ! subset loop
    end select ! SIMULATION_TYPE

  enddo   ! end of main time loop

  !
  !---- end of time iteration loop
  !

  if (myrank == 0) then
    duration_of_time_loop_in_seconds = wtime() - start_time_of_time_loop

    write(IMAIN,*)
    write(IMAIN,*) 'Total duration of the time loop in seconds = ',sngl(duration_of_time_loop_in_seconds),' s'
    write(IMAIN,*) 'Total number of time steps = ',NSTEP
    write(IMAIN,*) 'Average duration of a time step of the time loop = ',sngl(duration_of_time_loop_in_seconds / NSTEP),' s'
    write(IMAIN,*) 'Total number of spectral elements in the mesh = ',NSPEC
    write(IMAIN,*) '    of which ',NSPEC - count(ispec_is_PML),' are regular elements'
    write(IMAIN,*) '    and ',count(ispec_is_PML),' are PML elements.'
    write(IMAIN,*) 'Average duration of the calculation per spectral element = ', &
                         sngl(duration_of_time_loop_in_seconds / (NSTEP * NSPEC)),' s'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! Transfer fields from GPU card to host for further analysis
  if (GPU_MODE) call it_transfer_from_GPU()


  ! frees undo_attenuation buffers
  if (SIMULATION_TYPE == 3) then
    if (any_acoustic) then
      deallocate(b_potential_acoustic_buffer)
    endif
    if (any_elastic) then
      deallocate(b_displ_elastic_buffer)
      if (APPROXIMATE_HESS_KL) deallocate(b_accel_elastic_buffer)
    endif
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !----  formats
  400 format(/1x,41('=')/,' =  T i m e  e v o l u t i o n  l o o p  ='/1x,41('=')/)

  end subroutine iterate_time_undoatt
