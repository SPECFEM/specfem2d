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


  subroutine iterate_time()

  use constants, only: IMAIN,NOISE_SAVE_EVERYWHERE
  use specfem_par
  use specfem_par_movie, only: MOVIE_SIMULATION

  implicit none

  ! local parameters

! time the code in order to compute element weights, which can then be used to define typical weights
! for the SCOTCH domain decomposer (this is crucial, in particular in the presence of PML elements
! and/or for mixed simulations, for instance fluid/solid and so on)
  logical, parameter :: TIME_THE_COST_TO_COMPUTE_WEIGHTS_FOR_THE_DOMAIN_DECOMPOSER = .false.

  ! time
  character(len=8) :: datein
  character(len=10) :: timein
  character(len=5) :: zone
  integer, dimension(8) :: time_values
  integer :: year,mon,day,hr,minutes,timestamp

  ! timing
  double precision,external :: wtime
  double precision :: start_time_of_time_loop,duration_of_time_loop_in_seconds

  if (myrank == 0) write(IMAIN,400) ! Write = T i m e  e v o l u t i o n  l o o p =
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

  ! synchronize all processes to make sure everybody is ready to start time loop
  call synchronize_all()

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Starting time iteration loop ...'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! time loop increments begin/end
  it_begin = 1
  it_end = NSTEP

  ! initialize variables for writing seismograms
  seismo_offset = it_begin - 1
  seismo_current = 0

  ! timing safety checks
  if (TIME_THE_COST_TO_COMPUTE_WEIGHTS_FOR_THE_DOMAIN_DECOMPOSER) then
    if (NPROC /= 1) &
      call exit_MPI(myrank,'timing for element weights should be done in serial mode')
    if (NSTEP < 1000) &
      call exit_MPI(myrank,'timing for element weights should be done with at least 1000 time steps')
    if (NSTEP > 10000) &
      call exit_MPI(myrank,'timing for element weights does not need to be done with more than 10000 time steps')
    if (NTSTEP_BETWEEN_OUTPUT_INFO < NSTEP / 5) &
      call exit_MPI(myrank,'timing for element weights should be done with NTSTEP_BETWEEN_OUTPUT_INFO not smaller than NSTEP / 5')
    if (SIMULATION_TYPE /= 1) &
      call exit_MPI(myrank,'timing for element weights should be done with SIMULATION_TYPE = 1')
    if (.not. P_SV) &
      call exit_MPI(myrank,'timing for element weights should be done with P_SV set to true')
    if (GPU_MODE) &
      call exit_MPI(myrank,'timing for element weights should be done with GPU_MODE turned off')
    if (save_ASCII_seismograms) &
      call exit_MPI(myrank,'timing for element weights should be done with save_ASCII_seismograms turned off')
    if (save_binary_seismograms_single) &
      call exit_MPI(myrank,'timing for element weights should be done with save_binary_seismograms_single turned off')
    if (save_binary_seismograms_double) &
      call exit_MPI(myrank,'timing for element weights should be done with save_binary_seismograms_double turned off')
    if (output_color_image) &
      call exit_MPI(myrank,'timing for element weights should be done with output_color_image turned off')
    if (output_postscript_snapshot) &
      call exit_MPI(myrank,'timing for element weights should be done with output_postscript_snapshot turned off')
    if (output_wavefield_dumps) &
      call exit_MPI(myrank,'timing for element weights should be done with output_wavefield_dumps turned off')
    if (OUTPUT_ENERGY) &
      call exit_MPI(myrank,'timing for element weights should be done with OUTPUT_ENERGY turned off')
    if (COMPUTE_INTEGRATED_ENERGY_FIELD) &
      call exit_MPI(myrank,'timing for element weights should be done with COMPUTE_INTEGRATED_ENERGY_FIELD turned off')
  endif

  ! timing
  start_time_of_time_loop = wtime()

! *********************************************************
! ************* MAIN LOOP OVER THE TIME STEPS *************
! *********************************************************

  do it = it_begin,it_end
    ! compute current time
    current_timeval = (it-1) * DT

    ! display time step and max of norm of displacement
    if (mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then
      call check_stability()
    endif

    do i_stage = 1, NSTAGE_TIME_SCHEME

      ! updates wavefields
      select case(time_stepping_scheme)
      case (1)
        ! Newmark time scheme
        ! forward wavefield (acoustic/elastic/poroelastic)
        call update_displ_Newmark()
        ! reconstructed/backward wavefields
        if (SIMULATION_TYPE == 3 .and. .not. NO_BACKWARD_RECONSTRUCTION) call update_displ_Newmark_backward()

      case (2,3)
        ! LDDRK, RK4
        if (any_acoustic) then
          if (.not. GPU_MODE) then
            potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL
            if (SIMULATION_TYPE == 3 .and. .not. NO_BACKWARD_RECONSTRUCTION) b_potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL
          endif
        endif
        if (any_elastic) then
          if (.not. GPU_MODE) then
            accel_elastic(:,:) = 0._CUSTOM_REAL
            if (SIMULATION_TYPE == 3 .and. .not. NO_BACKWARD_RECONSTRUCTION) b_accel_elastic(:,:) = 0._CUSTOM_REAL
          endif
        endif
        if (any_poroelastic) then
          if (.not. GPU_MODE) then
            accels_poroelastic(:,:) = 0._CUSTOM_REAL
            accelw_poroelastic(:,:) = 0._CUSTOM_REAL
            if (SIMULATION_TYPE == 3 .and. .not. NO_BACKWARD_RECONSTRUCTION) then
              b_accels_poroelastic(:,:) = 0._CUSTOM_REAL
              b_accelw_poroelastic(:,:) = 0._CUSTOM_REAL
            endif
          endif
        endif

      case (4)
        ! symplectic PEFRL
        ! forward wavefield (acoustic/elastic/poroelastic)
        call update_displ_symplectic()
        ! reconstructed/backward wavefields
        if (SIMULATION_TYPE == 3 .and. .not. NO_BACKWARD_RECONSTRUCTION) call update_displ_symplectic_backward()

      case default
        call stop_the_code('Error time scheme not implemente yet in iterate_time()')
      end select


      ! note: the order of the computations for acoustic and elastic domains is crucial for coupled simulations.
      !
      !       the coupling terms involve continuity of traction, i.e., the exchange of pressure from fluid
      !       to the solid traction (added in the elastic domain),
      !       and continuity of the normal displacement, i.e., the exchange of the normal component of displacement from solid
      !       to the fluid (added in the acoustic domain).
      !
      !       for forward simulations:
      !         pressure p becomes
      !           p = - \partial_t^2 chi
      !         due to the definition of potential chi and displacement u:
      !           u = 1/rho grad(chi)
      !
      !         given n is the normal to the fluid-solid interface,
      !         the coupling term for elastic media involves: n * T = - p_fluid n = \partial_t^2 chi n
      !         the coupling term for acoustic media involves: 1/rho grad(chi) * n = n * u
      !
      !         where T is the stress tensor, p_fluid pressure, u displacement, chi the acoustic scalar potential,
      !         and rho is density.
      !
      !         \partial_t^2 chi is array potential_dot_dot_acoustic and only available after the acoustic domain update.
      !         thus, by using a Newmark time scheme the update order becomes important (see Chaljub & Valette 2004):
      !         1. acoustic domain update (using displ_accel for the coupling)
      !         2. elastic domain update  (using potential_dot_dot_acoustic for the coupling)
      !         (3. poroelastic domain update with same expression as for elastic domains)
      !
      !       for kernel simulations:
      !         using Lagrange multiplier optimization as in Luo et al. (2013), one finds that once the definition of the forward
      !         acoustic potential is chosen as u = 1/rho grad(chi), then the coupling of the adjoint
      !         wavefields between fluid-solid interfaces becomes different:
      !
      !         the coupling term for elastic media involves: n * T^adj = - p_fluid^adj n = - chi^adj n
      !         the coupling term of acoustic media involves: 1/rho grad(chi^adj) * n = - n * \partial_t^2 u^adj
      !
      !         where T^adj is the adjoint stress tensor, p_fluid^adj the adjoint pressure in the fluid.
      !
      !         note that the adjoint pressure p^adj for the adjoint wavefield corresponds to
      !           p^adj = chi^adj
      !         with an adjoint potential definition chi^adj related to adjoint acceleration (instead of displacement)
      !           \partial_t^2 u^adj = - 1/rho grad(chi^adj)
      !
      !         since elastic media now couples with potential_acoustic rather than potential_dot_dot_acoustic, the order
      !         of the domain updates switches to:
      !         1. elastic domain update  (using -potential_acoustic for the coupling) to have \partial_t^2 u^adj at time n+1
      !         (1b. poroelastic domain update uses same terms as elastic domain)
      !         2. acoustic domain update (using -accel_elastic for the coupling)
      !
      !       for purely adjoint simulations:
      !         if we choose a pure adjoint simulation, without forward wavefield propagation for kernel computations,
      !         then we are free to choose and define the adjoint acoustic potential.
      !         the coupling terms between adjoint wavefields then depend on the chosen potential definition, i.e.,
      !         (i)  if we choose a displacement potential:
      !                u^adj = 1/rho grad(chi^adj)
      !              then the adjoint pressure is p^adj = - potential_dot_dot_acoustic and we have the same coupling terms
      !              as for the forward case.
      !
      !         (ii) if we choose an acceleration potential:
      !                \partial_t^2 u^adj = - 1/rho grad(chi^adj)
      !              then the adjoint pressure is p^adj = chi^adj and we have the same coupling terms as for kernel simulations.
      !
      !       here, we choose to have (i) for a pure adjoint simulations, thus the coupling needs only a special treatment
      !       for kernel simulations.
      !
      !       and, if there is no coupling between different media, then the ordering is obviously not important
      !       as the domains can be computed independently from each other.

      select case(SIMULATION_TYPE)
      case (1,2)
        ! forward/adjoint simulations
        ! acoustic domains
        if (ACOUSTIC_SIMULATION) then
          if (.not. GPU_MODE) then
            call compute_forces_viscoacoustic_main()
          else
            ! on GPU
            if (any_acoustic) call compute_forces_viscoacoustic_GPU(.false.)
          endif
        endif
        ! elastic domains
        if (ELASTIC_SIMULATION) then
          if (.not. GPU_MODE) then
            call compute_forces_viscoelastic_main()
          else
            ! on GPU
            if (any_elastic) call compute_forces_viscoelastic_GPU(.false.)
          endif
        endif
        ! poroelastic domains
        if (POROELASTIC_SIMULATION) then
          if (.not. GPU_MODE) then
            call compute_forces_poroelastic_main()
          else
            ! on GPU
            if (any_poroelastic) call exit_MPI(myrank,'poroelastic not implemented in GPU MODE yet')
          endif
        endif

      case (3)
        ! kernel simulations
        ! forward (adjoint) wavefields use switched update ordering
        ! poroelastic domains
        if (POROELASTIC_SIMULATION) then
          if (.not. GPU_MODE) then
            call compute_forces_poroelastic_main()
          else
            ! on GPU
            if (any_poroelastic) call exit_MPI(myrank,'poroelastic not implemented in GPU MODE yet')
          endif
        endif
        ! elastic domains
        if (ELASTIC_SIMULATION) then
          if (.not. GPU_MODE) then
            call compute_forces_viscoelastic_main()
          else
            ! on GPU
            if (any_elastic) call compute_forces_viscoelastic_GPU(.false.)
          endif
        endif
        ! acoustic domains
        if (.not. GPU_MODE) then
          call compute_forces_viscoacoustic_main()
        else
          ! on GPU
          if (any_acoustic) call compute_forces_viscoacoustic_GPU(.false.)
        endif

        ! backward/reconstructed wavefields
        if (.not. NO_BACKWARD_RECONSTRUCTION) then
          ! note: the reconstruction of the forward wavefield uses the last wavefield and backpropagates it in time.
          !       backpropagation can use the same coupling terms as in the forward propagation case, thus the same
          !       order as for forward simulations.
          if (.not. GPU_MODE) then
            ! acoustic domains
            call compute_forces_viscoacoustic_main_backward()
            ! elastic domains
            call compute_forces_viscoelastic_main_backward()
            ! poroelastic domains
            call compute_forces_poroelastic_main_backward()
          else
            ! on GPU
            if (coupled_acoustic_elastic) then
              ! coupled simulations need a re-ordering of the domain updates
              ! due to a change of the coupling terms for adjoint wavefields
              !
              ! we don't need to call this for un-coupled simulation, as both wavefields will be taken take of in the
              ! GPU routines. only when the coupling terms for adjoint simulations are involved, we need to separate
              ! forward/adjoint and backward wavefield computations.

              ! here, we still need to compute reconstructed/backward wavefields
              ! acoustic domains
              if (any_acoustic) call compute_forces_viscoacoustic_GPU(.true.)  ! only reconstructed/backward wavefields
              ! elastic domains
              if (any_elastic) call compute_forces_viscoelastic_GPU(.true.)   ! only reconstructed/backward wavefields
              ! poroelastic domains not implemented on GPUs yet...
              if (any_poroelastic) call exit_MPI(myrank,'poroelastic coupling not implemented in GPU MODE yet')
            endif
          endif
        endif
      end select

    enddo ! NSTAGE_TIME_SCHEME (LDDRK or RK)

    ! reads in lastframe for adjoint/kernels calculation
    if (SIMULATION_TYPE == 3 .and. it == 1 .and. .not. NO_BACKWARD_RECONSTRUCTION) then
      call read_forward_arrays()
    endif

    if (NO_BACKWARD_RECONSTRUCTION .and. (SAVE_FORWARD .or. SIMULATION_TYPE == 3) ) then
      call manage_no_backward_reconstruction_io()
    endif

    if (OUTPUT_ENERGY) then
      call compute_and_output_energy()
    endif

    ! computes integrated_energy_field
    if (COMPUTE_INTEGRATED_ENERGY_FIELD) then
      ! compute the field int_0^t v^2 dt and write it on file if needed
      call compute_integrated_energy_field_and_output()
    endif

    ! loop on all the receivers to compute and store the seismograms
    call write_seismograms()

    ! kernels calculation
    if (SIMULATION_TYPE == 3) then
      call compute_kernels()
    endif

    ! display results at given time steps
    if (MOVIE_SIMULATION) then
      call write_movie_output(.true.)
    endif

    ! first step of noise tomography, i.e., save a surface movie at every time step
    if (NOISE_TOMOGRAPHY == 1) then
      call noise_save_surface_movie()
    endif

    ! noise simulations
    if (NOISE_SAVE_EVERYWHERE) then
      select case (NOISE_TOMOGRAPHY)
      case (2)
        ! stores complete wavefield for reconstruction
        call noise_save_surface_movie()
      case (3)
        ! reconstructs forward wavefield based on complete wavefield storage
        call noise_read_wavefield()
      end select
    endif

  enddo ! end of the main time loop

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

  call date_and_time(datein,timein,zone,time_values)
  year = time_values(1)
  mon = time_values(2)
  day = time_values(3)
  hr = time_values(5)
  minutes = time_values(6)
  call convtime(timestamp,year,mon,day,hr,minutes)

  if (myrank == 0) then
   write(IMAIN,*)
   write(IMAIN,*) 'Total duration of the timeloop in seconds, measured using '
   write(IMAIN,*)  'date and time of the system : ',sngl(timestamp*60.d0 + time_values(7) + &
                   time_values(8)/1000.d0 - timestamp_seconds_start),' s'
  endif

! *********************************************************
! ************* END MAIN LOOP OVER THE TIME STEPS *********
! *********************************************************!

  ! Transfer fields from GPU card to host for further analysis
  if (GPU_MODE) call it_transfer_from_GPU()

  !----  formats
  400 format(/1x,41('=')/,' =  T i m e  e v o l u t i o n  l o o p  ='/1x,41('=')/)

  end subroutine iterate_time

!
!----------------------------------------------------------------------------------------
!

  subroutine it_transfer_from_GPU()

! transfers fields on GPU back onto CPU

  use constants, only: TWO,FOUR_THIRDS
  use specfem_par
  use specfem_par_gpu, only: Mesh_pointer,NGLOB_AB,NSPEC_AB

  implicit none

  ! Fields transfer for imaging
  ! acoustic domains
  if (any_acoustic ) then
    call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic,potential_dot_acoustic, &
                                        potential_dot_dot_acoustic,Mesh_pointer)
  endif

  ! elastic domains
  if (any_elastic) then
    call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ_elastic,veloc_elastic,accel_elastic,Mesh_pointer)
  endif

  ! finishes kernel calculations
  if (SIMULATION_TYPE == 3) then
    ! Kernel transfer
    ! acoustic domains
    if (any_acoustic) then
      call transfer_kernels_ac_to_host(Mesh_pointer,rho_ac_kl,kappa_ac_kl,NSPEC_AB)

      ! note: acoustic kernels in CPU-version add (delta * NTSTEP_BETWEEN_COMPUTE_KERNELS) factors at each computation.
      !       for the GPU-kernels, the NSTEP_** is still missing... adding it here
      rho_ac_kl(:,:,:) = rho_ac_kl(:,:,:) * NTSTEP_BETWEEN_COMPUTE_KERNELS
      kappa_ac_kl(:,:,:) = kappa_ac_kl(:,:,:) * NTSTEP_BETWEEN_COMPUTE_KERNELS

      rhop_ac_kl(:,:,:) = rho_ac_kl(:,:,:) + kappa_ac_kl(:,:,:)
      alpha_ac_kl(:,:,:) = TWO *  kappa_ac_kl(:,:,:)
    endif

    ! elastic domains
    if (any_elastic) then
      ! note: kernel values in the elastic case will be multiplied with material properties
      !       in save_adjoint_kernels.f90
      call transfer_kernels_el_to_host(Mesh_pointer,rho_kl,mu_kl,kappa_kl,NSPEC_AB)
    endif
  endif

  end subroutine it_transfer_from_GPU


!
!----------------------------------------------------------------------------------------
!

  subroutine manage_no_backward_reconstruction_io()

  use constants
  use specfem_par
  use specfem_par_gpu

  implicit none

  ! EB EB June 2018 : this routine is saving/reading the frames of the forward
  ! wavefield that are necessary to compute the sensitivity kernels

  ! safety check
  if (.not. NO_BACKWARD_RECONSTRUCTION) return

  ! writes in frame for further adjoint/kernels calculation
  if (SAVE_FORWARD) then

    if (mod(NSTEP-it+1,NTSTEP_BETWEEN_COMPUTE_KERNELS) == 0 .or. it == NSTEP ) then

      call save_forward_arrays_no_backward()

    endif

  endif ! SAVE_FORWARD

  if (SIMULATION_TYPE == 3) then

    if (mod(it,NTSTEP_BETWEEN_COMPUTE_KERNELS) == 0 .or. it == 1 .or. it == NSTEP) then

      call read_forward_arrays_no_backward()

      ! For the first iteration we need to add a transfer, to initiate the
      ! transfer on the GPU
      if (no_backward_iframe == 1) &
        call read_forward_arrays_no_backward()

      ! If the wavefield is requested at it=1, then we enforce another transfer
      if (no_backward_iframe == 2 .and. NTSTEP_BETWEEN_COMPUTE_KERNELS == 1) &
        call read_forward_arrays_no_backward()

    endif

  endif ! SIMULATION == 3

  end subroutine  manage_no_backward_reconstruction_io
