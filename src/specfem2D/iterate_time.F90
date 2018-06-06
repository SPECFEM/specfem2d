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

#ifdef USE_MPI
  use mpi
#endif

  use constants, only: IMAIN,NOISE_SAVE_EVERYWHERE
  use specfem_par
  use specfem_par_gpu
  use specfem_par_noise, only: NOISE_TOMOGRAPHY

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

  real :: start_time_of_time_loop,finish_time_of_time_loop,duration_of_time_loop_in_seconds

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

  ! initialize variables for writing seismograms
  seismo_offset = 0
  seismo_current = 0

  if (TIME_THE_COST_TO_COMPUTE_WEIGHTS_FOR_THE_DOMAIN_DECOMPOSER) then
    if (NPROC /= 1) call exit_MPI(myrank,'timing for element weights should be done in serial mode')
    if (NSTEP < 1000) call exit_MPI(myrank,'timing for element weights should be done with at least 1000 time steps')
    if (NSTEP > 10000) call exit_MPI(myrank,'timing for element weights does not need to be done with more than 10000 time steps')
    if (NSTEP_BETWEEN_OUTPUT_INFO < NSTEP / 5) &
       call exit_MPI(myrank,'timing for element weights should be done with NSTEP_BETWEEN_OUTPUT_INFO not smaller than NSTEP / 5')
    if (SIMULATION_TYPE /= 1) call exit_MPI(myrank,'timing for element weights should be done with SIMULATION_TYPE = 1')
    if (.not. P_SV) call exit_MPI(myrank,'timing for element weights should be done with P_SV set to true')
    if (GPU_MODE) call exit_MPI(myrank,'timing for element weights should be done with GPU_MODE turned off')
    if (save_ASCII_seismograms) &
       call exit_MPI(myrank,'timing for element weights should be done with save_ASCII_seismograms turned off')
    if (save_binary_seismograms_single) &
       call exit_MPI(myrank,'timing for element weights should be done with save_binary_seismograms_single turned off')
    if (save_binary_seismograms_double) &
       call exit_MPI(myrank,'timing for element weights should be done with save_binary_seismograms_double turned off')
    if (output_color_image) call exit_MPI(myrank,'timing for element weights should be done with output_color_image turned off')
    if (output_postscript_snapshot) &
       call exit_MPI(myrank,'timing for element weights should be done with output_postscript_snapshot turned off')
    if (output_wavefield_dumps) &
       call exit_MPI(myrank,'timing for element weights should be done with output_wavefield_dumps turned off')
    if (OUTPUT_ENERGY) call exit_MPI(myrank,'timing for element weights should be done with OUTPUT_ENERGY turned off')
    if (COMPUTE_INTEGRATED_ENERGY_FIELD) &
       call exit_MPI(myrank,'timing for element weights should be done with COMPUTE_INTEGRATED_ENERGY_FIELD turned off')
  endif

  call cpu_time(start_time_of_time_loop)

! *********************************************************
! ************* MAIN LOOP OVER THE TIME STEPS *************
! *********************************************************

  do it = 1,NSTEP
    ! compute current time
    timeval = (it-1) * deltat

    ! display time step and max of norm of displacement
    if (mod(it,NSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then
      call check_stability()
    endif

    do i_stage = 1, stage_time_scheme

      ! updates wavefields using Newmark time scheme
      ! forward wavefield (acoustic/elastic/poroelastic)
      call update_displ_Newmark()
      ! reconstructed/backward wavefields
      if (SIMULATION_TYPE == 3 .and. .not. NO_BACKWARD_RECONSTRUCTION) call update_displ_Newmark_backward()

      ! acoustic domains
      if (ACOUSTIC_SIMULATION) then
        if (.not. GPU_MODE) then
          call compute_forces_viscoacoustic_main()
          if (SIMULATION_TYPE == 3 .and. .not. NO_BACKWARD_RECONSTRUCTION) call compute_forces_viscoacoustic_main_backward()
        else
          ! on GPU
          if (any_acoustic) call compute_forces_viscoacoustic_GPU(.false.)
        endif
      endif

      ! elastic domains
      if (ELASTIC_SIMULATION) then
        if (.not. GPU_MODE) then
          call compute_forces_viscoelastic_main()
          if (SIMULATION_TYPE == 3 .and. .not. NO_BACKWARD_RECONSTRUCTION) call compute_forces_viscoelastic_main_backward()
        else
          ! on GPU
          if (any_elastic) call compute_forces_viscoelastic_GPU()
        endif
      endif

      ! poroelastic domains
      if (POROELASTIC_SIMULATION) then
        if (.not. GPU_MODE) then
          call compute_forces_poroelastic_main()
          if (SIMULATION_TYPE == 3) call compute_forces_poroelastic_main_backward()
        else
          ! on GPU
          if (any_poroelastic) call exit_MPI(myrank,'poroelastic not implemented in GPU MODE yet')
        endif
      endif

    enddo ! stage_time_scheme (LDDRK or RK)

    ! reads in lastframe for adjoint/kernels calculation
    if (SIMULATION_TYPE == 3 .and. it == 1 .and. .not. NO_BACKWARD_RECONSTRUCTION) then
      call read_forward_arrays()
    endif

    if (NO_BACKWARD_RECONSTRUCTION .and. (SAVE_FORWARD .or. SIMULATION_TYPE == 3) ) then
      call manage_no_backward_reconstruction_io()
    endif

    ! noise simulations
    select case (NOISE_TOMOGRAPHY)
    case (1)
      ! stores generating wavefield
      call save_surface_movie_noise()
    case (2)
      ! stores complete wavefield for reconstruction
      if (NOISE_SAVE_EVERYWHERE) call save_surface_movie_noise()
    case (3)
      ! reconstructs forward wavefield based on complete wavefield storage
      if (NOISE_SAVE_EVERYWHERE) call read_wavefield_noise()
    end select

    if (output_energy) then
      call compute_and_output_energy()
    endif

    ! computes integrated_energy_field
    if (COMPUTE_INTEGRATED_ENERGY_FIELD) then
      call it_compute_integrated_energy_field_and_output() ! compute the field int_0^t v^2 dt and write it on file if needed
    endif

    ! loop on all the receivers to compute and store the seismograms
    call write_seismograms()

    ! kernels calculation
    if (SIMULATION_TYPE == 3 .and. mod(it,NSTEP_BETWEEN_COMPUTE_KERNELS) == 0) then
      call compute_kernels()
    endif

    ! display results at given time steps
    call write_movie_output(.true.)

  enddo ! end of the main time loop

  call cpu_time(finish_time_of_time_loop)

  if (myrank == 0) then
    duration_of_time_loop_in_seconds = finish_time_of_time_loop - start_time_of_time_loop
    write(IMAIN,*)
    write(IMAIN,*) 'Total duration of the time loop in seconds = ',duration_of_time_loop_in_seconds,' s'
    write(IMAIN,*) 'Total number of time steps = ',NSTEP
    write(IMAIN,*) 'Average duration of a time step of the time loop = ',duration_of_time_loop_in_seconds / real(NSTEP),' s'
    write(IMAIN,*) 'Total number of spectral elements in the mesh = ',NSPEC
    write(IMAIN,*) '    of which ',NSPEC - count(ispec_is_PML),' are regular elements'
    write(IMAIN,*) '    and ',count(ispec_is_PML),' are PML elements.'
    write(IMAIN,*) 'Average duration of the calculation per spectral element = ', &
                         duration_of_time_loop_in_seconds / real(NSTEP * NSPEC),' s'
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
  use specfem_par_gpu

  implicit none

  ! local parameters
  integer :: i,j,ispec
  real(kind=CUSTOM_REAL) :: rhol,mul,kappal

  ! Fields transfer for imaging
  ! acoustic domains
  if (any_acoustic ) then
    call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic,potential_dot_acoustic, &
                                        potential_dot_dot_acoustic,Mesh_pointer)
  endif

  ! elastic domains
  if (any_elastic) then
    call transfer_fields_el_from_device(NDIM*NGLOB_AB,tmp_displ_2D,tmp_veloc_2D,tmp_accel_2D,Mesh_pointer)

    if (P_SV) then
      ! P-SV waves
      displ_elastic(1,:) = tmp_displ_2D(1,:)
      displ_elastic(2,:) = tmp_displ_2D(2,:)

      veloc_elastic(1,:) = tmp_veloc_2D(1,:)
      veloc_elastic(2,:) = tmp_veloc_2D(2,:)

      accel_elastic(1,:) = tmp_accel_2D(1,:)
      accel_elastic(2,:) = tmp_accel_2D(2,:)
    else
      ! SH waves
      displ_elastic(1,:) = tmp_displ_2D(1,:)
      veloc_elastic(1,:) = tmp_veloc_2D(1,:)
      accel_elastic(1,:) = tmp_accel_2D(1,:)
    endif
  endif

  ! finishes kernel calculations
  if (SIMULATION_TYPE == 3) then
    ! Kernel transfer
    ! acoustic domains
    if (any_acoustic) then
      call transfer_kernels_ac_to_host(Mesh_pointer,rho_ac_kl,kappa_ac_kl,NSPEC_AB)
      rho_ac_kl(:,:,:) = rho_ac_kl(:,:,:) * NSTEP_BETWEEN_COMPUTE_KERNELS
      kappa_ac_kl(:,:,:) = kappa_ac_kl(:,:,:) * NSTEP_BETWEEN_COMPUTE_KERNELS
      rhop_ac_kl(:,:,:) = rho_ac_kl(:,:,:) + kappa_ac_kl(:,:,:)
      alpha_ac_kl(:,:,:) = TWO *  kappa_ac_kl(:,:,:)
      print *
      print *,'Maximum value of rho prime kernel : ',maxval(rhop_ac_kl)
      print *,'Maximum value of alpha kernel : ',maxval(alpha_ac_kl)
      print *
    endif

    ! elastic domains
    if (any_elastic) then
      call transfer_kernels_el_to_host(Mesh_pointer,rho_kl,mu_kl,kappa_kl,NSPEC_AB)

      ! Multiply each kernel point with the local coefficient
      do ispec = 1, nspec
        if (ispec_is_elastic(ispec)) then
          do j = 1, NGLLZ
            do i = 1, NGLLX
              if (.not. assign_external_model) then
                mul = poroelastcoef(2,1,kmato(ispec))
                kappal = poroelastcoef(3,1,kmato(ispec)) - FOUR_THIRDS * mul
                rhol = density(1,kmato(ispec))
              else
                rhol = rhoext(i,j,ispec)
                mul = rhoext(i,j,ispec)*vsext(i,j,ispec)*vsext(i,j,ispec)
                kappal = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec) - FOUR_THIRDS * mul
              endif

              rho_kl(i,j,ispec) = - rhol * rho_kl(i,j,ispec)
              mu_kl(i,j,ispec) =  - TWO * mul * mu_kl(i,j,ispec)
              kappa_kl(i,j,ispec) = - kappal * kappa_kl(i,j,ispec)
            enddo
          enddo
        endif
      enddo
    endif
  endif

end subroutine it_transfer_from_GPU

!
!----------------------------------------------------------------------------------------
!

subroutine it_compute_integrated_energy_field_and_output()

  ! compute int_0^t v^2 dt and write it on file if needed

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,IIN,MAX_STRING_LEN,OUTPUT_FILES

  use specfem_par, only: myrank,it,coord,nspec,ibool,integrated_kinetic_energy_field,max_kinetic_energy_field, &
                         integrated_potential_energy_field,max_potential_energy_field,kinetic_effective_duration_field, &
                         potential_effective_duration_field,total_integrated_energy_field,max_total_energy_field, &
                         total_effective_duration_field,NSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP

  implicit none

  ! local variables
  integer :: ispec,iglob
!!!! DK DK commenting this out for now because "call execute_command_line" is Fortran 2008
!!!! DK DK and some older compilers do not support it yet. We can probably put it back in a few years.
  !integer :: statMkdir
  character(len=MAX_STRING_LEN)  :: filename

  !! ABAB Uncomment to write the velocity profile in acoustic part
  !real(kind=CUSTOM_REAL) :: cpl,kappal
  !double precision :: rhol
  !double precision :: lambdal_unrelaxed_elastic
  !! ABAB

  ! computes maximum energy and integrated energy fields
  call compute_energy_fields()

  ! Create directories
  if (it == 1) then
    if (myrank == 0) then
!!!! DK DK commenting this out for now because "call execute_command_line" is Fortran 2008
!!!! DK DK and some older compilers do not support it yet. We can probably put it back in a few years.
     !call execute_command_line('mkdir -p '//trim(OUTPUT_FILES)//'energyFields',wait = .true.,cmdstat = statMkdir)
     !if (statMkdir /= 0) call exit_MPI(myrank,'Impossible to create '//trim(OUTPUT_FILES)//'energyFields')

     !call execute_command_line('mkdir -p '//trim(OUTPUT_FILES)//'energyFields/kinetic',wait = .true.,cmdstat = statMkdir)
     !if (statMkdir /= 0) call exit_MPI(myrank,'Impossible to create '//trim(OUTPUT_FILES)//'energyFields/kinetic')

     !call execute_command_line('mkdir -p '//trim(OUTPUT_FILES)//'energyFields/potential',wait = .true.,cmdstat = statMkdir)
     !if (statMkdir /= 0) call exit_MPI(myrank,'Impossible to create '//trim(OUTPUT_FILES)//'energyFields/potential')

     !call execute_command_line('mkdir -p '//trim(OUTPUT_FILES)//'energyFields/total',wait = .true.,cmdstat = statMkdir)
     !if (statMkdir /= 0) call exit_MPI(myrank,'Impossible to create '//trim(OUTPUT_FILES)//'energyFields/total')
    endif
    !call synchronize_all() ! Wait for first proc to create directories
  endif

  if (mod(it,NSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then
    ! write integrated kinetic energy field in external file
    write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'energyFields/kinetic/integrated_kinetic_energy_field',myrank
    open(unit=IIN,file=trim(filename),status='unknown',action='write')
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(integrated_kinetic_energy_field(ispec),4)
    enddo
    close(IIN)

    ! write max kinetic energy field in external file
    write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'energyFields/kinetic/max_kinetic_energy_field',myrank
    open(unit=IIN,file=trim(filename),status='unknown',action='write')
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(max_kinetic_energy_field(ispec),4)
    enddo
    close(IIN)

    ! write integrated potential energy field in external file
    write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'energyFields/potential/integrated_potential_energy_field',myrank
    open(unit=IIN,file=trim(filename),status='unknown',action='write')
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(integrated_potential_energy_field(ispec),4)
    enddo
    close(IIN)

    ! write max potential energy field in external file
    write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'energyFields/potential/max_potential_energy_field',myrank
    open(unit=IIN,file=trim(filename),status='unknown',action='write')
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(max_potential_energy_field(ispec),4)
    enddo
    close(IIN)

    ! write potential effective duration field in external file
    write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'energyFields/potential/potential_effective_duration_field',myrank
    open(unit=IIN,file=trim(filename),status='unknown',action='write')
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(potential_effective_duration_field(ispec),4)
    enddo
    close(IIN)

   ! write kinetic effective duration field in external file
    write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'energyFields/kinetic/kinetic_effective_duration_field',myrank
    open(unit=IIN,file=trim(filename),status='unknown',action='write')
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(kinetic_effective_duration_field(ispec),4)
    enddo
    close(IIN)

    ! write total integrated energy field in external file
    write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'energyFields/total/total_integrated_energy_field',myrank
    open(unit=IIN,file=trim(filename),status='unknown',action='write')
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(total_integrated_energy_field(ispec),4)
    enddo
    close(IIN)

    ! write max total energy field in external file
    write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'energyFields/total/max_total_energy_field',myrank
    open(unit=IIN,file=trim(filename),status='unknown',action='write')
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(max_total_energy_field(ispec),4)
    enddo
    close(IIN)

    ! write total effective duration field in external file
    write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'energyFields/total/total_effective_duration_field',myrank
    open(unit=IIN,file=trim(filename),status='unknown',action='write')
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(total_effective_duration_field(ispec),4)
    enddo
    close(IIN)
  endif

  ! ABAB Uncomment to write the velocity profile in the acoustic part in file
  !
  !  write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'velocities',myrank
  !  open(unit=IIN,file=trim(filename),status='unknown',action='write')
  !
  !  if (mod(it,NSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then
  !    ! loop over spectral elements
  !    do ispec = 1,nspec
  !      if (ispec_is_acoustic(ispec)) then
  !        ! get density of current spectral element
  !        lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
  !        rhol  = density(1,kmato(ispec))
  !        kappal = lambdal_unrelaxed_elastic
  !        cpl = sqrt(kappal/rhol)
  !
  !        !--- if external medium, get density of current grid point
  !        if (assign_external_model) then
  !          cpl = vpext(2,2,ispec)
  !        endif
  !        iglob = ibool(2,2,ispec)
  !        write(IIN,*) real(coord(2,iglob),4),cpl
  !      endif
  !    enddo
  !  endif
  !  close(IIN)

end subroutine it_compute_integrated_energy_field_and_output

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

    if (mod(NSTEP-it+1,NSTEP_BETWEEN_COMPUTE_KERNELS) == 0 .or. it == NSTEP ) then

      call save_forward_arrays_no_backward()

    endif

  endif ! SAVE_FORWARD

  if (SIMULATION_TYPE == 3) then

    if (mod(it,NSTEP_BETWEEN_COMPUTE_KERNELS) == 0 .or. it == 1 .or. it == NSTEP) then

      call read_forward_arrays_no_backward()

      ! For the first iteration we need to add a transfer, to initiate the
      ! transfer on the GPU
      if (no_backward_iframe == 1) call read_forward_arrays_no_backward()

      ! If the wavefield is requested at it=1, then we enforce another transfer
      if (no_backward_iframe == 2 .and. NSTEP_BETWEEN_COMPUTE_KERNELS == 1) &
                                   call read_forward_arrays_no_backward()

    endif

  endif ! SIMULATION == 3

end subroutine  manage_no_backward_reconstruction_io
