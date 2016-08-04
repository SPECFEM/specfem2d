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

subroutine iterate_time()

#ifdef USE_MPI
  use mpi
#endif

  use constants,only: IMAIN,NOISE_SAVE_EVERYWHERE
  use specfem_par
  use specfem_par_gpu
  use specfem_par_noise,only: NOISE_TOMOGRAPHY

  implicit none

  ! local parameters
  ! time
  character(len=8) :: datein
  character(len=10) :: timein
  character(len=5) :: zone
  integer, dimension(8) :: time_values
  integer :: year,mon,day,hr,minutes,timestamp

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

! *********************************************************
! ************* MAIN LOOP OVER THE TIME STEPS *************
! *********************************************************

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
      if (SIMULATION_TYPE == 3) call update_displ_Newmark_backward()

      ! acoustic domains
      if (ACOUSTIC_SIMULATION) then
        if (.not. GPU_MODE) then
          call compute_forces_acoustic_main()
          if (SIMULATION_TYPE == 3) call compute_forces_acoustic_main_backward()
        else
          ! on GPU
          if (any_acoustic) call compute_forces_acoustic_GPU()
        endif
      endif

      ! gravitoacoustic domains
      if (GRAVITOACOUSTIC_SIMULATION) then
        if (.not. GPU_MODE) then
          call compute_forces_gravitoacoustic_main()
        else
          ! on GPU
          if (any_gravitoacoustic) call exit_MPI(myrank,'gravitoacoustic not implemented in GPU MODE yet')
        endif
      endif

      ! elastic domains
      if (ELASTIC_SIMULATION) then
        if (.not. GPU_MODE) then
          call compute_forces_viscoelastic_main()
          if (SIMULATION_TYPE == 3) call compute_forces_viscoelastic_main_backward()
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
    if (SIMULATION_TYPE == 3 .and. it == 1) then
      call read_forward_arrays()
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

    ! computes kinetic and potential energy
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
    if (SIMULATION_TYPE == 3) then
      call compute_kernels()
    endif

    ! display results at given time steps
    call write_movie_output()

  enddo ! end of the main time loop

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

  use constants,only: TWO,FOUR_THIRDS
  use specfem_par
  use specfem_par_gpu

  implicit none

  ! local parameters
  integer :: i,j,ispec,iglob
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

      rhop_ac_kl(:,:,:) = rho_ac_kl(:,:,:) + kappa_ac_kl(:,:,:)
      alpha_ac_kl(:,:,:) = TWO *  kappa_ac_kl(:,:,:)
    endif

    ! elastic domains
    if (any_elastic) then
      call transfer_kernels_el_to_host(Mesh_pointer,rho_kl,mu_kl,kappa_kl,NSPEC_AB)

      ! Multiply each kernel point with the local coefficient
      do ispec = 1, nspec
        if (ispec_is_elastic(ispec)) then
          do j = 1, NGLLZ
            do i = 1, NGLLX
              iglob = ibool(i,j,ispec)
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

  use constants,only:CUSTOM_REAL,NGLLX,NGLLZ,IIN,MAX_STRING_LEN

  use specfem_par,only: myrank,it,coord,nspec,ibool,integrated_kinetic_energy_field,max_kinetic_energy_field, &
                        integrated_potential_energy_field,max_potential_energy_field,kinetic_effective_duration_field, &
                        potential_effective_duration_field,total_integrated_energy_field,max_total_energy_field, &
                        total_effective_duration_field,NSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP
                        !poroelastcoef,kmato,density,assign_external_model,vpext,ispec_is_acoustic ! ABAB

  implicit none

  ! local variables
  integer :: ispec,iglob,statMkdir,system !i,j
  character(len=MAX_STRING_LEN)  :: filename

  !! ABAB Uncomment to write the velocity profile in acoustic part
  !real(kind=CUSTOM_REAL) :: cpl,kappal
  !double precision :: rhol
  !double precision :: lambdal_unrelaxed_elastic
  !! ABAB

  ! computes maximum energy and integrated energy fields
  call compute_energy_fields()

  ! Create directories
  ! TODO replace by execute_command_line or makedirqq (change compilation)
  if (it == 1) then
    if (myrank == 0) then
      statMkdir = system('mkdir -p ./OUTPUT_FILES/energyFields')
      if(statMkdir .ne. 0) then
        call exit_MPI(myrank,'Impossible to create ./OUTPUT_FILES/energyFields')
      endif
     statMkdir = system('mkdir -p ./OUTPUT_FILES/energyFields/kinetic')
     if(statMkdir .ne. 0) then
       call exit_MPI(myrank,'Impossible to create ./OUTPUT_FILES/energyFields/kinetic')
     endif
     statMkdir = system('mkdir -p ./OUTPUT_FILES/energyFields/potential')
     if(statMkdir .ne. 0) then
       call exit_MPI(myrank,'Impossible to create ./OUTPUT_FILES/energyFields/potential')
     endif
     statMkdir = system('mkdir -p ./OUTPUT_FILES/energyFields/total')
     if(statMkdir .ne. 0) then
       call exit_MPI(myrank,'Impossible to create ./OUTPUT_FILES/energyFields/total')
     endif
     !if(.not. makedirqq("./OUTPUT_FILES/energyFields/total")) then ! ABAB This would be better
     !  call exit_MPI(myrank,'Impossible to create ./OUTPUT_FILES/energyFields/total')
     !endif
    endif
  endif

  call synchronize_all() ! Wait for first proc to create directories

  ! write integrated kinetic energy field in external file

  write(filename,"('./OUTPUT_FILES/energyFields/kinetic/integrated_kinetic_energy_field',i5.5)") myrank
  open(unit=IIN,file=trim(filename),status='unknown',action='write')

  if (mod(it,NSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(integrated_kinetic_energy_field(ispec),4)
    enddo
  endif
  close(IIN)

  ! write max kinetic energy field in external file

  write(filename,"('./OUTPUT_FILES/energyFields/kinetic/max_kinetic_energy_field',i5.5)") myrank
  open(unit=IIN,file=trim(filename),status='unknown',action='write')

  if (mod(it,NSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(max_kinetic_energy_field(ispec),4)
    enddo
  endif
  close(IIN)

  ! write integrated potential energy field in external file

  write(filename,"('./OUTPUT_FILES/energyFields/potential/integrated_potential_energy_field',i5.5)") myrank
  open(unit=IIN,file=trim(filename),status='unknown',action='write')

  if (mod(it,NSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(integrated_potential_energy_field(ispec),4)
    enddo
  endif
  close(IIN)

  ! write max potential energy field in external file

  write(filename,"('./OUTPUT_FILES/energyFields/potential/max_potential_energy_field',i5.5)") myrank
  open(unit=IIN,file=trim(filename),status='unknown',action='write')

  if (mod(it,NSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(max_potential_energy_field(ispec),4)
    enddo
  endif
  close(IIN)

  ! write potential effective duration field in external file

  write(filename,"('./OUTPUT_FILES/energyFields/potential/potential_effective_duration_field',i5.5)") myrank
  open(unit=IIN,file=trim(filename),status='unknown',action='write')

  if (mod(it,NSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(potential_effective_duration_field(ispec),4)
    enddo
  endif
  close(IIN)

  ! write kinetic effective duration field in external file

  write(filename,"('./OUTPUT_FILES/energyFields/kinetic/kinetic_effective_duration_field',i5.5)") myrank
  open(unit=IIN,file=trim(filename),status='unknown',action='write')

  if (mod(it,NSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(kinetic_effective_duration_field(ispec),4)
    enddo
  endif
  close(IIN)

  ! write total integrated energy field in external file

  write(filename,"('./OUTPUT_FILES/energyFields/total/total_integrated_energy_field',i5.5)") myrank
  open(unit=IIN,file=trim(filename),status='unknown',action='write')

  if (mod(it,NSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(total_integrated_energy_field(ispec),4)
    enddo
  endif
  close(IIN)

  ! write max total energy field in external file

  write(filename,"('./OUTPUT_FILES/energyFields/total/max_total_energy_field',i5.5)") myrank
  open(unit=IIN,file=trim(filename),status='unknown',action='write')

  if (mod(it,NSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(max_total_energy_field(ispec),4)
    enddo
  endif
  close(IIN)

  ! write total effective duration field in external file

  write(filename,"('./OUTPUT_FILES/energyFields/total/total_effective_duration_field',i5.5)") myrank
  open(unit=IIN,file=trim(filename),status='unknown',action='write')

  if (mod(it,NSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(total_effective_duration_field(ispec),4)
    enddo
  endif
  close(IIN)

  ! ABAB Uncomment to write the velocity profile in the acoustic part in file
  !
  !  write(filename,"('./OUTPUT_FILES/velocities',i5.5)") myrank
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

