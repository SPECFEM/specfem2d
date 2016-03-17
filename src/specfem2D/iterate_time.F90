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

  use specfem_par
  use specfem_par_gpu
  use specfem_par_noise,only: NOISE_TOMOGRAPHY,save_everywhere

  implicit none

  ! local parameters
  integer :: ier

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
      call update_displacement_scheme()

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
        else
          ! on GPU
          if (any_poroelastic) call exit_MPI(myrank,'poroelastic not implemented in GPU MODE yet')
        endif
      endif

    enddo ! stage_time_scheme (LDDRK or RK)

    ! reads in lastframe for adjoint/kernels calculation
    if (SIMULATION_TYPE == 3 .and. it == 1) then
      call it_read_forward_arrays()
    endif

    ! noise simulations
    if (NOISE_TOMOGRAPHY == 1) then
      call save_surface_movie_noise()
    else if (NOISE_TOMOGRAPHY == 2 .and. save_everywhere) then
      call save_surface_movie_noise()
    else if (NOISE_TOMOGRAPHY == 3 .and. save_everywhere) then
      if (it == 1) open(unit=500,file='OUTPUT_FILES/NOISE_TOMOGRAPHY/phi',access='direct', &
                        recl=nglob*CUSTOM_REAL,action='write',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error retrieving noise ensemble forward wavefield')

      ! safety check
      if (P_SV) then
        ! P-SV case
        call exit_MPI(myrank,'P-SV case not yet implemented.')
      else
        ! SH case
        read(unit=500,rec=NSTEP-it+1) b_displ_elastic(2,:)
      endif
    endif

    ! computes kinetic and potential energy
    if (output_energy) then
      call it_compute_and_output_energy()
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

subroutine it_read_forward_arrays()

! restores last time snapshot saved for backward/reconstruction of wavefields
! note: this is done here after the Newmark time scheme, otherwise the indexing for sources
!          and adjoint sources will become more complicated
!          that is, index it for adjoint sources will match index NSTEP - 1 for backward/reconstructed wavefields

  use specfem_par
  use specfem_par_gpu

  implicit none

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: outputname

  ! acoustic medium
  if (any_acoustic) then
    write(outputname,'(a,i6.6,a)') 'lastframe_acoustic',myrank,'.bin'
    open(unit=55,file='OUTPUT_FILES/'//trim(outputname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file OUTPUT_FILES/lastframe_acoustic**.bin')

    read(55) b_potential_acoustic
    read(55) b_potential_dot_acoustic
    read(55) b_potential_dot_dot_acoustic

    close(55)

    if (GPU_MODE) then
      ! transfers fields onto GPU
      call transfer_b_fields_ac_to_device(NGLOB_AB,b_potential_acoustic, &
                                          b_potential_dot_acoustic,      &
                                          b_potential_dot_dot_acoustic,  &
                                          Mesh_pointer)
    else
      ! free surface for an acoustic medium
      call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic)
    endif
  endif

  ! elastic medium
  if (any_elastic) then
    write(outputname,'(a,i6.6,a)') 'lastframe_elastic',myrank,'.bin'
    open(unit=55,file='OUTPUT_FILES/'//trim(outputname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file OUTPUT_FILES/lastframe_elastic**.bin')

    read(55) b_displ_elastic
    read(55) b_veloc_elastic
    read(55) b_accel_elastic
    close(55)

    !SH (membrane) waves
    if (.not. P_SV) then
      ! only index array(1,:) contains SH wavefield
      b_displ_elastic(2,:) = 0._CUSTOM_REAL
      b_veloc_elastic(2,:) = 0._CUSTOM_REAL
      b_accel_elastic(2,:) = 0._CUSTOM_REAL
    endif

    if (GPU_MODE) then
      ! prepares wavefields for transfering
      if (P_SV) then
        tmp_displ_2D(1,:) = b_displ_elastic(1,:)
        tmp_displ_2D(2,:) = b_displ_elastic(2,:)
        tmp_veloc_2D(1,:) = b_veloc_elastic(1,:)
        tmp_veloc_2D(2,:) = b_veloc_elastic(2,:)
        tmp_accel_2D(1,:) = b_accel_elastic(1,:)
        tmp_accel_2D(2,:) = b_accel_elastic(2,:)
      else
        ! SH waves
        tmp_displ_2D(1,:) = b_displ_elastic(1,:)
        tmp_displ_2D(2,:) = 0._CUSTOM_REAL
        tmp_veloc_2D(1,:) = b_veloc_elastic(1,:)
        tmp_veloc_2D(2,:) = 0._CUSTOM_REAL
        tmp_accel_2D(1,:) = b_accel_elastic(1,:)
        tmp_accel_2D(2,:) = 0._CUSTOM_REAL
      endif
      call transfer_b_fields_to_device(NDIM*NGLOB_AB,tmp_displ_2D,tmp_veloc_2D,tmp_accel_2D,Mesh_pointer)
    endif
  endif

  ! poroelastic medium
  if (any_poroelastic) then
    write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_s',myrank,'.bin'
    open(unit=55,file='OUTPUT_FILES/'//trim(outputname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file OUTPUT_FILES/lastframe_poroelastic_s**.bin')

    read(55) b_displs_poroelastic
    read(55) b_velocs_poroelastic
    read(55) b_accels_poroelastic
    close(55)

    write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_w',myrank,'.bin'
    open(unit=56,file='OUTPUT_FILES/'//trim(outputname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file OUTPUT_FILES/lastframe_poroelastic_w**.bin')

    read(56) b_displw_poroelastic
    read(56) b_velocw_poroelastic
    read(56) b_accelw_poroelastic
    close(56)

    ! safety check
    if (GPU_MODE) then
      stop 'GPU_MODE error: sorry, reading lastframe from poroelastic simulation not implemented yet'
    endif
  endif

end subroutine it_read_forward_arrays

!
!----------------------------------------------------------------------------------------
!

subroutine it_compute_and_output_energy()

  use constants,only: IOUT_ENERGY,CUSTOM_REAL

  use specfem_par,only: GPU_MODE,myrank,it,deltat,kinetic_energy,potential_energy,t0

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: kinetic_energy_total,potential_energy_total

  ! safety check
  if (GPU_MODE) stop 'Error computing energy for output is not implemented on GPUs yet'

  ! computes energy
  call compute_energy()

  ! computes total for all processes
  call sum_all_cr(kinetic_energy,kinetic_energy_total)
  call sum_all_cr(potential_energy,potential_energy_total)

  ! saves kinetic, potential and total energy for this time step in external file
  if (myrank == 0) then
    write(IOUT_ENERGY,*) real(dble(it-1)*deltat - t0,4),real(kinetic_energy_total,4), &
                         real(potential_energy_total,4),real(kinetic_energy_total + potential_energy_total,4)
  endif

end subroutine it_compute_and_output_energy

!
!----------------------------------------------------------------------------------------
!

subroutine it_compute_integrated_energy_field_and_output()
  ! compute int_0^t v^2 dt and write it on file if needed

  use constants,only:CUSTOM_REAL,NGLLX,NGLLZ,IIN

  use specfem_par,only: myrank,it,coord,nspec,ibool,integrated_energy_field,max_energy_field, &
                        NSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP

  implicit none

  ! local variables
  integer :: ispec,iglob!,i,j
  character(len=256)  :: filename

  ! computes maximum energy and integrated energy field
  call compute_energy_fields()

  ! write integrated energy field in external file

  write(filename,"('./OUTPUT_FILES/integrated_energy_field',i5.5)") myrank
  open(unit=IIN,file=trim(filename),status='unknown',action='write')

  if (mod(it,NSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(integrated_energy_field(ispec),4)
    enddo
  endif
  close(IIN)

  ! write max energy field in external file

  write(filename,"('./OUTPUT_FILES/max_energy_field',i5.5)") myrank
  open(unit=IIN,file=trim(filename),status='unknown',action='write')

  if (mod(it,NSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(max_energy_field(ispec),4)
    enddo
  endif
  close(IIN)

end subroutine it_compute_integrated_energy_field_and_output

