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

  use constants, only: IMAIN
  use specfem_par
  use specfem_par_noise, only: NOISE_TOMOGRAPHY

  implicit none

#ifdef USE_MPI
  include "precision.h"
#endif

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_potential_acoustic_buffer,b_potential_dot_dot_acoustic_buffer
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_displ_elastic_buffer,b_accel_elastic_buffer
  double precision :: sizeval

  integer :: i,j
  integer :: it_temp,seismo_current_temp
  integer :: ier
  integer, parameter :: it_begin = 1

  ! time
  character(len=8) :: datein
  character(len=10) :: timein
  character(len=5) :: zone
  integer, dimension(8) :: time_values
  integer :: year,mon,day,hr,minutes,timestamp

  ! checks if anything to do
  if (.not. UNDO_ATTENUATION) return

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
  if (GPU_MODE ) call exit_MPI(myrank,'for undo_attenuation, GPU_MODE is not supported')
  if (time_stepping_scheme /= 1 ) call exit_MPI(myrank,'for undo_attenuation, only Newmark scheme has implemented ')
  if (any_gravitoacoustic ) call exit_MPI(myrank,'undo_attenuation has not implemented for gravitoacoustic yet')
  if (any_poroelastic ) call exit_MPI(myrank,'undo_attenuation has not implemented for poroelastic simulation yet')
  if (NOISE_TOMOGRAPHY /= 0 ) call exit_MPI(myrank,'for undo_attenuation, NOISE_TOMOGRAPHY is not supported')
  if (AXISYM ) call exit_MPI(myrank,'Just axisymmetric FORWARD simulations are possible so far')

  ! number of time subsets for time loop
  NSUBSET_ITERATIONS = ceiling( dble(NSTEP)/dble(NT_DUMP_ATTENUATION) )

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
        sizeval = 2 * dble(nglob) * dble(NT_DUMP_ATTENUATION) * dble(CUSTOM_REAL) / 1024.d0 / 1024.d0
        write(IMAIN,*) '  size of acoustic wavefield buffer per slice      = ', sngl(sizeval),'MB'
      endif
      if (any_elastic) then
        ! buffer(3,nglob,NT_DUMP_ATTENUATION) in MB
        sizeval = 2 * dble(NDIM) * dble(nglob) * dble(NT_DUMP_ATTENUATION) * dble(CUSTOM_REAL) / 1024.d0 / 1024.d0
        write(IMAIN,*) '  size of elastic wavefield buffer per slice       = ', sngl(sizeval),'MB'
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! allocates buffers
  if (SIMULATION_TYPE == 3) then
    if (any_acoustic) then
      allocate(b_potential_acoustic_buffer(nglob,NT_DUMP_ATTENUATION),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'error allocating b_potential_acoustic')
      allocate(b_potential_dot_dot_acoustic_buffer(nglob,NT_DUMP_ATTENUATION),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'error allocating b_potential_dot_dot_acoustic')
    endif

    if (any_elastic) then
      allocate(b_displ_elastic_buffer(NDIM,nglob,NT_DUMP_ATTENUATION),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'error allocating b_displ_elastic')
      allocate(b_accel_elastic_buffer(NDIM,nglob,NT_DUMP_ATTENUATION),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'error allocating b_accel_elastic')

      allocate(b_e1(NGLLX,NGLLZ,nspec_ATT,N_SLS), &
               b_e11(NGLLX,NGLLZ,nspec_ATT,N_SLS), &
               b_e13(NGLLX,NGLLZ,nspec_ATT,N_SLS),stat=ier)
      if (ier /= 0) stop 'Error allocating attenuation arrays'
      b_e1(:,:,:,:) = 0._CUSTOM_REAL
      b_e11(:,:,:,:) = 0._CUSTOM_REAL
      b_e13(:,:,:,:) = 0._CUSTOM_REAL
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

  ! initializes time increments
  it = 0

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
      do it_of_this_subset = 1, it_subset_end

        it = it + 1

        ! compute current time
        timeval = (it-1) * deltat

        ! display time step and max of norm of displacement
        if (mod(it,NSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then
          call check_stability()
        endif

        do i_stage = 1, stage_time_scheme ! is equal to 1 if Newmark because only one stage then
          ! update_displacement_newmark
          call update_displ_Newmark()

          ! main solver for the acoustic elements
          call compute_forces_acoustic_main()

          ! main solver for the elastic elements
          call compute_forces_viscoelastic_main()
        enddo

        ! computes kinetic and potential energy
        if (output_energy) then
          call compute_and_output_energy()
        endif

        ! loop on all the receivers to compute and store the seismograms
        call write_seismograms()

        ! display results at given time steps
        call write_movie_output()

      enddo ! subset loop

    case (3)
      ! kernel simulations

      ! intermediate storage of it and seismo_current positions
      it_temp = it
      seismo_current_temp = seismo_current

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

      ! subset loop
      do it_of_this_subset = 1, it_subset_end

        it = it + 1
        do i_stage = 1, stage_time_scheme
          ! backward_inner_loop
          ! update_displacement_newmark
          call update_displ_acoustic_backward()
          call update_displ_elastic_backward()

          ! main solver for the acoustic elements
          call compute_forces_acoustic_main_backward()

          ! main solver for the elastic elements
          call compute_forces_viscoelastic_main_backward()
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

      ! resets current it and seismo_current positions
      it = it_temp
      seismo_current = seismo_current_temp

      ! adjoint wavefield simulation
      do it_of_this_subset = 1, it_subset_end
        ! reads backward/reconstructed wavefield from buffers
        ! note: uses wavefield at corresponding time (NSTEP - it + 1 ), i.e. we have now time-reversed wavefields
        if (any_acoustic) then
          do j = 1,NGLOB
            b_potential_acoustic(j) = b_potential_acoustic_buffer(j,it_subset_end-it_of_this_subset+1)
            b_potential_dot_dot_acoustic(j) = b_potential_dot_dot_acoustic_buffer(j,it_subset_end-it_of_this_subset+1)
          enddo
        endif

        if (any_elastic) then
          do i = 1,NDIM
            do j = 1,NGLOB
              b_displ_elastic(i,j) = b_displ_elastic_buffer(i,j,it_subset_end-it_of_this_subset+1)
              b_accel_elastic(i,j) = b_accel_elastic_buffer(i,j,it_subset_end-it_of_this_subset+1)
            enddo
          enddo
        endif

        ! safety stop
        if (NOISE_TOMOGRAPHY == 3) stop 'Undo_attenuation for noise kernels not implemented yet'

        it = it + 1

        ! compute current time
        timeval = (it-1) * deltat

        ! display time step and max of norm of displacement
        if (mod(it,NSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then
          call check_stability()
        endif

        do i_stage = 1, stage_time_scheme
          ! adjoint
          ! update_displacement_newmark
          call update_displ_acoustic_forward()

          if (any_acoustic) then
            !ZN here we remove the trick introduced by Luoyang to stabilized the adjoint simulation
            !ZN However in order to keep the current code be consistent, we still keep potential_acoustic_adj_coupling
            !ZN the final goal should remove the *adj_coupling
            potential_acoustic_adj_coupling(:) = potential_acoustic(:)
          endif

          call update_displ_elastic_forward()

!daniel TODO: not sure if the following below is correct or needs to switch orders
!             usually one computes first the updated pressure and afterwards computes the elastic domain
!             and its coupling terms...please make sure...

          ! main solver for the elastic elements
          call compute_forces_viscoelastic_main()

          ! for coupling with adjoint wavefields, stores temporary old wavefield
          if (coupled_acoustic_elastic) then
            ! handles adjoint runs coupling between adjoint potential and adjoint elastic wavefield
            ! adjoint definition: \partial_t^2 \bfs^\dagger = - \frac{1}{\rho} \bfnabla \phi^\dagger
#ifdef FORCE_VECTORIZATION
            do i = 1,NDIM * nglob_elastic
              accel_elastic_adj_coupling(i,1) = - accel_elastic(i,1)
            enddo
#else
            accel_elastic_adj_coupling(:,:) = - accel_elastic(:,:)
#endif
          endif

          ! main solver for the acoustic elements
          call compute_forces_acoustic_main()
        enddo

        ! computes kinetic and potential energy
        if (output_energy) then
          call compute_and_output_energy()
        endif

        ! loop on all the receivers to compute and store the seismograms
        call write_seismograms()

        ! kernels calculation
        call compute_kernels()

        ! display results at given time steps
        call write_movie_output()

      enddo ! subset loop
    end select ! SIMULATION_TYPE

  enddo   ! end of main time loop

  !
  !---- end of time iteration loop
  !

  ! frees undo_attenuation buffers
  if (SIMULATION_TYPE == 3) then
    if (any_acoustic) then
      deallocate(b_potential_acoustic_buffer,b_potential_dot_dot_acoustic_buffer)
    endif
    if (any_elastic) then
      deallocate(b_displ_elastic_buffer,b_accel_elastic_buffer)
      deallocate(b_e1,b_e11,b_e13)
    endif
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !----  formats
  400 format(/1x,41('=')/,' =  T i m e  e v o l u t i o n  l o o p  ='/1x,41('=')/)

  end subroutine iterate_time_undoatt
