
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
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and Inria at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

subroutine iterate_time()

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par
  implicit none

  integer i,j,ispec,iglob,it_temp

#ifdef USE_MPI
  include "precision.h"
#endif


  if (myrank == 0) write(IOUT,400) ! Write = T i m e  e v o l u t i o n  l o o p =
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

  do it = 1,NSTEP
    ! compute current time
    timeval = (it-1)*deltat

    do i_stage=1, stage_time_scheme

      if( GPU_MODE ) then
        call update_displacement_precondition_newmark_GPU()
      endif

      if( .not. GPU_MODE ) then

        if( any_acoustic ) then
          ! free surface for an acoustic medium
          if( nelem_acoustic_surface > 0 ) then
            call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                               potential_acoustic)
          endif
          if( time_stepping_scheme == 1 ) then
            call update_displacement_precondition_newmark_acoustic(deltat,deltatover2,deltatsquareover2,&
                                                                   potential_dot_dot_acoustic,potential_dot_acoustic,&
                                                                   potential_acoustic,potential_acoustic_old, &
                                                                   PML_BOUNDARY_CONDITIONS)
          else
#ifdef FORCE_VECTORIZATION
            do i = 1,nglob_acoustic
              potential_dot_dot_acoustic(i) = 0._CUSTOM_REAL
            enddo
#else
            potential_dot_dot_acoustic = 0._CUSTOM_REAL
#endif
          endif

          if( SIMULATION_TYPE == 3 ) then
            !Since we do not do anything in PML region in case of backward simulation, thus we set
            !PML_BOUNDARY_CONDITIONS = .false.
            if( time_stepping_scheme == 1 ) then
              call update_displacement_precondition_newmark_acoustic(b_deltat,b_deltatover2,b_deltatsquareover2,&
                                                                     b_potential_dot_dot_acoustic,b_potential_dot_acoustic,&
                                                                     b_potential_acoustic,b_potential_acoustic_old, &
                                                                     .false.)
            else
#ifdef FORCE_VECTORIZATION
              do i = 1,nglob_acoustic
                b_potential_dot_dot_acoustic(i) = 0._CUSTOM_REAL
              enddo
#else
              b_potential_dot_dot_acoustic = 0._CUSTOM_REAL
#endif
            endif
          endif
        endif

        if( any_elastic ) then
          if( time_stepping_scheme == 1 ) then
            if( SIMULATION_TYPE == 3 ) then
#ifdef FORCE_VECTORIZATION
              do i = 1,3*nglob_elastic
                accel_elastic_adj_coupling(i,1) = accel_elastic(i,1)
              enddo
#else
              accel_elastic_adj_coupling = accel_elastic
#endif
            endif

            if( time_stepping_scheme == 1 ) then
              call update_displacement_precondition_newmark_elastic(deltat,deltatover2,deltatsquareover2,&
                                                                    accel_elastic,veloc_elastic,&
                                                                    displ_elastic,displ_elastic_old,&
                                                                    PML_BOUNDARY_CONDITIONS)
            else
#ifdef FORCE_VECTORIZATION
              do i = 1,3*nglob_elastic
                accel_elastic(i,1) = 0._CUSTOM_REAL
              enddo
#else
              accel_elastic = 0._CUSTOM_REAL
#endif
            endif
          endif

          if( SIMULATION_TYPE == 3 ) then
            !Since we do not do anything in PML region in case of backward simulation, thus we set
            !PML_BOUNDARY_CONDITIONS = .false.
            if( time_stepping_scheme == 1 ) then
              call update_displacement_precondition_newmark_elastic(b_deltat,b_deltatover2,b_deltatsquareover2,&
                                                                    b_accel_elastic,b_veloc_elastic,&
                                                                    b_displ_elastic,b_displ_elastic_old,&
                                                                    .false.)
            else
#ifdef FORCE_VECTORIZATION
              do i = 1,3*nglob_elastic
                b_accel_elastic(i,1) = 0._CUSTOM_REAL
              enddo
#else
              b_accel_elastic = 0._CUSTOM_REAL
#endif
            endif
          endif
        endif

        if( AXISYM ) then
          call enforce_zero_radial_displacements_on_the_axis()
        endif

        if( any_poroelastic ) then
          if( SIMULATION_TYPE == 3 ) then
            accels_poroelastic_adj_coupling = accels_poroelastic
            accelw_poroelastic_adj_coupling = accelw_poroelastic
          endif

          if( time_stepping_scheme == 1 ) then
            call update_displacement_precondition_newmark_poroelastic(deltat,deltatover2,deltatsquareover2,&
                                                                      accels_poroelastic,velocs_poroelastic,&
                                                                      displs_poroelastic,accelw_poroelastic,&
                                                                      velocw_poroelastic,displw_poroelastic)
          else
            accels_poroelastic = 0._CUSTOM_REAL
            accelw_poroelastic = 0._CUSTOM_REAL
          endif

          if( SIMULATION_TYPE == 3 ) then
            if( time_stepping_scheme == 1 ) then
              !PML did not implemented for poroelastic simulation
              call update_displacement_precondition_newmark_poroelastic(b_deltat,b_deltatover2,b_deltatsquareover2,&
                                                                        b_accels_poroelastic,b_velocs_poroelastic,&
                                                                        b_displs_poroelastic,b_accelw_poroelastic,&
                                                                        b_velocw_poroelastic,b_displw_poroelastic)
            else
              b_accels_poroelastic = 0._CUSTOM_REAL
              b_accelw_poroelastic = 0._CUSTOM_REAL
            endif
          endif
        endif
! *********************************************************
! ************* main solver for the acoustic elements
! *********************************************************
        if( any_acoustic ) then
          call compute_forces_acoustic(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                       potential_acoustic,potential_acoustic_old,PML_BOUNDARY_CONDITIONS)

          if( SIMULATION_TYPE == 3 ) then
            if( PML_BOUNDARY_CONDITIONS ) then
              it_temp = NSTEP-it+1
              call rebuild_value_on_PML_interface_acoustic(it_temp)
            endif

            call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                               b_potential_acoustic)
            call compute_forces_acoustic_backward(b_potential_dot_dot_acoustic,b_potential_acoustic)

            if( PML_BOUNDARY_CONDITIONS ) then
              it_temp = NSTEP-it+1
              call rebuild_value_on_PML_interface_acoustic(it_temp)
            endif
          endif

        endif ! end of test if any acoustic element

        ! *********************************************************
        ! ************* add acoustic forcing at a rigid boundary
        ! *********************************************************
        if( any_acoustic ) then
          if( ACOUSTIC_FORCING ) then
            call add_acoustic_forcing_at_rigid_boundary(potential_dot_dot_acoustic)
          endif
        endif ! end of test if any acoustic element

        ! *********************************************************
        ! ************* add coupling with the elastic side
        ! *********************************************************
        if( coupled_acoustic_elastic ) then
          if( SIMULATION_TYPE == 1 ) then
            call compute_coupling_acoustic_el(displ_elastic,displ_elastic_old,potential_dot_dot_acoustic, &
                                              PML_BOUNDARY_CONDITIONS)
          endif

          if( SIMULATION_TYPE == 3 ) then
            accel_elastic_adj_coupling2 = - accel_elastic_adj_coupling
            call compute_coupling_acoustic_el(accel_elastic_adj_coupling2,displ_elastic_old,potential_dot_dot_acoustic,&
                                              PML_BOUNDARY_CONDITIONS)

            call compute_coupling_acoustic_el_backward(b_displ_elastic,b_potential_dot_dot_acoustic)
          endif
        endif

        ! *********************************************************
        ! ************* add coupling with the poroelastic side
        ! *********************************************************
        if( coupled_acoustic_poro) then
          if( SIMULATION_TYPE == 1 ) then
            call compute_coupling_acoustic_po()
          endif

          if( SIMULATION_TYPE == 3 ) then
            call compute_coupling_acoustic_po()
            call compute_coupling_acoustic_po_backward()
          endif
        endif

        ! ************************************************************************************
        ! ************************************ add force source
        ! ************************************************************************************
        if( any_acoustic ) then
          if( .not. initialfield ) then
            if( SIMULATION_TYPE == 1 ) then
              call compute_add_sources_acoustic(potential_dot_dot_acoustic,it,i_stage)
            endif

            if( SIMULATION_TYPE == 3 ) then   ! adjoint and backward wavefield
              call compute_add_sources_acoustic_adjoint()
              call compute_add_sources_acoustic(b_potential_dot_dot_acoustic,NSTEP-it+1,stage_time_scheme-i_stage+1)
            endif ! SIMULATION_TYPE == 3 adjoint wavefield
          endif ! if not using an initial field
        endif !if(any_acoustic)

        ! ************************************************************************************
        ! ********** assembling potential_dot_dot or b_potential_dot_dot for acoustic elements
        ! ************************************************************************************
#ifdef USE_MPI
        if( nproc > 1 .and. any_acoustic .and. ninterface_acoustic > 0 ) then
          call assemble_MPI_vector_ac(potential_dot_dot_acoustic)

          if( time_stepping_scheme == 2 ) then
            if( i_stage==1 .and. it == 1 .and. (.not. initialfield) ) then
              potential_dot_acoustic_temp = potential_dot_acoustic
              call assemble_MPI_vector_ac(potential_dot_acoustic)
            endif
          endif

          if( SIMULATION_TYPE == 3) then
            call assemble_MPI_vector_ac(b_potential_dot_dot_acoustic)
          endif
        endif
#endif

        if( PML_BOUNDARY_CONDITIONS ) then
          if( any_acoustic .and. nglob_interface > 0 ) then
            if( SAVE_FORWARD .and. SIMULATION_TYPE == 1 ) then
              do i = 1, nglob_interface
                write(72)potential_dot_dot_acoustic(point_interface(i)),&
                         potential_dot_acoustic(point_interface(i)),&
                         potential_acoustic(point_interface(i))
              enddo
            endif

            if( SIMULATION_TYPE == 3 ) then
              do i = 1, nglob_interface
                b_potential_dot_dot_acoustic(point_interface(i)) = pml_interface_history_potential_dot_dot(i,NSTEP-it+1)
              enddo
            endif
          endif
        endif
! ************************************************************************************
! ************* multiply by the inverse of the mass matrix and update velocity
! ************************************************************************************

        if( any_acoustic ) then
          ! free surface for an acoustic medium
          if( nelem_acoustic_surface > 0 ) then
            call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                               potential_acoustic)
            if( SIMULATION_TYPE == 3 ) then
              call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                                 b_potential_acoustic)
            endif
          endif

          if( time_stepping_scheme == 1 ) then

            !! DK DK this should be vectorized
            potential_dot_dot_acoustic = potential_dot_dot_acoustic * rmass_inverse_acoustic
            potential_dot_acoustic = potential_dot_acoustic + deltatover2*potential_dot_dot_acoustic

            if( SIMULATION_TYPE == 3 ) then
            !! DK DK this should be vectorized
              b_potential_dot_dot_acoustic = b_potential_dot_dot_acoustic * rmass_inverse_acoustic
              b_potential_dot_acoustic = b_potential_dot_acoustic + b_deltatover2*b_potential_dot_dot_acoustic
            endif

            ! update the potential field (use a new array here) for coupling terms
            potential_acoustic_adj_coupling = potential_acoustic + deltat*potential_dot_acoustic + &
                                              deltatsquareover2*potential_dot_dot_acoustic
          endif

          if( time_stepping_scheme == 2 ) then
            !! DK DK this should be vectorized
            potential_dot_dot_acoustic = potential_dot_dot_acoustic * rmass_inverse_acoustic
            potential_dot_acoustic_LDDRK = alpha_LDDRK(i_stage) * potential_dot_acoustic_LDDRK + &
                                         deltat * potential_dot_dot_acoustic
            potential_acoustic_LDDRK = alpha_LDDRK(i_stage) * potential_acoustic_LDDRK + &
                                       deltat*potential_dot_acoustic

            if( i_stage==1 .and. it == 1 .and. (.not. initialfield) ) then
              !! DK DK this should be vectorized
              potential_dot_acoustic_temp = potential_dot_acoustic_temp + &
                                            beta_LDDRK(i_stage) * potential_dot_acoustic_LDDRK
              potential_dot_acoustic = potential_dot_acoustic_temp
            else
              potential_dot_acoustic = potential_dot_acoustic + beta_LDDRK(i_stage) * potential_dot_acoustic_LDDRK
            endif

            !! DK DK this should be vectorized
            potential_acoustic = potential_acoustic + beta_LDDRK(i_stage) * potential_acoustic_LDDRK
          endif

          if( time_stepping_scheme == 3 ) then
            !! DK DK this should be vectorized
            potential_dot_dot_acoustic = potential_dot_dot_acoustic * rmass_inverse_acoustic
            potential_dot_dot_acoustic_rk(:,i_stage) = deltat * potential_dot_dot_acoustic(:)
            potential_dot_acoustic_rk(:,i_stage) = deltat * potential_dot_acoustic(:)

            if( i_stage==1 .or. i_stage==2 .or. i_stage==3 ) then
              if( i_stage == 1 )weight_rk = 0.5d0
              if( i_stage == 2 )weight_rk = 0.5d0
              if( i_stage == 3 )weight_rk = 1.0d0

              if( i_stage==1 ) then
!! DK DK this should be vectorized
                potential_dot_acoustic_init_rk = potential_dot_acoustic
                potential_acoustic_init_rk = potential_acoustic
              endif
!! DK DK this should be vectorized
              potential_dot_acoustic(:) = potential_dot_acoustic_init_rk(:) + &
                                          weight_rk * potential_dot_dot_acoustic_rk(:,i_stage)
              potential_acoustic(:) = potential_acoustic_init_rk(:) + weight_rk * potential_dot_acoustic_rk(:,i_stage)
            else if( i_stage==4 ) then
!! DK DK this should be vectorized
              potential_dot_acoustic(:) = potential_dot_acoustic_init_rk(:) + &
                                          1.0d0 / 6.0d0 * ( potential_dot_dot_acoustic_rk(:,1) + &
                                                            2.0d0 * potential_dot_dot_acoustic_rk(:,2) + &
                                                            2.0d0 * potential_dot_dot_acoustic_rk(:,3) + &
                                                            potential_dot_dot_acoustic_rk(:,4) )

!! DK DK this should be vectorized
              potential_acoustic(:) = potential_acoustic_init_rk(:) + &
                                      1.0d0 / 6.0d0 * ( potential_dot_acoustic_rk(:,1) + &
                                                        2.0d0 * potential_dot_acoustic_rk(:,2) + &
                                                        2.0d0 * potential_dot_acoustic_rk(:,3) + &
                                                        potential_dot_acoustic_rk(:,4) )
            endif
          endif
        endif ! of if(any_acoustic)
      else ! GPU_MODE
        if(any_acoustic) call compute_forces_acoustic_GPU()
      endif ! GPU_MODE

! *********************************************************
! ************* main solver for the gravitoacoustic elements
! *********************************************************
! only SIMULATION_TYPE == 1, time_stepping_scheme == 1, and no PML or STACEY yet
! NO MIX OF ACOUSTIC AND GRAVITOACOUTIC ELEMENTS
! NO COUPLING TO ELASTIC AND POROELASTIC SIDES
! *********************************************************
      if( .not. GPU_MODE ) then
        if( (any_gravitoacoustic) ) then
          if( time_stepping_scheme==1 ) then
            ! Newmark time scheme
            !! DK DK this should be vectorized
            potential_gravitoacoustic = potential_gravitoacoustic + deltat*potential_dot_gravitoacoustic + &
                                        deltatsquareover2*potential_dot_dot_gravitoacoustic
            potential_dot_gravitoacoustic = potential_dot_gravitoacoustic + &
                                            deltatover2*potential_dot_dot_gravitoacoustic
            potential_gravito = potential_gravito + deltat*potential_dot_gravito + &
                                deltatsquareover2*potential_dot_dot_gravito
            potential_dot_gravito = potential_dot_gravito + deltatover2*potential_dot_dot_gravito
          else
            stop 'Only time_stepping_scheme=1 for gravitoacoustic'
          endif

          potential_dot_dot_gravitoacoustic = ZERO
          potential_dot_dot_gravito = ZERO

! Impose displacements from boundary forcing here
! because at this step the displacement (potentials) values
! are already equal to value at n+1
! equivalent to free surface condition
! the contour integral u.n is computed after compute_forces_gravitoacoustic
! *********************************************************
! ** impose displacement from acoustic forcing at a rigid boundary
! ** force potential_dot_dot_gravito by displacement
! *********************************************************
          if( ACOUSTIC_FORCING ) then
            call add_acoustic_forcing_at_rigid_boundary_gravitoacoustic()
          endif ! end ACOUSTIC_FORCING !

! free surface for a gravitoacoustic medium
!!! to be coded !!!
!      if ( nelem_acoustic_surface > 0 ) then
!        call enforce_acoustic_free_surface(potential_dot_dot_gravitoacoustic,potential_dot_gravitoacoustic, &
!                                          potential_gravitoacoustic)

!        if(SIMULATION_TYPE == 3) then ! Adjoint calculation
!          call enforce_acoustic_free_surface(b_potential_dot_dot_gravitoacoustic,b_potential_dot_gravitoacoustic, &
!                                            b_potential_gravitoacoustic)
!        endif
!      endif

! *********************************************************
! ************* compute forces for the gravitoacoustic elements
! *********************************************************

          call compute_forces_gravitoacoustic(potential_dot_dot_gravitoacoustic,potential_dot_gravitoacoustic, &
                       potential_gravitoacoustic, potential_dot_dot_gravito, &
                       potential_gravito,.false.,PML_BOUNDARY_CONDITIONS)

          if( (mod(it,100)==0) ) then
            iglob=iglobzero
            write(*,*)it,Nsql,gravityl, &
                      maxval(potential_dot_dot_gravito),potential_dot_dot_gravito(iglob), &
                      maxval(potential_dot_dot_gravitoacoustic),potential_dot_dot_gravitoacoustic(iglob)
          endif
        endif ! end of test if any gravitoacoustic element

! *********************************************************
! ************* add coupling with the elastic side
! *********************************************************

! *********************************************************
! ************* add coupling with the poroelastic side
! *********************************************************

! ************************************************************************************
! ************************************ add force source
! ************************************************************************************

! assembling potential_dot_dot for gravitoacoustic elements
!#ifdef USE_MPI
!    if ( nproc > 1 .and. any_acoustic .and. ninterface_acoustic > 0) then
!      call assemble_MPI_vector_ac(potential_dot_dot_gravitoacoustic)
!
!    endif
!
!#endif

! ************************************************************************************
! ************* multiply by the inverse of the mass matrix and update velocity
! ************************************************************************************

        if( (any_gravitoacoustic) ) then
          if( time_stepping_scheme == 1 ) then
          !! DK DK this should be vectorized
          potential_dot_dot_gravitoacoustic = potential_dot_dot_gravitoacoustic * rmass_inverse_gravitoacoustic
          potential_dot_gravitoacoustic = potential_dot_gravitoacoustic + &
                                          deltatover2*potential_dot_dot_gravitoacoustic

!! line below already done in compute_forces_gravitoacoustic, because necessary
!! for the computation of potential_dot_dot_gravitoacoustic
!      potential_dot_dot_gravito = potential_dot_dot_gravito * rmass_inverse_gravito
          potential_dot_gravito = potential_dot_gravito + deltatover2*potential_dot_dot_gravito
        else
          stop 'Only time_stepping_scheme = 1 implemented for gravitoacoustic case'
        endif

! free surface for an acoustic medium
!      if ( nelem_acoustic_surface > 0 ) then
!        call enforce_acoustic_free_surface(potential_dot_dot_gravitoacoustic,potential_dot_gravitoacoustic, &
!                                        potential_gravitoacoustic)
!
!        if(SIMULATION_TYPE == 3) then
!          call enforce_acoustic_free_surface(b_potential_dot_dot_gravitoacoustic,b_potential_dot_gravitoacoustic, &
!                                          b_potential_gravitoacoustic)
!        endif
!
!      endif
!
      ! update the potential field (use a new array here) for coupling terms
!      potential_gravitoacoustic_adj_coupling = potential_gravitoacoustic &
!                          + deltat*potential_dot_gravitoacoustic &
!                          + deltatsquareover2*potential_dot_dot_gravitoacoustic

        endif ! of if(any_gravitoacoustic)
      else ! GPU_MODE
        if ((any_gravitoacoustic)) call exit_mpi('gravitoacoustic not implemented in GPU MODE yet')
      endif

! *********************************************************
! ************* main solver for the elastic elements
! *********************************************************

      if(.not. GPU_MODE ) then
        if(any_elastic) then
          if( SIMULATION_TYPE == 1 ) then
            call compute_forces_viscoelastic(accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old, &
                                             x_source(1),z_source(1),f0(1),v0x_left(1,it),v0z_left(1,it), &
                                             v0x_right(1,it),v0z_right(1,it),v0x_bot(1,it),v0z_bot(1,it), &
                                             t0x_left(1,it),t0z_left(1,it),t0x_right(1,it),t0z_right(1,it), &
                                             t0x_bot(1,it),t0z_bot(1,it),count_left,count_right,count_bottom, &
                                             PML_BOUNDARY_CONDITIONS,e1,e11,e13)
          endif

          if( SIMULATION_TYPE == 3 ) then
            !ZN currently we do not support plane wave source in adjoint inversion
            call compute_forces_viscoelastic(accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old, &
                                             x_source(1),z_source(1),f0(1),v0x_left(1,it),v0z_left(1,it), &
                                             v0x_right(1,it),v0z_right(1,it),v0x_bot(1,it),v0z_bot(1,it), &
                                             t0x_left(1,it),t0z_left(1,it),t0x_right(1,it),t0z_right(1,it), &
                                             t0x_bot(1,it),t0z_bot(1,it),count_left,count_right,count_bottom, &
                                             PML_BOUNDARY_CONDITIONS,e1,e11,e13)

            if( PML_BOUNDARY_CONDITIONS ) then
              it_temp = NSTEP-it+1
              call rebuild_value_on_PML_interface_viscoelastic(it_temp)
            endif

            call compute_forces_viscoelastic_backward(b_accel_elastic,b_displ_elastic,b_displ_elastic_old, &
                                                      PML_BOUNDARY_CONDITIONS,e1,e11,e13)
            if( PML_BOUNDARY_CONDITIONS ) then
              it_temp = NSTEP-it+1
              call rebuild_value_on_PML_interface_viscoelastic(it_temp)
            endif
          endif

        endif !if(any_elastic)

        ! *********************************************************
        ! ************* add coupling with the acoustic side
        ! *********************************************************
        if( coupled_acoustic_elastic ) then
          if( SIMULATION_TYPE == 1 ) then
            call compute_coupling_viscoelastic_ac()
          endif

          if( SIMULATION_TYPE == 3 ) then
            call compute_coupling_viscoelastic_ac()
            call compute_coupling_viscoelastic_ac_backward()
          endif
        endif
        ! ****************************************************************************
        ! ************* add coupling with the poroelastic side
        ! ****************************************************************************
        if( coupled_acoustic_elastic ) then
          if( SIMULATION_TYPE == 1 ) then
            call compute_coupling_viscoelastic_po()
          endif

          if( SIMULATION_TYPE == 3 ) then
            call compute_coupling_viscoelastic_po()
            call compute_coupling_viscoelastic_po_backward()
          endif
        endif

        if( AXISYM ) then
          call enforce_zero_radial_displacements_on_the_axis()
        endif

        ! ************************************************************************************
        ! ************************************ add force source
        ! ************************************************************************************
        if( any_elastic ) then
          if( .not. initialfield ) then
            if( SIMULATION_TYPE == 1 ) then
              call compute_add_sources_viscoelastic(accel_elastic,it,i_stage)
            endif

            if( SIMULATION_TYPE == 3 ) then   ! adjoint and backward wavefield
              call compute_add_sources_viscoelastic_adjoint()
              call compute_add_sources_viscoelastic(b_accel_elastic,NSTEP-it+1,stage_time_scheme-i_stage+1)
            endif ! SIMULATION_TYPE == 3 ! adjoint and backward wavefield

            !<NOISE_TOMOGRAPHY
            ! inject wavefield sources for noise simulations
            if( NOISE_TOMOGRAPHY == 1 ) then
              call  add_point_source_noise()
            else if( NOISE_TOMOGRAPHY == 2 ) then
              call add_surface_movie_noise(accel_elastic)
            else if( NOISE_TOMOGRAPHY == 3 ) then
              if( .not. save_everywhere ) then
                call add_surface_movie_noise(b_accel_elastic)
              endif
            endif
            !>NOISE_TOMOGRAPHY
          endif ! if not using an initial field
        endif !if(any_elastic)

        if( AXISYM ) then
          call enforce_zero_radial_displacements_on_the_axis()
        endif
        ! ************************************************************************************
        ! ************************************ assembling accel_elastic for elastic elements
        ! ************************************************************************************
#ifdef USE_MPI
        if( nproc > 1 .and. any_elastic .and. ninterface_elastic > 0 ) then
          if( time_stepping_scheme == 2 ) then
            if( i_stage==1 .and. it == 1 .and. (.not. initialfield) ) then
              veloc_elastic_LDDRK_temp = veloc_elastic
              call assemble_MPI_vector_el(veloc_elastic)
            endif
          endif

          call assemble_MPI_vector_el(accel_elastic)

          if( SIMULATION_TYPE == 3 ) then
            call assemble_MPI_vector_el(b_accel_elastic)
          endif
        endif
#endif

        if( PML_BOUNDARY_CONDITIONS ) then
          if( any_elastic .and. nglob_interface > 0 ) then
            if( SAVE_FORWARD .and. SIMULATION_TYPE == 1)then
              do i = 1, nglob_interface
                write(71)accel_elastic(1,point_interface(i)),accel_elastic(2,point_interface(i)),&
                         accel_elastic(3,point_interface(i)),&
                         veloc_elastic(1,point_interface(i)),veloc_elastic(2,point_interface(i)),&
                         veloc_elastic(3,point_interface(i)),&
                         displ_elastic(1,point_interface(i)),displ_elastic(2,point_interface(i)),&
                         displ_elastic(3,point_interface(i))
              enddo
            endif

            if( SIMULATION_TYPE == 3 ) then
              do i = 1, nglob_interface
                b_accel_elastic(1,point_interface(i)) = pml_interface_history_accel(1,i,NSTEP-it+1)
                b_accel_elastic(2,point_interface(i)) = pml_interface_history_accel(2,i,NSTEP-it+1)
                b_accel_elastic(3,point_interface(i)) = pml_interface_history_accel(3,i,NSTEP-it+1)
              enddo
            endif
          endif
        endif

        ! ************************************************************************************
        ! ************* multiply by the inverse of the mass matrix and update velocity
        ! ************************************************************************************
        if( any_elastic ) then
          !! DK DK this should be vectorized
          accel_elastic(1,:) = accel_elastic(1,:) * rmass_inverse_elastic_one
          accel_elastic(2,:) = accel_elastic(2,:) * rmass_inverse_elastic_one
          accel_elastic(3,:) = accel_elastic(3,:) * rmass_inverse_elastic_three

          if( time_stepping_scheme == 1 ) then
            !! DK DK this should be vectorized
            veloc_elastic = veloc_elastic + deltatover2 * accel_elastic
          endif

          if( time_stepping_scheme == 2 ) then
            !! DK DK this should be vectorized
            veloc_elastic_LDDRK = alpha_LDDRK(i_stage) * veloc_elastic_LDDRK + deltat * accel_elastic
            displ_elastic_LDDRK = alpha_LDDRK(i_stage) * displ_elastic_LDDRK + deltat * veloc_elastic
            if( i_stage==1 .and. it == 1 .and. (.not. initialfield) ) then
              veloc_elastic_LDDRK_temp = veloc_elastic_LDDRK_temp + beta_LDDRK(i_stage) * veloc_elastic_LDDRK
              veloc_elastic = veloc_elastic_LDDRK_temp
            else
              veloc_elastic = veloc_elastic + beta_LDDRK(i_stage) * veloc_elastic_LDDRK
            endif
            displ_elastic = displ_elastic + beta_LDDRK(i_stage) * displ_elastic_LDDRK
          endif

          if( time_stepping_scheme == 3 ) then
            !! DK DK this should be vectorized
            accel_elastic_rk(1,:,i_stage) = deltat * accel_elastic(1,:)
            accel_elastic_rk(2,:,i_stage) = deltat * accel_elastic(2,:)
            accel_elastic_rk(3,:,i_stage) = deltat * accel_elastic(3,:)

            veloc_elastic_rk(1,:,i_stage) = deltat * veloc_elastic(1,:)
            veloc_elastic_rk(2,:,i_stage) = deltat * veloc_elastic(2,:)
            veloc_elastic_rk(3,:,i_stage) = deltat * veloc_elastic(3,:)

            if( i_stage==1 .or. i_stage==2 .or. i_stage==3 ) then

              if(i_stage == 1)weight_rk = 0.5d0
              if(i_stage == 2)weight_rk = 0.5d0
              if(i_stage == 3)weight_rk = 1.0d0

              if( i_stage==1 ) then
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

            else if( i_stage==4 ) then
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

          if( SIMULATION_TYPE == 3 ) then
            !! DK DK this should be vectorized
            b_accel_elastic(1,:) = b_accel_elastic(1,:) * rmass_inverse_elastic_one(:)
            b_accel_elastic(2,:) = b_accel_elastic(2,:) * rmass_inverse_elastic_one(:)
            b_accel_elastic(3,:) = b_accel_elastic(3,:) * rmass_inverse_elastic_three(:)

            b_veloc_elastic = b_veloc_elastic + b_deltatover2*b_accel_elastic
          endif

        endif !if(any_elastic)

      else ! GPU_MODE
        if(any_elastic)  call compute_forces_elastic_GPU()
      endif

! ******************************************************************************************************************
! ************* main solver for the poroelastic elements: first the solid (u_s) then the fluid (w)
! ******************************************************************************************************************
      if( .not. GPU_MODE) then
        if( any_poroelastic ) then
          !--------------------------------------------------------------------------------------------
          ! implement viscous attenuation for poroelastic media
          !--------------------------------------------------------------------------------------------
          if( ATTENUATION_PORO_FLUID_PART ) then
            call compute_attenuation_poro_fluid_part()
          endif

          if( SIMULATION_TYPE == 3 ) then
            ! if inviscid fluid, comment the reading and uncomment the zeroing
            !     read(23,rec=NSTEP-it+1) b_viscodampx
            !     read(24,rec=NSTEP-it+1) b_viscodampz
            b_viscodampx(:) = ZERO
            b_viscodampz(:) = ZERO
          endif

          call compute_forces_poro_solid(f0(1))
          call compute_forces_poro_fluid(f0(1))

          if( SAVE_FORWARD .and. SIMULATION_TYPE == 1 ) then
            ! if inviscid fluid, comment
            !     write(23,rec=it) b_viscodampx
            !     write(24,rec=it) b_viscodampz
          endif

        endif !if(any_poroelastic) then

        ! *********************************************************
        ! ************* add coupling with the acoustic side
        ! *********************************************************
        if( coupled_acoustic_poro ) then
          call compute_coupling_poro_ac()
        endif

        ! **********************************************************
        ! ************* add coupling with the elastic side
        ! **********************************************************
        if( coupled_elastic_poro ) then
          call compute_coupling_poro_viscoelastic()
        endif

        ! ***********************************************************
        ! ******************************** add force source
        ! ***********************************************************
        if( any_poroelastic ) then
          ! --- add the source if it is a collocated force
          if( .not. initialfield ) then
            if( SIMULATION_TYPE == 1 ) then
              call compute_add_sources_poro(accels_poroelastic,accelw_poroelastic,it,i_stage)
            endif

            if( SIMULATION_TYPE == 3 ) then   ! adjoint and backward wavefield
!ZNZN the add force source for adjoint simulation in poro medium are inside compute_forces_poro_solid
!ZNZN and compute_forces_poro_fluid,
              call compute_add_sources_poro(b_accels_poroelastic,b_accelw_poroelastic,NSTEP-it+1,stage_time_scheme-i_stage+1)
            endif
          endif
        endif
        ! ***********************************************************
        ! ******************************** ! assembling accels_proelastic & accelw_poroelastic for poroelastic elements
        ! ***********************************************************
#ifdef USE_MPI
        if( nproc > 1 .and. any_poroelastic .and. ninterface_poroelastic > 0 ) then
          call assemble_MPI_vector_po(accels_poroelastic,accelw_poroelastic)
          if( SIMULATION_TYPE == 3 ) then
            call assemble_MPI_vector_po(b_accels_poroelastic,b_accelw_poroelastic)
          endif
        endif
#endif
        ! ************************************************************************************
        ! ************* multiply by the inverse of the mass matrix and update velocity
        ! ************************************************************************************
        if( any_poroelastic ) then
          if( time_stepping_scheme == 1 ) then
            accels_poroelastic(1,:) = accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
            accels_poroelastic(2,:) = accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)
            velocs_poroelastic = velocs_poroelastic + deltatover2*accels_poroelastic

            accelw_poroelastic(1,:) = accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
            accelw_poroelastic(2,:) = accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)
            velocw_poroelastic = velocw_poroelastic + deltatover2*accelw_poroelastic

            if( SIMULATION_TYPE == 3 ) then
              b_accels_poroelastic(1,:) = b_accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
              b_accels_poroelastic(2,:) = b_accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)
              b_velocs_poroelastic = b_velocs_poroelastic + b_deltatover2*b_accels_poroelastic

              b_accelw_poroelastic(1,:) = b_accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
              b_accelw_poroelastic(2,:) = b_accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)
              b_velocw_poroelastic = b_velocw_poroelastic + b_deltatover2*b_accelw_poroelastic
            endif
          endif

          if( time_stepping_scheme == 2 ) then
            accels_poroelastic(1,:) = accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
            accels_poroelastic(2,:) = accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)

            velocs_poroelastic_LDDRK = alpha_LDDRK(i_stage) * velocs_poroelastic_LDDRK + deltat * accels_poroelastic
            displs_poroelastic_LDDRK = alpha_LDDRK(i_stage) * displs_poroelastic_LDDRK + deltat * velocs_poroelastic
            velocs_poroelastic = velocs_poroelastic + beta_LDDRK(i_stage) * velocs_poroelastic_LDDRK
            displs_poroelastic = displs_poroelastic + beta_LDDRK(i_stage) * displs_poroelastic_LDDRK

            accelw_poroelastic(1,:) = accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
            accelw_poroelastic(2,:) = accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)

            velocw_poroelastic_LDDRK = alpha_LDDRK(i_stage) * velocw_poroelastic_LDDRK + deltat * accelw_poroelastic
            displw_poroelastic_LDDRK = alpha_LDDRK(i_stage) * displw_poroelastic_LDDRK + deltat * velocw_poroelastic
            velocw_poroelastic = velocw_poroelastic + beta_LDDRK(i_stage) * velocw_poroelastic_LDDRK
            displw_poroelastic = displw_poroelastic + beta_LDDRK(i_stage) * displw_poroelastic_LDDRK
          endif

          if( time_stepping_scheme == 3 ) then
            accels_poroelastic(1,:) = accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
            accels_poroelastic(2,:) = accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)

            accels_poroelastic_rk(1,:,i_stage) = deltat * accels_poroelastic(1,:)
            accels_poroelastic_rk(2,:,i_stage) = deltat * accels_poroelastic(2,:)
            velocs_poroelastic_rk(1,:,i_stage) = deltat * velocs_poroelastic(1,:)
            velocs_poroelastic_rk(2,:,i_stage) = deltat * velocs_poroelastic(2,:)

            accelw_poroelastic(1,:) = accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
            accelw_poroelastic(2,:) = accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)

            accelw_poroelastic_rk(1,:,i_stage) = deltat * accelw_poroelastic(1,:)
            accelw_poroelastic_rk(2,:,i_stage) = deltat * accelw_poroelastic(2,:)
            velocw_poroelastic_rk(1,:,i_stage) = deltat * velocw_poroelastic(1,:)
            velocw_poroelastic_rk(2,:,i_stage) = deltat * velocw_poroelastic(2,:)

            if( i_stage==1 .or. i_stage==2 .or. i_stage==3 ) then
              if( i_stage == 1 )weight_rk = 0.5d0
              if( i_stage == 2 )weight_rk = 0.5d0
              if( i_stage == 3 )weight_rk = 1.0d0

              if( i_stage==1 ) then
                velocs_poroelastic_initial_rk(1,:) = velocs_poroelastic(1,:)
                velocs_poroelastic_initial_rk(2,:) = velocs_poroelastic(2,:)
                displs_poroelastic_initial_rk(1,:) = displs_poroelastic(1,:)
                displs_poroelastic_initial_rk(2,:) = displs_poroelastic(2,:)

                velocw_poroelastic_initial_rk(1,:) = velocw_poroelastic(1,:)
                velocw_poroelastic_initial_rk(2,:) = velocw_poroelastic(2,:)
                displw_poroelastic_initial_rk(1,:) = displw_poroelastic(1,:)
                displw_poroelastic_initial_rk(2,:) = displw_poroelastic(2,:)
              endif

              velocs_poroelastic(1,:) = velocs_poroelastic_initial_rk(1,:) + weight_rk * accels_poroelastic_rk(1,:,i_stage)
              velocs_poroelastic(2,:) = velocs_poroelastic_initial_rk(2,:) + weight_rk * accels_poroelastic_rk(2,:,i_stage)
              displs_poroelastic(1,:) = displs_poroelastic_initial_rk(1,:) + weight_rk * velocs_poroelastic_rk(1,:,i_stage)
              displs_poroelastic(2,:) = displs_poroelastic_initial_rk(2,:) + weight_rk * velocs_poroelastic_rk(2,:,i_stage)

              velocw_poroelastic(1,:) = velocw_poroelastic_initial_rk(1,:) + weight_rk * accelw_poroelastic_rk(1,:,i_stage)
              velocw_poroelastic(2,:) = velocw_poroelastic_initial_rk(2,:) + weight_rk * accelw_poroelastic_rk(2,:,i_stage)
              displw_poroelastic(1,:) = displw_poroelastic_initial_rk(1,:) + weight_rk * velocw_poroelastic_rk(1,:,i_stage)
              displw_poroelastic(2,:) = displw_poroelastic_initial_rk(2,:) + weight_rk * velocw_poroelastic_rk(2,:,i_stage)

            else if( i_stage==4 ) then

              velocs_poroelastic(1,:) = velocs_poroelastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
              (accels_poroelastic_rk(1,:,1) + 2.0d0 * accels_poroelastic_rk(1,:,2) + &
              2.0d0 * accels_poroelastic_rk(1,:,3) + accels_poroelastic_rk(1,:,4))

              velocs_poroelastic(2,:) = velocs_poroelastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
              (accels_poroelastic_rk(2,:,1) + 2.0d0 * accels_poroelastic_rk(2,:,2) + &
               2.0d0 * accels_poroelastic_rk(2,:,3) + accels_poroelastic_rk(2,:,4))

              displs_poroelastic(1,:) = displs_poroelastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
              (velocs_poroelastic_rk(1,:,1) + 2.0d0 * velocs_poroelastic_rk(1,:,2) + &
               2.0d0 * velocs_poroelastic_rk(1,:,3) + velocs_poroelastic_rk(1,:,4))

              displs_poroelastic(2,:) = displs_poroelastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
              (velocs_poroelastic_rk(2,:,1) + 2.0d0 * velocs_poroelastic_rk(2,:,2) + &
               2.0d0 * velocs_poroelastic_rk(2,:,3) + velocs_poroelastic_rk(2,:,4))

              velocw_poroelastic(1,:) = velocw_poroelastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
              (accelw_poroelastic_rk(1,:,1) + 2.0d0 * accelw_poroelastic_rk(1,:,2) + &
               2.0d0 * accelw_poroelastic_rk(1,:,3) + accelw_poroelastic_rk(1,:,4))

              velocw_poroelastic(2,:) = velocw_poroelastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
              (accelw_poroelastic_rk(2,:,1) + 2.0d0 * accelw_poroelastic_rk(2,:,2) + &
               2.0d0 * accelw_poroelastic_rk(2,:,3) + accelw_poroelastic_rk(2,:,4))

              displw_poroelastic(1,:) = displw_poroelastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
              (velocw_poroelastic_rk(1,:,1) + 2.0d0 * velocw_poroelastic_rk(1,:,2) + &
               2.0d0 * velocw_poroelastic_rk(1,:,3) + velocw_poroelastic_rk(1,:,4))

              displw_poroelastic(2,:) = displw_poroelastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
              (velocw_poroelastic_rk(2,:,1) + 2.0d0 * velocw_poroelastic_rk(2,:,2) + &
               2.0d0 * velocw_poroelastic_rk(2,:,3) + velocw_poroelastic_rk(2,:,4))
            endif
          endif
        endif !if(any_poroelastic)

        !*******************************************************************************
        ! assembling the displacements on the elastic-poro boundaries
        !*******************************************************************************

        !
        ! Explanation of the code below, from Christina Morency and Yang Luo, January 2012:
        !
        ! Coupled elastic-poroelastic simulations imply continuity of traction and
        ! displacement at the interface.
        ! For the traction we pass on both sides n*(T + Te)/2 , that is, the average
        ! between the total stress (from the poroelastic part) and the elastic stress.
        ! For the displacement, we enforce its continuity in the assembling stage,
        ! realizing that continuity of displacement correspond to the continuity of
        ! the acceleration we have:
        !
        ! accel_elastic = rmass_inverse_elastic * force_elastic
        ! accels_poroelastic = rmass_s_inverse_poroelastic * force_poroelastic
        !
        ! Therefore, continuity of acceleration gives
        !
        ! accel = (force_elastic + force_poroelastic)/
        !     (1/rmass_inverse_elastic + 1/rmass_inverse_poroelastic)
        !
        ! Then
        !
        ! accel_elastic = accel
        ! accels_poroelastic = accel
        ! accelw_poroelastic = 0
        !
        ! From there, the velocity and displacement are updated.
        ! Note that force_elastic and force_poroelastic are the right hand sides of
        ! the equations we solve, that is, the acceleration terms before the
        ! division by the inverse of the mass matrices. This is why in the code below
        ! we first need to recover the accelerations (which are then
        ! the right hand sides forces) and the velocities before the update.
        !
        ! This implementation highly helped stability especially with unstructured meshes.
        !
        if( coupled_elastic_poro ) then
          call compute_coupling_poro_viscoelastic_for_stabilization()
        endif
     else !GPU_MODE
       if(any_poroelastic)   call exit_mpi('poroelastic not implemented in GPU MODE yet')
     endif
   enddo !LDDRK or RK

! ********************************************************************************************
!                       reading lastframe for adjoint/kernels calculation
! ********************************************************************************************
   if( it == 1 .and. SIMULATION_TYPE == 3 ) then
     ! acoustic medium
     if( any_acoustic ) then
       write(outputname,'(a,i6.6,a)') 'lastframe_acoustic',myrank,'.bin'
       open(unit=55,file='OUTPUT_FILES/'//outputname,status='old',action='read',form='unformatted')

         read(55) b_potential_acoustic
         read(55) b_potential_dot_acoustic
         read(55) b_potential_dot_dot_acoustic

       close(55)

       if( GPU_MODE ) then
         ! transfers fields onto GPU
         call transfer_b_fields_ac_to_device(NGLOB_AB,b_potential_acoustic, &
                                             b_potential_dot_acoustic,      &
                                             b_potential_dot_dot_acoustic,  &
                                             Mesh_pointer)
       else
         ! free surface for an acoustic medium
         if( nelem_acoustic_surface > 0 ) then
           call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                              b_potential_acoustic)
         endif
       endif
     endif

     ! elastic medium
     if( any_elastic ) then
       write(outputname,'(a,i6.6,a)') 'lastframe_elastic',myrank,'.bin'
       open(unit=55,file='OUTPUT_FILES/'//outputname,status='old',action='read',form='unformatted')
       if( p_sv ) then !P-SV waves
           read(55) b_displ_elastic
           read(55) b_veloc_elastic
           read(55) b_accel_elastic
         if( GPU_MODE ) then
           b_displ_2D(1,:) = b_displ_elastic(1,:)
           b_displ_2D(2,:) = b_displ_elastic(2,:)
           b_veloc_2D(1,:) = b_veloc_elastic(1,:)
           b_veloc_2D(2,:) = b_veloc_elastic(2,:)
           b_accel_2D(1,:) = b_accel_elastic(1,:)
           b_accel_2D(2,:) = b_accel_elastic(2,:)
           call transfer_b_fields_to_device(NDIM*NGLOB_AB,b_displ_2D,b_veloc_2D,b_accel_2D,Mesh_pointer)
         endif

       else !SH (membrane) waves

           read(55) b_displ_elastic
           read(55) b_veloc_elastic
           read(55) b_accel_elastic

         b_displ_elastic(1,:) = 0._CUSTOM_REAL
         b_displ_elastic(3,:) = 0._CUSTOM_REAL
         b_veloc_elastic(1,:) = 0._CUSTOM_REAL
         b_veloc_elastic(3,:) = 0._CUSTOM_REAL
         b_accel_elastic(1,:) = 0._CUSTOM_REAL
         b_accel_elastic(3,:) = 0._CUSTOM_REAL
       endif
       close(55)
     endif

     ! poroelastic medium
     if(any_poroelastic) then
       write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_s',myrank,'.bin'
       open(unit=55,file='OUTPUT_FILES/'//outputname,status='old',action='read',form='unformatted')
       write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_w',myrank,'.bin'
       open(unit=56,file='OUTPUT_FILES/'//outputname,status='old',action='read',form='unformatted')
       do j=1,nglob
         read(55) (b_displs_poroelastic(i,j), i=1,NDIM), &
                  (b_velocs_poroelastic(i,j), i=1,NDIM), &
                  (b_accels_poroelastic(i,j), i=1,NDIM)
         read(56) (b_displw_poroelastic(i,j), i=1,NDIM), &
                  (b_velocw_poroelastic(i,j), i=1,NDIM), &
                  (b_accelw_poroelastic(i,j), i=1,NDIM)
       enddo
       close(55)
       close(56)
     endif
   endif ! if(it == 1 .and. SIMULATION_TYPE == 3)

!<NOISE_TOMOGRAPHY
   if( NOISE_TOMOGRAPHY == 1 ) then
      call save_surface_movie_noise()
   else if( NOISE_TOMOGRAPHY == 2 .and. save_everywhere ) then
      call save_surface_movie_noise()
   else if( NOISE_TOMOGRAPHY == 3 .and. save_everywhere ) then
     if( it==1 ) open(unit=500,file='OUTPUT_FILES/NOISE_TOMOGRAPHY/phi',access='direct', &
                      recl=nglob*CUSTOM_REAL,action='write',iostat=ios)
     if( ios /= 0 ) write(*,*) 'Error retrieving ensemble forward wavefield.'
     if( p_sv ) then
       call exit_mpi('P-SV case not yet implemented.')
     else
       read(unit=500,rec=NSTEP-it+1) b_displ_elastic(2,:)
     endif
   endif
!>NOISE_TOMOGRAPHY

! ********************************************************************************************
!                                      kernels calculation
! ********************************************************************************************
   !*******************************************************************************
   ! GPU_MODE
   !*******************************************************************************
   if( GPU_MODE ) then
     ! Kernel calculation
     if(SIMULATION_TYPE == 3 ) then
       if( any_acoustic ) call compute_kernels_acoustic_cuda(Mesh_pointer,deltatf)
       if( any_elastic ) call compute_kernels_elastic_cuda(Mesh_pointer,deltatf)

       if( APPROXIMATE_HESS_KL ) then
         ! computes contribution to density and bulk modulus kernel
         call compute_kernels_hess_cuda(Mesh_pointer,any_elastic,any_acoustic)
       endif

       ! Kernel transfer
       if( it == NSTEP ) then
         if( any_acoustic ) then
           call transfer_kernels_ac_to_host(Mesh_pointer,rho_ac_kl,kappa_ac_kl,NSPEC_AB)
          rhop_ac_kl(:,:,:) = rho_ac_kl(:,:,:) + kappa_ac_kl(:,:,:)
          alpha_ac_kl(:,:,:) = TWO *  kappa_ac_kl(:,:,:)

         endif

         if( any_elastic ) then
           call transfer_kernels_el_to_host(Mesh_pointer,rho_kl,mu_kl,kappa_kl,NSPEC_AB)
           ! Multiply each kernel point with the local coefficient
           do ispec = 1, nspec
             if( elastic(ispec) ) then
               do j = 1, NGLLZ
                 do i = 1, NGLLX
                   iglob = ibool(i,j,ispec)
                   if( .not. assign_external_model ) then
                     mul_global(iglob) = poroelastcoef(2,1,kmato(ispec))
                     kappal_global(iglob) = poroelastcoef(3,1,kmato(ispec)) - &
                                            4._CUSTOM_REAL*mul_global(iglob)/3._CUSTOM_REAL
                     rhol_global(iglob) = density(1,kmato(ispec))
                   else
                     rhol_global(iglob)   = rhoext(i,j,ispec)
                     mul_global(iglob)    = rhoext(i,j,ispec)*vsext(i,j,ispec)*vsext(i,j,ispec)
                     kappal_global(iglob) = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec) - &
                                            4._CUSTOM_REAL*mul_global(iglob)/3._CUSTOM_REAL
                   endif

                   rho_kl(i,j,ispec) = - rhol_global(iglob) * rho_kl(i,j,ispec)
                   mu_kl(i,j,ispec) =  - TWO * mul_global(iglob) * mu_kl(i,j,ispec)
                   kappa_kl(i,j,ispec) = - kappal_global(iglob) * kappa_kl(i,j,ispec)
                 enddo
               enddo
             endif
           enddo
         endif  !!End elastic
       endif  !! End NSTEP
     endif  !! End Sim 3

     if( mod(it-1,subsamp_seismos) == 0 .and. SIMULATION_TYPE == 1 ) then
       seismo_current = seismo_current + 1
       if( nrecloc > 0 ) then
          if (USE_TRICK_FOR_BETTER_PRESSURE) then

           call compute_seismograms_cuda(Mesh_pointer,seismotype,sisux,sisuz,seismo_current,&
                                                       NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos, &
                                                       any_elastic_glob,any_acoustic_glob,1)
          else
           call compute_seismograms_cuda(Mesh_pointer,seismotype,sisux,sisuz,seismo_current,&
                                                       NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,&
                                                       any_elastic_glob,any_acoustic_glob,0)
          endif
       endif
     endif

     ! Fields transfer for imaging
     if( (output_color_image .and. ( (mod(it,NSTEP_BETWEEN_OUTPUT_IMAGES) == 0 .or. it == 5)) .or. it == NSTEP) ) then
       if( any_acoustic ) &
         call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic,potential_dot_acoustic, &
                                             potential_dot_dot_acoustic,Mesh_pointer)
       if( any_elastic ) then
         call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ_2D,veloc_2D,accel_2D,Mesh_pointer)
         displ_elastic(1,:) = displ_2D(1,:)
         veloc_elastic(1,:) = veloc_2D(1,:)
         accel_elastic(1,:) = accel_2D(1,:)
         displ_elastic(3,:) = displ_2D(2,:)
         veloc_elastic(3,:) = veloc_2D(2,:)
         accel_elastic(3,:) = accel_2D(2,:)
      endif
     endif !If transfer
   endif !If GPU Mode

   !*******************************************************************************
   ! CPU_MODE
   !*******************************************************************************
   if( .not. GPU_MODE ) then
     !----  compute kinetic and potential energy
     if( output_energy ) then
       call compute_energy()
#ifdef USE_MPI
       call MPI_REDUCE(kinetic_energy, kinetic_energy_total, 1, CUSTOM_MPI_TYPE, MPI_SUM, 0, MPI_COMM_WORLD, ier)
       call MPI_REDUCE(potential_energy, potential_energy_total, 1, CUSTOM_MPI_TYPE, MPI_SUM, 0, MPI_COMM_WORLD, ier)
#else
       kinetic_energy_total = kinetic_energy
       potential_energy_total = potential_energy
#endif

       ! save kinetic, potential and total energy for this time step in external file
       if(myrank == 0) write(IOUT_ENERGY,*) real(dble(it-1)*deltat - t0,4),real(kinetic_energy_total,4), &
                       real(potential_energy_total,4),real(kinetic_energy_total + potential_energy_total,4)
     endif

     !----  display time step and max of norm of displacement
     if( mod(it,NSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP ) then
       call check_stability()
     endif

     !---- loop on all the receivers to compute and store the seismograms
     if( mod(it-1,subsamp_seismos) == 0 ) then
       call write_seismograms()
     endif

     ! kernels calculation
     if( SIMULATION_TYPE == 3 ) then
       if( any_acoustic ) then
         call compute_kernels_ac()
       endif !if(any_acoustic)

       if( any_elastic ) then
         call compute_kernels_el()
       endif

       if(any_poroelastic) then
         call compute_kernels_po()
       endif ! if(any_poroelastic)
     endif ! if(SIMULATION_TYPE == 3)
   endif !Not GPU_MODE

   !
   !----  display results at given time steps
   !
   if( mod(it,NSTEP_BETWEEN_OUTPUT_IMAGES) == 0 .or. it == 5 .or. it == NSTEP ) then
     ! write kernel files
     if(SIMULATION_TYPE == 3 .and. it == NSTEP) then
        call save_adjoint_kernels()
     endif

     !<NOISE_TOMOGRAPHY
     if(.not. GPU_MODE ) then
       if(NOISE_TOMOGRAPHY == 3 .and. output_wavefields_noise) then

         !load ensemble forward source
         inquire(unit=500,exist=ex,opened=od)
         if (.not. od) &
           open(unit=500,file='OUTPUT_FILES/NOISE_TOMOGRAPHY/eta',access='direct', &
           recl=nglob*CUSTOM_REAL,action='write',iostat=ios)
         read(unit=500,rec=it) surface_movie_y_noise

         !load product of fwd, adj wavefields
         call spec2glob(nspec,nglob,ibool,rho_kl,noise_output_rhokl)

         !write text file
         noise_output_array(1,:) = surface_movie_y_noise(:) * mask_noise(:)
         noise_output_array(2,:) = b_displ_elastic(2,:)
         noise_output_array(3,:) = accel_elastic(2,:)
         noise_output_array(4,:) = rho_k(:)
         noise_output_array(5,:) = noise_output_rhokl(:)
         write(noise_output_file,"('OUTPUT_FILES/snapshot_all_',i6.6)") it
         call snapshots_noise(noise_output_ncol,nglob,noise_output_file,noise_output_array)
       endif
     endif
     !>NOISE_TOMOGRAPHY

     ! ********************************************************************************************
     ! output_postscript_snapshot
     ! ********************************************************************************************
     if( output_postscript_snapshot ) then
       call write_postscript_snapshot()
     endif

     ! ********************************************************************************************
     ! display color image
     ! ********************************************************************************************
     if( output_color_image ) then
       call write_color_image_snaphot()
     endif  ! of display images at a given time step

     ! ********************************************************************************************
     ! dump the full (local) wavefield to a file
     ! note: in the case of MPI, in the future it would be more convenient to output a single file
     !       rather than one for each myrank
     ! ********************************************************************************************
     if( output_wavefield_dumps ) then
       call write_wavefield_dumps()
     endif  ! of display wavefield dumps at a given time step
   endif

   !----  save temporary or final seismograms
   ! suppress seismograms if we generate traces of the run for analysis with "ParaVer", because time consuming
   if( mod(it,NSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP ) then
     call write_seismograms_to_file(x_source(1),z_source(1))
     seismo_offset = seismo_offset + seismo_current
     seismo_current = 0
   endif  ! of display images at a given time step

  enddo ! end of the main time loop

! *********************************************************
! ************* END MAIN LOOP OVER THE TIME STEPS *********
! *********************************************************!
  !----  formats
  400 format(/1x,41('=')/,' =  T i m e  e v o l u t i o n  l o o p  ='/1x,41('=')/)

 end subroutine
