
subroutine iterate_time()

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par

  implicit none

integer i,j,ispec,i_source,iglob,irec,k

#ifdef USE_MPI
  include "precision.h"
#endif


if (myrank == 0) write(IOUT,400)



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


     call update_displacement_precondition_newmark()

      if (.NOT. GPU_MODE) then

      if (AXISYM) then
        do ispec=1,nspec
          if (elastic(ispec) .and. is_on_the_axis(ispec)) then
            do j = 1,NGLLZ
              do i = 1,NGLJ
                if (abs(coord(1,ibool(i,j,ispec))) < TINYVAL) then
                  displ_elastic(1,ibool(i,j,ispec))=ZERO
                endif
              enddo
            enddo
          endif
        enddo
      endif

    if(any_acoustic) then
      ! free surface for an acoustic medium
      if ( nelem_acoustic_surface > 0 ) then
        call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                          potential_acoustic)

        if(SIMULATION_TYPE == 3) then ! Adjoint calculation
          call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                            b_potential_acoustic)
        endif
      endif

! *********************************************************
! ************* compute forces for the acoustic elements
! *********************************************************

      call compute_forces_acoustic(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                   potential_acoustic,potential_acoustic_old,PML_BOUNDARY_CONDITIONS)

      if( SIMULATION_TYPE == 3 ) then

       if(PML_BOUNDARY_CONDITIONS)then
          do ispec = 1,nspec
            do i = 1, NGLLX
              do j = 1, NGLLZ
                if(acoustic(ispec) .and. is_pml(ispec))then
                  b_potential_dot_dot_acoustic(ibool(i,j,ispec)) = 0.
                  b_potential_dot_acoustic(ibool(i,j,ispec)) = 0.
                  b_potential_acoustic(ibool(i,j,ispec)) = 0.
                endif
               enddo
            enddo
          enddo
       endif

       if(PML_BOUNDARY_CONDITIONS)then
         if(any_acoustic .and. nglob_interface > 0)then
           do i = 1, nglob_interface
             b_potential_dot_acoustic(point_interface(i)) = pml_interface_history_potential_dot(i,NSTEP-it+1)
             b_potential_acoustic(point_interface(i)) = pml_interface_history_potential(i,NSTEP-it+1)
           enddo
         endif
       endif

       call compute_forces_acoustic_backward(b_potential_dot_dot_acoustic,b_potential_acoustic)

       if(PML_BOUNDARY_CONDITIONS)then
          do ispec = 1,nspec
            do i = 1, NGLLX
              do j = 1, NGLLZ
                if(acoustic(ispec) .and. is_pml(ispec))then
                  b_potential_dot_acoustic(ibool(i,j,ispec)) = 0.
                  b_potential_acoustic(ibool(i,j,ispec)) = 0.
                endif
               enddo
            enddo
          enddo
       endif

       if(PML_BOUNDARY_CONDITIONS)then
         if(any_acoustic .and. nglob_interface > 0)then
           do i = 1, nglob_interface
             b_potential_dot_acoustic(point_interface(i)) = pml_interface_history_potential_dot(i,NSTEP-it+1)
             b_potential_acoustic(point_interface(i)) = pml_interface_history_potential(i,NSTEP-it+1)
           enddo
         endif
       endif

      endif


    ! *********************************************************
    ! ************* add acoustic forcing at a rigid boundary
    ! *********************************************************
    if(ACOUSTIC_FORCING) then
      ! loop on all the forced edges

     do inum = 1,nelem_acforcing

        ispec = numacforcing(inum)

        !--- left acoustic forcing boundary
        if(codeacforcing(IEDGE4,inum)) then

           i = 1

           do j = 1,NGLLZ

              ! acoustic spectral element
              if(acoustic(ispec)) then
                 iglob = ibool(i,j,ispec)

                 xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
                 zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xgamma**2 + zgamma**2)
                 nx = - zgamma / jacobian1D
                 nz = + xgamma / jacobian1D

                 weight = jacobian1D * wzgll(j)

          ! define displacement components which will force the boundary

            if(PML_BOUNDARY_CONDITIONS) then
              if(is_PML(ispec)) then
              displ_x = 0
              displ_z = 0
              else
              call acoustic_forcing_boundary()
              endif

            else

            call acoustic_forcing_boundary()
            endif

          ! compute dot product
          displ_n = displ_x*nx + displ_z*nz

          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n

              endif  !end of acoustic
           enddo

        endif  !  end of left acoustic forcing boundary

        !--- right acoustic forcing boundary
        if(codeacforcing(IEDGE2,inum)) then

           i = NGLLX

           do j = 1,NGLLZ

              ! acoustic spectral element
              if(acoustic(ispec)) then
                 iglob = ibool(i,j,ispec)

                 xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
                 zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xgamma**2 + zgamma**2)
                 nx = + zgamma / jacobian1D
                 nz = - xgamma / jacobian1D

                 weight = jacobian1D * wzgll(j)

          ! define displacement components which will force the boundary

            if(PML_BOUNDARY_CONDITIONS) then

              if(is_PML(ispec)) then
              displ_x = 0
              displ_z = 0
              else
              call acoustic_forcing_boundary()
              endif

            else

            call acoustic_forcing_boundary()

            endif

          ! compute dot product
          displ_n = displ_x*nx + displ_z*nz

          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n

              endif  !end of acoustic
           enddo

        endif  !  end of right acoustic forcing boundary

        !--- bottom acoustic forcing boundary
        if(codeacforcing(IEDGE1,inum)) then

           j = 1

           do i = 1,NGLLX

              ! acoustic spectral element
              if(acoustic(ispec)) then
                 iglob = ibool(i,j,ispec)

                 xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                 zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xxi**2 + zxi**2)
                 nx = + zxi / jacobian1D
                 nz = - xxi / jacobian1D

                 weight = jacobian1D * wxgll(i)

          ! define displacement components which will force the boundary

            if(PML_BOUNDARY_CONDITIONS) then
              if(is_PML(ispec)) then
              displ_x = 0
              displ_z = 0
              else
              call acoustic_forcing_boundary()
              endif

            else
            call acoustic_forcing_boundary()

            endif

          ! compute dot product
          displ_n = displ_x*nx + displ_z*nz

          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n

              endif  !end of acoustic
           enddo

        endif  !  end of bottom acoustic forcing boundary

        !--- top acoustic forcing boundary
        if(codeacforcing(IEDGE3,inum)) then

           j = NGLLZ

           do i = 1,NGLLX

              ! acoustic spectral element
              if(acoustic(ispec)) then
                 iglob = ibool(i,j,ispec)

                 xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                 zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xxi**2 + zxi**2)
                 nx = - zxi / jacobian1D
                 nz = + xxi / jacobian1D

                 weight = jacobian1D * wxgll(i)

          ! define displacement components which will force the boundary

            if(PML_BOUNDARY_CONDITIONS) then

              if(is_PML(ispec)) then
              displ_x = 0
              displ_z = 0
              else
              call acoustic_forcing_boundary()
              endif

            else

            call acoustic_forcing_boundary()

            endif

          ! compute dot product
          displ_n = displ_x*nx + displ_z*nz

          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n

              endif  !end of acoustic
           enddo

        endif  !  end of top acoustic forcing boundary

     enddo

     endif ! of if ACOUSTIC_FORCING

    endif ! end of test if any acoustic element

! *********************************************************
! ************* add coupling with the elastic side
! *********************************************************

    if(coupled_acoustic_elastic) then

      if(SIMULATION_TYPE == 1)then
        call compute_coupling_acoustic_el(displ_elastic,displ_elastic_old,potential_dot_dot_acoustic, PML_BOUNDARY_CONDITIONS)
      endif

      if(SIMULATION_TYPE == 3)then

        accel_elastic_adj_coupling2 = - accel_elastic_adj_coupling

        call compute_coupling_acoustic_el(accel_elastic_adj_coupling2,displ_elastic_old,potential_dot_dot_acoustic,&
                                          PML_BOUNDARY_CONDITIONS)

        call compute_coupling_acoustic_el(b_displ_elastic,b_displ_elastic_old,b_potential_dot_dot_acoustic,.false.)

      endif

    endif

! *********************************************************
! ************* add coupling with the poroelastic side
! *********************************************************

    if(coupled_acoustic_poro) then

      ! loop on all the coupling edges
      do inum = 1,num_fluid_poro_edges

        ! get the edge of the acoustic element
        ispec_acoustic = fluid_poro_acoustic_ispec(inum)
        iedge_acoustic = fluid_poro_acoustic_iedge(inum)

        ! get the corresponding edge of the poroelastic element
        ispec_poroelastic = fluid_poro_poroelastic_ispec(inum)
        iedge_poroelastic = fluid_poro_poroelastic_iedge(inum)

        ! implement 1D coupling along the edge
        do ipoin1D = 1,NGLLX

          ! get point values for the poroelastic side, which matches our side in the inverse direction
          i = ivalue_inverse(ipoin1D,iedge_poroelastic)
          j = jvalue_inverse(ipoin1D,iedge_poroelastic)
          iglob = ibool(i,j,ispec_poroelastic)

          displ_x = displs_poroelastic(1,iglob)
          displ_z = displs_poroelastic(2,iglob)

          phil = porosity(kmato(ispec_poroelastic))
          displw_x = displw_poroelastic(1,iglob)
          displw_z = displw_poroelastic(2,iglob)

          if(SIMULATION_TYPE == 3) then
            b_displ_x = b_displs_poroelastic(1,iglob)
            b_displ_z = b_displs_poroelastic(2,iglob)

            b_displw_x = b_displw_poroelastic(1,iglob)
            b_displw_z = b_displw_poroelastic(2,iglob)

            ! new definition of adjoint displacement and adjoint potential
            displ_x = -accels_poroelastic_adj_coupling(1,iglob)
            displ_z = -accels_poroelastic_adj_coupling(2,iglob)

            displw_x = -accelw_poroelastic_adj_coupling(1,iglob)
            displw_z = -accelw_poroelastic_adj_coupling(2,iglob)
          endif

          ! get point values for the acoustic side
          ! get point values for the acoustic side
          i = ivalue(ipoin1D,iedge_acoustic)
          j = jvalue(ipoin1D,iedge_acoustic)
          iglob = ibool(i,j,ispec_acoustic)

          ! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
          ! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
          ! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
          ! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
          ! Blackwell Science, page 110, equation (4.60).
          if(iedge_acoustic == ITOP)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = - zxi / jacobian1D
            nz = + xxi / jacobian1D
            weight = jacobian1D * wxgll(i)
          else if(iedge_acoustic == IBOTTOM)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D
            weight = jacobian1D * wxgll(i)
          else if(iedge_acoustic ==ILEFT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = - zgamma / jacobian1D
            nz = + xgamma / jacobian1D
            weight = jacobian1D * wzgll(j)
          else if(iedge_acoustic ==IRIGHT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = + zgamma / jacobian1D
            nz = - xgamma / jacobian1D
            weight = jacobian1D * wzgll(j)
          endif

          ! compute dot product [u_s + w]*n
          displ_n = (displ_x + displw_x)*nx + (displ_z + displw_z)*nz

          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n

          if(SIMULATION_TYPE == 3) then
            b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) &
                   + weight*((b_displ_x + b_displw_x)*nx + (b_displ_z + b_displw_z)*nz)
          endif

        enddo

      enddo

    endif


! ************************************************************************************
! ************************************ add force source
! ************************************************************************************

    if(any_acoustic) then

      ! --- add the source
      if(.not. initialfield) then

        do i_source=1,NSOURCES
          ! if this processor core carries the source and the source element is acoustic
          if (is_proc_source(i_source) == 1 .and. acoustic(ispec_selected_source(i_source))) then

            ! collocated force
            ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
            ! the sign is negative because pressure p = - Chi_dot_dot therefore we need
            ! to add minus the source to Chi_dot_dot to get plus the source in pressure
            if(source_type(i_source) == 1) then
              if(SIMULATION_TYPE == 1) then
                ! forward wavefield
                do j = 1,NGLLZ
                  do i = 1,NGLLX
                    iglob = ibool(i,j,ispec_selected_source(i_source))
                    hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
                    potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                            - source_time_function(i_source,it,i_stage)*hlagrange &
                                               / kappastore(i,j,ispec_selected_source(i_source))
                  enddo
                enddo
              else
                ! backward wavefield
                do j = 1,NGLLZ
                  do i = 1,NGLLX
                    iglob = ibool(i,j,ispec_selected_source(i_source))
                    hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
                    b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) &
                                          - source_time_function(i_source,NSTEP-it+1,stage_time_scheme-i_stage+1)*hlagrange &
                                            / kappastore(i,j,ispec_selected_source(i_source))
                  enddo
                enddo
              endif

            ! moment tensor
            else if(source_type(i_source) == 2) then
              call exit_MPI('cannot have moment tensor source in acoustic element')

            endif
          endif ! if this processor core carries the source and the source element is acoustic
        enddo ! do i_source=1,NSOURCES

        if(SIMULATION_TYPE == 3) then   ! adjoint wavefield
          irec_local = 0
          do irec = 1,nrec
            !   add the source (only if this proc carries the source)
            if (myrank == which_proc_receiver(irec)) then

              irec_local = irec_local + 1
              if (acoustic(ispec_selected_rec(irec))) then
                ! add source array
                do j=1,NGLLZ
                  do i=1,NGLLX
                    iglob = ibool(i,j,ispec_selected_rec(irec))
                    potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                  + adj_sourcearrays(irec_local,NSTEP-it+1,1,i,j) / kappastore(i,j,ispec_selected_rec(irec))
                  enddo
                enddo
              endif ! if element acoustic

            endif ! if this processor core carries the adjoint source
          enddo ! irec = 1,nrec
        endif ! SIMULATION_TYPE == 3 adjoint wavefield

      endif ! if not using an initial field

    endif !if(any_acoustic)


! assembling potential_dot_dot for acoustic elements
#ifdef USE_MPI
    if ( nproc > 1 .and. any_acoustic .and. ninterface_acoustic > 0) then
      call assemble_MPI_vector_ac(potential_dot_dot_acoustic)

     if(time_stepping_scheme == 2)then
      if(i_stage==1 .and. it == 1 .and. (.not. initialfield))then
       potential_dot_acoustic_temp = potential_dot_acoustic
       call assemble_MPI_vector_ac(potential_dot_acoustic)
      endif
     endif

      if ( SIMULATION_TYPE == 3) then
        call assemble_MPI_vector_ac(b_potential_dot_dot_acoustic)

      endif

    endif

#endif

      if(PML_BOUNDARY_CONDITIONS .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1)then
       if(any_acoustic .and. nglob_interface > 0)then
        do i = 1, nglob_interface
          write(72)potential_dot_dot_acoustic(point_interface(i)),&
                   potential_dot_acoustic(point_interface(i)),&
                   potential_acoustic(point_interface(i))
        enddo
       endif
      endif

     if(SIMULATION_TYPE == 3)then
       if(PML_BOUNDARY_CONDITIONS)then
         if(any_acoustic .and. nglob_interface > 0)then
           do i = 1, nglob_interface
             b_potential_dot_dot_acoustic(point_interface(i)) = pml_interface_history_potential_dot_dot(i,NSTEP-it+1)
           enddo
         endif
       endif
     endif

! ************************************************************************************
! ************* multiply by the inverse of the mass matrix and update velocity
! ************************************************************************************

    if(any_acoustic) then
      if(time_stepping_scheme == 1)then
!! DK DK this should be vectorized
      potential_dot_dot_acoustic = potential_dot_dot_acoustic * rmass_inverse_acoustic
      potential_dot_acoustic = potential_dot_acoustic + deltatover2*potential_dot_dot_acoustic
      endif

      if(time_stepping_scheme == 2)then

!! DK DK this should be vectorized
        potential_dot_dot_acoustic = potential_dot_dot_acoustic * rmass_inverse_acoustic

        potential_dot_acoustic_LDDRK = alpha_LDDRK(i_stage) * potential_dot_acoustic_LDDRK &
                                       + deltat * potential_dot_dot_acoustic

        potential_acoustic_LDDRK = alpha_LDDRK(i_stage) * potential_acoustic_LDDRK &
                                   + deltat*potential_dot_acoustic

        if(i_stage==1 .and. it == 1 .and. (.not. initialfield))then
!! DK DK this should be vectorized
        potential_dot_acoustic_temp = potential_dot_acoustic_temp &
                                      + beta_LDDRK(i_stage) * potential_dot_acoustic_LDDRK
        potential_dot_acoustic = potential_dot_acoustic_temp
        else
        potential_dot_acoustic = potential_dot_acoustic + beta_LDDRK(i_stage) * potential_dot_acoustic_LDDRK
        endif

!! DK DK this should be vectorized
        potential_acoustic = potential_acoustic + beta_LDDRK(i_stage) * potential_acoustic_LDDRK

      endif

      if(time_stepping_scheme == 3)then

!! DK DK this should be vectorized
        potential_dot_dot_acoustic = potential_dot_dot_acoustic * rmass_inverse_acoustic

        potential_dot_dot_acoustic_rk(:,i_stage) = deltat * potential_dot_dot_acoustic(:)
        potential_dot_acoustic_rk(:,i_stage) = deltat * potential_dot_acoustic(:)

        if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then

        if(i_stage == 1)weight_rk = 0.5d0
        if(i_stage == 2)weight_rk = 0.5d0
        if(i_stage == 3)weight_rk = 1.0d0

        if(i_stage==1)then

!! DK DK this should be vectorized
        potential_dot_acoustic_init_rk = potential_dot_acoustic
        potential_acoustic_init_rk = potential_acoustic

        endif

!! DK DK this should be vectorized
        potential_dot_acoustic(:) = potential_dot_acoustic_init_rk(:) + weight_rk * potential_dot_dot_acoustic_rk(:,i_stage)
        potential_acoustic(:) = potential_acoustic_init_rk(:) + weight_rk * potential_dot_acoustic_rk(:,i_stage)

        else if(i_stage==4)then

!! DK DK this should be vectorized
        potential_dot_acoustic(:) = potential_dot_acoustic_init_rk(:) + 1.0d0 / 6.0d0 * &
        (potential_dot_dot_acoustic_rk(:,1) + 2.0d0 * potential_dot_dot_acoustic_rk(:,2) + &
         2.0d0 * potential_dot_dot_acoustic_rk(:,3) + potential_dot_dot_acoustic_rk(:,4))

!! DK DK this should be vectorized
        potential_acoustic(:) = potential_acoustic_init_rk(:) + 1.0d0 / 6.0d0 * &
        (potential_dot_acoustic_rk(:,1) + 2.0d0 * potential_dot_acoustic_rk(:,2) + &
         2.0d0 * potential_dot_acoustic_rk(:,3) + potential_dot_acoustic_rk(:,4))

        endif

      endif

      if(SIMULATION_TYPE == 3)then
!! DK DK this should be vectorized
        b_potential_dot_dot_acoustic = b_potential_dot_dot_acoustic * rmass_inverse_acoustic
        b_potential_dot_acoustic = b_potential_dot_acoustic + b_deltatover2*b_potential_dot_dot_acoustic
      endif


      ! free surface for an acoustic medium
      if ( nelem_acoustic_surface > 0 ) then
        call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                        potential_acoustic)

        if(SIMULATION_TYPE == 3) then
          call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                          b_potential_acoustic)
        endif

      endif

      ! update the potential field (use a new array here) for coupling terms
      potential_acoustic_adj_coupling = potential_acoustic &
                          + deltat*potential_dot_acoustic &
                          + deltatsquareover2*potential_dot_dot_acoustic

    endif ! of if(any_acoustic)


   else ! GPU_MODE

    if(any_acoustic) call compute_forces_acoustic_GPU()

   endif


! *********************************************************
! ************* main solver for the gravitoacoustic elements
! *********************************************************
! only SIMULATION_TYPE == 1, time_stepping_scheme == 1, and no PML or STACEY yet
! NO MIX OF ACOUSTIC AND GRAVITOACOUTIC ELEMENTS
! NO COUPLING TO ELASTIC AND POROELASTIC SIDES
! *********************************************************
!-----------------------------------------

 if (.NOT. GPU_MODE) then

    if ((any_gravitoacoustic)) then

      if(time_stepping_scheme==1)then
      ! Newmark time scheme
!! DK DK this should be vectorized
      potential_gravitoacoustic = potential_gravitoacoustic &
                          + deltat*potential_dot_gravitoacoustic &
                          + deltatsquareover2*potential_dot_dot_gravitoacoustic
      potential_dot_gravitoacoustic = potential_dot_gravitoacoustic &
                              + deltatover2*potential_dot_dot_gravitoacoustic
      potential_gravito = potential_gravito &
                          + deltat*potential_dot_gravito &
                          + deltatsquareover2*potential_dot_dot_gravito
      potential_dot_gravito = potential_dot_gravito &
                              + deltatover2*potential_dot_dot_gravito
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
    if(ACOUSTIC_FORCING) then

      ! loop on all the forced edges

     do inum = 1,nelem_acforcing

        ispec = numacforcing(inum)

        !--- left acoustic forcing boundary
        if(codeacforcing(IEDGE4,inum)) then

           i = 1

           do j = 1,NGLLZ

              ! acoustic spectral element
              if(gravitoacoustic(ispec)) then
                 iglob = ibool(i,j,ispec)
                 xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
                 zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xgamma**2 + zgamma**2)
                 nx = - zgamma / jacobian1D
                 nz = + xgamma / jacobian1D

                 weight = jacobian1D * wzgll(j)

          ! define displacement components which will force the boundary

            if(PML_BOUNDARY_CONDITIONS) then

              if(is_PML(ispec)) then
              displ_x = 0
              displ_z = 0
              else
              call acoustic_forcing_boundary()
              endif

            else

            call acoustic_forcing_boundary()

            endif

! compute displacement at this point
        ! derivative along x
        tempx1l = 0._CUSTOM_REAL
        do k = 1,NGLLX
          hp1 = hprime_xx(i,k)
          iglob = ibool(k,j,ispec)
          tempx1l = tempx1l + potential_gravitoacoustic(iglob)*hp1
        enddo

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + potential_gravitoacoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

          ! if external density model
          if(assign_external_model)then
            if(CUSTOM_REAL == SIZE_REAL) then
              rhol = sngl(rhoext(i,j,ispec))
              gravityl = sngl(gravityext(i,j,ispec))
            else
              rhol = rhoext(i,j,ispec)
              gravityl = gravityext(i,j,ispec)
            endif
          endif

! impose potential_gravito in order to have z displacement equal to forced
! value
          iglob = ibool(i,j,ispec)
          displ_n = displ_x*nx + displ_z*nz
        if (abs(nz) > TINYVAL) then
          potential_gravito(iglob) = (rhol*displ_n - &
          (tempx1l*xizl + tempx2l*gammazl)*nz - (tempx1l*xixl + tempx2l*gammaxl)*nx)/ &
          (0.0 - gravityl*nz)
        else
          write(*,*) 'STOP : forcing surface element along z',i,j,ispec,iglob,nx,nz
          stop
        endif

          ! compute dot product
          displ_n = displ_x*nx + displ_z*nz
          potential_dot_dot_gravito(iglob) = potential_dot_dot_gravito(iglob) - rhol*weight*displ_n

              endif  !end of gravitoacoustic
           enddo

        endif  !  end of left acoustic forcing boundary

        !--- right acoustic forcing boundary
        if(codeacforcing(IEDGE2,inum)) then

           i = NGLLX

           do j = 1,NGLLZ

              ! acoustic spectral element
              if(gravitoacoustic(ispec)) then
                 iglob = ibool(i,j,ispec)
                 xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
                 zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xgamma**2 + zgamma**2)
                 nx = + zgamma / jacobian1D
                 nz = - xgamma / jacobian1D

                 weight = jacobian1D * wzgll(j)

          ! define displacement components which will force the boundary

            if(PML_BOUNDARY_CONDITIONS) then

              if(is_PML(ispec)) then
              displ_x = 0
              displ_z = 0
              else
              call acoustic_forcing_boundary()
              endif

            else

            call acoustic_forcing_boundary()

            endif

! compute displacement at this point
        ! derivative along x
        tempx1l = 0._CUSTOM_REAL
        do k = 1,NGLLX
          hp1 = hprime_xx(i,k)
          iglob = ibool(k,j,ispec)
          tempx1l = tempx1l + potential_gravitoacoustic(iglob)*hp1
        enddo

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + potential_gravitoacoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

          ! if external density model
          if(assign_external_model)then
            if(CUSTOM_REAL == SIZE_REAL) then
              rhol = sngl(rhoext(i,j,ispec))
              gravityl = sngl(gravityext(i,j,ispec))
            else
              rhol = rhoext(i,j,ispec)
              gravityl = gravityext(i,j,ispec)
            endif
          endif

! impose potential_gravito in order to have z displacement equal to forced
! value
          iglob = ibool(i,j,ispec)
          displ_n = displ_x*nx + displ_z*nz
        if (abs(nz) > TINYVAL) then
          potential_gravito(iglob) = (rhol*displ_n - &
          (tempx1l*xizl + tempx2l*gammazl)*nz - (tempx1l*xixl + tempx2l*gammaxl)*nx)/ &
          (0.0-gravityl*nz)
        else
          write(*,*) 'STOP : forcing surface element along z',i,j,ispec,iglob,nx,nz
          stop
        endif

          ! compute dot product
          displ_n = displ_x*nx + displ_z*nz
          potential_dot_dot_gravito(iglob) = potential_dot_dot_gravito(iglob) - rhol*weight*displ_n

              endif  !end of gravitoacoustic
           enddo

        endif  !  end of right acoustic forcing boundary

        !--- bottom acoustic forcing boundary
        if(codeacforcing(IEDGE1,inum)) then

           j = 1

           do i = 1,NGLLX

              ! acoustic spectral element
              if(gravitoacoustic(ispec)) then
                 iglob = ibool(i,j,ispec)
                 xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                 zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xxi**2 + zxi**2)
                 nx = + zxi / jacobian1D
                 nz = - xxi / jacobian1D

                 weight = jacobian1D * wxgll(i)

          ! define displacement components which will force the boundary

            if(PML_BOUNDARY_CONDITIONS) then

              if(is_PML(ispec)) then
              displ_x = 0
              displ_z = 0
              else
              call acoustic_forcing_boundary()
              endif

            else

            call acoustic_forcing_boundary()

            endif

! compute displacement at this point
        ! derivative along x
        tempx1l = 0._CUSTOM_REAL
        do k = 1,NGLLX
          hp1 = hprime_xx(i,k)
          iglob = ibool(k,j,ispec)
          tempx1l = tempx1l + potential_gravitoacoustic(iglob)*hp1
        enddo

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + potential_gravitoacoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

          ! if external density model
          if(assign_external_model)then
            if(CUSTOM_REAL == SIZE_REAL) then
              rhol = sngl(rhoext(i,j,ispec))
              gravityl = sngl(gravityext(i,j,ispec))
            else
              rhol = rhoext(i,j,ispec)
              gravityl = gravityext(i,j,ispec)
            endif
          endif

! impose potential_gravito in order to have z displacement equal to forced
! value
          iglob = ibool(i,j,ispec)
          displ_n = displ_x*nx + displ_z*nz
        if (abs(nz) > TINYVAL) then
          potential_gravito(iglob) = (rhol*displ_n - &
          (tempx1l*xizl + tempx2l*gammazl)*nz - (tempx1l*xixl + tempx2l*gammaxl)*nx)/ &
          (0.0 - gravityl*nz)
        else
          write(*,*) 'STOP : forcing surface element along z',i,j,ispec,iglob,nx,nz
          stop
        endif

          ! compute dot product
          displ_n = displ_x*nx + displ_z*nz
          potential_dot_dot_gravito(iglob) = potential_dot_dot_gravito(iglob) - rhol*weight*displ_n

              endif  !end of gravitoacoustic
           enddo

        endif  !  end of bottom acoustic forcing boundary

        !--- top acoustic forcing boundary
        if(codeacforcing(IEDGE3,inum)) then

           j = NGLLZ

           do i = 1,NGLLX

              ! acoustic spectral element
              if(gravitoacoustic(ispec)) then
                 iglob = ibool(i,j,ispec)
                 xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                 zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xxi**2 + zxi**2)
                 nx = - zxi / jacobian1D
                 nz = + xxi / jacobian1D

                 weight = jacobian1D * wxgll(i)

          ! define displacement components which will force the boundary

            if(PML_BOUNDARY_CONDITIONS) then

              if(is_PML(ispec)) then
              displ_x = 0
              displ_z = 0
              else
              call acoustic_forcing_boundary()
              endif

            else

            call acoustic_forcing_boundary()

            endif

! compute z displacement at this point
        ! derivative along x
        tempx1l = 0._CUSTOM_REAL
        do k = 1,NGLLX
          hp1 = hprime_xx(i,k)
          iglob = ibool(k,j,ispec)
          tempx1l = tempx1l + potential_gravitoacoustic(iglob)*hp1
        enddo

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + potential_gravitoacoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

          ! if external density model
          if(assign_external_model)then
            if(CUSTOM_REAL == SIZE_REAL) then
              rhol = sngl(rhoext(i,j,ispec))
              gravityl = sngl(gravityext(i,j,ispec))
            else
              rhol = rhoext(i,j,ispec)
              gravityl = gravityext(i,j,ispec)
              Nsql = Nsqext(i,j,ispec)
            endif
          endif

! impose potential_gravito in order to have z displacement equal to forced
! value on the boundary
!!!! Passe deux fois sur le meme iglob
!!!! Mais vrai pour tous les points partages entre deux elements
          iglob = ibool(i,j,ispec)
          displ_n = displ_x*nx + displ_z*nz
        if (abs(nz) > TINYVAL) then
          potential_gravito(iglob) = (rhol*displ_n - &
          (tempx1l*xizl + tempx2l*gammazl)*nz - (tempx1l*xixl + tempx2l*gammaxl)*nx)/ &
          (0.0 - gravityl*nz)
        else
          write(*,*) 'STOP : forcing surface element along z',i,j,ispec,iglob,nx,nz
          stop
        endif

          ! compute dot product
          displ_n = displ_x*nx + displ_z*nz
          potential_dot_dot_gravito(iglob) = potential_dot_dot_gravito(iglob) - rhol*weight*displ_n

              endif  !end of gravitoacoustic
           enddo

!       write(*,*) 'ispec detection =',ispec
!       if ((ispec==2000).and.(mod(it,100)==0)) then
       if ((ispec==800).and.(mod(it,100)==0)) then
!       if ((ispec==800)) then
       iglobzero=iglob
       write(*,*) ispec,it,Nsql,rhol,displ_n, &
       maxval(potential_dot_dot_gravito),potential_dot_dot_gravito(iglob), &
       maxval(potential_gravitoacoustic),potential_gravitoacoustic(iglob), &
       maxval(potential_gravito),potential_gravito(iglob)
       endif

        endif  !  end of top acoustic forcing boundary

     enddo

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

       if ((mod(it,100)==0)) then
         iglob=iglobzero
         write(*,*) it,Nsql,gravityl, &
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

    if((any_gravitoacoustic)) then
      if(time_stepping_scheme == 1)then
!! DK DK this should be vectorized

      potential_dot_dot_gravitoacoustic = potential_dot_dot_gravitoacoustic * rmass_inverse_gravitoacoustic
      potential_dot_gravitoacoustic = potential_dot_gravitoacoustic + deltatover2*potential_dot_dot_gravitoacoustic

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

   if (.NOT. GPU_MODE) then

    if(any_elastic) then

      call compute_forces_viscoelastic(accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old,x_source(1),z_source(1), &
               f0(1),v0x_left(1,it),v0z_left(1,it),v0x_right(1,it),v0z_right(1,it),v0x_bot(1,it),v0z_bot(1,it), &
               t0x_left(1,it),t0z_left(1,it),t0x_right(1,it),t0z_right(1,it),t0x_bot(1,it),t0z_bot(1,it), &
               count_left,count_right,count_bottom,PML_BOUNDARY_CONDITIONS,.false.)

      if(SIMULATION_TYPE == 3)then
       if(PML_BOUNDARY_CONDITIONS)then
          do ispec = 1,nspec
            do i = 1, NGLLX
              do j = 1, NGLLZ
                if(elastic(ispec) .and. is_pml(ispec))then
                  b_veloc_elastic(:,ibool(i,j,ispec)) = 0.
                  b_displ_elastic(:,ibool(i,j,ispec)) = 0.
                endif
               enddo
            enddo
          enddo
       endif

       if(PML_BOUNDARY_CONDITIONS)then
         if(any_elastic .and. nglob_interface > 0)then
           do i = 1, nglob_interface
             b_veloc_elastic(1,point_interface(i)) = pml_interface_history_veloc(1,i,NSTEP-it+1)
             b_veloc_elastic(2,point_interface(i)) = pml_interface_history_veloc(2,i,NSTEP-it+1)
             b_veloc_elastic(3,point_interface(i)) = pml_interface_history_veloc(3,i,NSTEP-it+1)
             b_displ_elastic(1,point_interface(i)) = pml_interface_history_displ(1,i,NSTEP-it+1)
             b_displ_elastic(2,point_interface(i)) = pml_interface_history_displ(2,i,NSTEP-it+1)
             b_displ_elastic(3,point_interface(i)) = pml_interface_history_displ(3,i,NSTEP-it+1)
           enddo
         endif
       endif

      call compute_forces_viscoelastic(b_accel_elastic,b_veloc_elastic,b_displ_elastic,&
               b_displ_elastic_old,x_source(1),z_source(1),f0(1),&
               v0x_left(1,it),v0z_left(1,it),v0x_right(1,it),v0z_right(1,it),v0x_bot(1,it),v0z_bot(1,it), &
               t0x_left(1,it),t0z_left(1,it),t0x_right(1,it),t0z_right(1,it),t0x_bot(1,it),t0z_bot(1,it), &
               count_left,count_right,count_bottom,.false.,.true.)

       if(PML_BOUNDARY_CONDITIONS)then
          do ispec = 1,nspec
            do i = 1, NGLLX
              do j = 1, NGLLZ
                if(elastic(ispec) .and. is_pml(ispec))then
                  b_veloc_elastic(:,ibool(i,j,ispec)) = 0.
                  b_displ_elastic(:,ibool(i,j,ispec)) = 0.
                endif
               enddo
            enddo
          enddo
       endif

       if(PML_BOUNDARY_CONDITIONS)then
         if(any_elastic .and. nglob_interface > 0)then
           do i = 1, nglob_interface
             b_veloc_elastic(1,point_interface(i)) = pml_interface_history_veloc(1,i,NSTEP-it+1)
             b_veloc_elastic(2,point_interface(i)) = pml_interface_history_veloc(2,i,NSTEP-it+1)
             b_veloc_elastic(3,point_interface(i)) = pml_interface_history_veloc(3,i,NSTEP-it+1)
             b_displ_elastic(1,point_interface(i)) = pml_interface_history_displ(1,i,NSTEP-it+1)
             b_displ_elastic(2,point_interface(i)) = pml_interface_history_displ(2,i,NSTEP-it+1)
             b_displ_elastic(3,point_interface(i)) = pml_interface_history_displ(3,i,NSTEP-it+1)
           enddo
         endif
       endif


      call compute_forces_viscoelastic_pre_kernel()

      endif


      if(anyabs .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. (.not. PML_BOUNDARY_CONDITIONS)) then
        !--- left absorbing boundary
        if(nspec_left >0) then
          do ispec = 1,nspec_left
            if(p_sv)then!P-SV waves
              do i=1,NGLLZ
                write(35) b_absorb_elastic_left(1,i,ispec,it)
              enddo
              do i=1,NGLLZ
                write(35) b_absorb_elastic_left(3,i,ispec,it)
              enddo
            else!SH (membrane) waves
              do i=1,NGLLZ
                write(35) b_absorb_elastic_left(2,i,ispec,it)
              enddo
            endif
          enddo
        endif

        !--- right absorbing boundary
        if(nspec_right >0) then
          do ispec = 1,nspec_right
            if(p_sv)then!P-SV waves
              do i=1,NGLLZ
                write(36) b_absorb_elastic_right(1,i,ispec,it)
              enddo
              do i=1,NGLLZ
                write(36) b_absorb_elastic_right(3,i,ispec,it)
              enddo
            else!SH (membrane) waves
              do i=1,NGLLZ
                write(36) b_absorb_elastic_right(2,i,ispec,it)
              enddo
            endif
          enddo
        endif

        !--- bottom absorbing boundary
        if(nspec_bottom >0) then
          do ispec = 1,nspec_bottom
            if(p_sv)then!P-SV waves
              do i=1,NGLLX
                write(37) b_absorb_elastic_bottom(1,i,ispec,it)
              enddo
              do i=1,NGLLX
                write(37) b_absorb_elastic_bottom(3,i,ispec,it)
              enddo
            else!SH (membrane) waves
              do i=1,NGLLX
                write(37) b_absorb_elastic_bottom(2,i,ispec,it)
              enddo
            endif
          enddo
        endif

        !--- top absorbing boundary
        if(nspec_top >0) then
          do ispec = 1,nspec_top
            if(p_sv)then!P-SV waves
              do i=1,NGLLX
                write(38) b_absorb_elastic_top(1,i,ispec,it)
              enddo
              do i=1,NGLLX
                write(38) b_absorb_elastic_top(3,i,ispec,it)
              enddo
            else!SH (membrane) waves
              do i=1,NGLLX
                write(38) b_absorb_elastic_top(2,i,ispec,it)
              enddo
            endif
          enddo
        endif

      endif ! if(anyabs .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1)

    endif !if(any_elastic)

! *********************************************************
! ************* add coupling with the acoustic side
! *********************************************************

    if(coupled_acoustic_elastic) call compute_coupling_viscoelastic_ac()

! ****************************************************************************
! ************* add coupling with the poroelastic side
! ****************************************************************************
    if(coupled_elastic_poro) then

      ! loop on all the coupling edges
      do inum = 1,num_solid_poro_edges

        ! get the edge of the elastic element
        ispec_elastic = solid_poro_elastic_ispec(inum)
        iedge_elastic = solid_poro_elastic_iedge(inum)

        ! get the corresponding edge of the poroelastic element
        ispec_poroelastic = solid_poro_poroelastic_ispec(inum)
        iedge_poroelastic = solid_poro_poroelastic_iedge(inum)

        ! implement 1D coupling along the edge
        do ipoin1D = 1,NGLLX

          ! get point values for the poroelastic side, which matches our side in the inverse direction
          i = ivalue_inverse(ipoin1D,iedge_poroelastic)
          j = jvalue_inverse(ipoin1D,iedge_poroelastic)
          iglob = ibool(i,j,ispec_poroelastic)

          ! get poroelastic domain paramters
          phil = porosity(kmato(ispec_poroelastic))
          tortl = tortuosity(kmato(ispec_poroelastic))
          !solid properties
          mul_s = poroelastcoef(2,1,kmato(ispec_poroelastic))
          kappal_s = poroelastcoef(3,1,kmato(ispec_poroelastic)) - 4._CUSTOM_REAL*mul_s/3._CUSTOM_REAL
          rhol_s = density(1,kmato(ispec_poroelastic))
          !fluid properties
          kappal_f = poroelastcoef(1,2,kmato(ispec_poroelastic))
          rhol_f = density(2,kmato(ispec_poroelastic))
          !frame properties
          mul_fr = poroelastcoef(2,3,kmato(ispec_poroelastic))
          kappal_fr = poroelastcoef(3,3,kmato(ispec_poroelastic)) - 4._CUSTOM_REAL*mul_fr/3._CUSTOM_REAL
          rhol_bar =  (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f
          !Biot coefficients for the input phi
          D_biot = kappal_s*(1._CUSTOM_REAL + phil*(kappal_s/kappal_f - 1._CUSTOM_REAL))
          H_biot = (kappal_s - kappal_fr)*(kappal_s - kappal_fr)/(D_biot - kappal_fr) + &
                kappal_fr + 4._CUSTOM_REAL*mul_fr/3._CUSTOM_REAL
          C_biot = kappal_s*(kappal_s - kappal_fr)/(D_biot - kappal_fr)
          M_biot = kappal_s*kappal_s/(D_biot - kappal_fr)
          mul_G = mul_fr
          lambdal_G = H_biot - 2._CUSTOM_REAL*mul_fr
          lambdalplus2mul_G = lambdal_G + TWO*mul_G

          ! derivative along x and along z for u_s and w
          dux_dxi = ZERO
          duz_dxi = ZERO

          dux_dgamma = ZERO
          duz_dgamma = ZERO

          dwx_dxi = ZERO
          dwz_dxi = ZERO

          dwx_dgamma = ZERO
          dwz_dgamma = ZERO

          if(SIMULATION_TYPE == 3) then
            b_dux_dxi = ZERO
            b_duz_dxi = ZERO

            b_dux_dgamma = ZERO
            b_duz_dgamma = ZERO

            b_dwx_dxi = ZERO
            b_dwz_dxi = ZERO

            b_dwx_dgamma = ZERO
            b_dwz_dgamma = ZERO
          endif

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displs_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displs_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displs_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displs_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)

            dwx_dxi = dwx_dxi + displw_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
            dwz_dxi = dwz_dxi + displw_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
            dwx_dgamma = dwx_dgamma + displw_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
            dwz_dgamma = dwz_dgamma + displw_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
            if(SIMULATION_TYPE == 3) then
              b_dux_dxi = b_dux_dxi + b_displs_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
              b_duz_dxi = b_duz_dxi + b_displs_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
              b_dux_dgamma = b_dux_dgamma + b_displs_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
              b_duz_dgamma = b_duz_dgamma + b_displs_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)

              b_dwx_dxi = b_dwx_dxi + b_displw_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
              b_dwz_dxi = b_dwz_dxi + b_displw_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
              b_dwx_dgamma = b_dwx_dgamma + b_displw_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
              b_dwz_dgamma = b_dwz_dgamma + b_displw_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
            endif
          enddo

          xixl = xix(i,j,ispec_poroelastic)
          xizl = xiz(i,j,ispec_poroelastic)
          gammaxl = gammax(i,j,ispec_poroelastic)
          gammazl = gammaz(i,j,ispec_poroelastic)

          ! derivatives of displacement
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

          dwx_dxl = dwx_dxi*xixl + dwx_dgamma*gammaxl
          dwx_dzl = dwx_dxi*xizl + dwx_dgamma*gammazl

          dwz_dxl = dwz_dxi*xixl + dwz_dgamma*gammaxl
          dwz_dzl = dwz_dxi*xizl + dwz_dgamma*gammazl

          if(SIMULATION_TYPE == 3) then
            b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
            b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl

            b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
            b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl

            b_dwx_dxl = b_dwx_dxi*xixl + b_dwx_dgamma*gammaxl
            b_dwx_dzl = b_dwx_dxi*xizl + b_dwx_dgamma*gammazl

            b_dwz_dxl = b_dwz_dxi*xixl + b_dwz_dgamma*gammaxl
            b_dwz_dzl = b_dwz_dxi*xizl + b_dwz_dgamma*gammazl
          endif
          ! compute stress tensor (include attenuation or anisotropy if needed)

          ! no attenuation
          sigma_xx = lambdalplus2mul_G*dux_dxl + lambdal_G*duz_dzl + C_biot*(dwx_dxl + dwz_dzl)
          sigma_xz = mul_G*(duz_dxl + dux_dzl)
          sigma_zz = lambdalplus2mul_G*duz_dzl + lambdal_G*dux_dxl + C_biot*(dwx_dxl + dwz_dzl)

          if(SIMULATION_TYPE == 3) then
            b_sigma_xx = lambdalplus2mul_G*b_dux_dxl + lambdal_G*b_duz_dzl + C_biot*(b_dwx_dxl + b_dwz_dzl)
            b_sigma_xz = mul_G*(b_duz_dxl + b_dux_dzl)
            b_sigma_zz = lambdalplus2mul_G*b_duz_dzl + lambdal_G*b_dux_dxl + C_biot*(b_dwx_dxl + b_dwz_dzl)
          endif
          ! get point values for the elastic domain, which matches our side in the inverse direction
          ii2 = ivalue(ipoin1D,iedge_elastic)
          jj2 = jvalue(ipoin1D,iedge_elastic)
          iglob = ibool(ii2,jj2,ispec_elastic)

          ! get elastic properties
          lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec_elastic))
          mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec_elastic))
          lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec_elastic))

          ! derivative along x and along z for u_s and w
          dux_dxi = ZERO
          duz_dxi = ZERO

          dux_dgamma = ZERO
          duz_dgamma = ZERO

          if(SIMULATION_TYPE == 3) then
            b_dux_dxi = ZERO
            b_duz_dxi = ZERO

            b_dux_dgamma = ZERO
            b_duz_dgamma = ZERO
          endif

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,jj2,ispec_elastic))*hprime_xx(ii2,k)
            duz_dxi = duz_dxi + displ_elastic(3,ibool(k,jj2,ispec_elastic))*hprime_xx(ii2,k)
            dux_dgamma = dux_dgamma + displ_elastic(1,ibool(ii2,k,ispec_elastic))*hprime_zz(jj2,k)
            duz_dgamma = duz_dgamma + displ_elastic(3,ibool(ii2,k,ispec_elastic))*hprime_zz(jj2,k)

            if(SIMULATION_TYPE == 3) then
              b_dux_dxi = b_dux_dxi + b_displ_elastic(1,ibool(k,jj2,ispec_elastic))*hprime_xx(ii2,k)
              b_duz_dxi = b_duz_dxi + b_displ_elastic(3,ibool(k,jj2,ispec_elastic))*hprime_xx(ii2,k)
              b_dux_dgamma = b_dux_dgamma + b_displ_elastic(1,ibool(ii2,k,ispec_elastic))*hprime_zz(jj2,k)
              b_duz_dgamma = b_duz_dgamma + b_displ_elastic(3,ibool(ii2,k,ispec_elastic))*hprime_zz(jj2,k)
            endif
          enddo

          xixl = xix(ii2,jj2,ispec_elastic)
          xizl = xiz(ii2,jj2,ispec_elastic)
          gammaxl = gammax(ii2,jj2,ispec_elastic)
          gammazl = gammaz(ii2,jj2,ispec_elastic)

          ! derivatives of displacement
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

          if(SIMULATION_TYPE == 3) then
            b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
            b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl

            b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
            b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl
          endif
          ! compute stress tensor
          ! full anisotropy
          if(kmato(ispec_elastic) == 2) then
            ! implement anisotropy in 2D
            if(assign_external_model) then
              c11 = c11ext(ii2,jj2,ispec_elastic)
              c13 = c13ext(ii2,jj2,ispec_elastic)
              c15 = c15ext(ii2,jj2,ispec_elastic)
              c33 = c33ext(ii2,jj2,ispec_elastic)
              c35 = c35ext(ii2,jj2,ispec_elastic)
              c55 = c55ext(ii2,jj2,ispec_elastic)
              c12 = c12ext(ii2,jj2,ispec_elastic)
              c23 = c23ext(ii2,jj2,ispec_elastic)
              c25 = c25ext(ii2,jj2,ispec_elastic)
            else
              c11 = anisotropy(1,kmato(ispec_elastic))
              c13 = anisotropy(2,kmato(ispec_elastic))
              c15 = anisotropy(3,kmato(ispec_elastic))
              c33 = anisotropy(4,kmato(ispec_elastic))
              c35 = anisotropy(5,kmato(ispec_elastic))
              c55 = anisotropy(6,kmato(ispec_elastic))
              c12 = anisotropy(7,kmato(ispec_elastic))
              c23 = anisotropy(8,kmato(ispec_elastic))
              c25 = anisotropy(9,kmato(ispec_elastic))
            endif

            sigma_xx = sigma_xx + c11*dux_dxl + c15*(duz_dxl + dux_dzl) + c13*duz_dzl
            sigma_zz = sigma_zz + c13*dux_dxl + c35*(duz_dxl + dux_dzl) + c33*duz_dzl
            sigma_xz = sigma_xz + c15*dux_dxl + c55*(duz_dxl + dux_dzl) + c35*duz_dzl
          else
            ! no attenuation
            sigma_xx = sigma_xx + lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
            sigma_xz = sigma_xz + mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
            sigma_zz = sigma_zz + lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl
          endif

          if(SIMULATION_TYPE == 3) then
            b_sigma_xx = b_sigma_xx + lambdaplus2mu_unrelaxed_elastic*b_dux_dxl + lambdal_unrelaxed_elastic*b_duz_dzl
            b_sigma_xz = b_sigma_xz + mul_unrelaxed_elastic*(b_duz_dxl + b_dux_dzl)
            b_sigma_zz = b_sigma_zz + lambdaplus2mu_unrelaxed_elastic*b_duz_dzl + lambdal_unrelaxed_elastic*b_dux_dxl
          endif

          ! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
          ! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
          ! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
          ! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
          ! Blackwell Science, page 110, equation (4.60).
          if(iedge_poroelastic == ITOP)then
            xxi = + gammaz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zxi = - gammax(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = - zxi / jacobian1D
            nz = + xxi / jacobian1D
            weight = jacobian1D * wxgll(i)
          else if(iedge_poroelastic == IBOTTOM)then
            xxi = + gammaz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zxi = - gammax(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D
            weight = jacobian1D * wxgll(i)
          else if(iedge_poroelastic ==ILEFT)then
            xgamma = - xiz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zgamma = + xix(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = - zgamma / jacobian1D
            nz = + xgamma / jacobian1D
            weight = jacobian1D * wzgll(j)
          else if(iedge_poroelastic ==IRIGHT)then
            xgamma = - xiz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zgamma = + xix(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = + zgamma / jacobian1D
            nz = - xgamma / jacobian1D
            weight = jacobian1D * wzgll(j)
          endif

          accel_elastic(1,iglob) = accel_elastic(1,iglob) - weight* &
                (sigma_xx*nx + sigma_xz*nz)/2.d0

          accel_elastic(3,iglob) = accel_elastic(3,iglob) - weight* &
                (sigma_xz*nx + sigma_zz*nz)/2.d0

          if(SIMULATION_TYPE == 3) then
            b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - weight* &
                (b_sigma_xx*nx + b_sigma_xz*nz)/2.d0

            b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) - weight* &
                (b_sigma_xz*nx + b_sigma_zz*nz)/2.d0
          endif !if(SIMULATION_TYPE == 3) then

        enddo

      enddo

    endif


! ************************************************************************************
! ************************************ add force source
! ************************************************************************************

    if(any_elastic) then

      ! --- add the source if it is a collocated force
      if(.not. initialfield) then

        do i_source=1,NSOURCES
          ! if this processor core carries the source and the source element is elastic
          if (is_proc_source(i_source) == 1 .and. elastic(ispec_selected_source(i_source))) then

            ! collocated force
            if(source_type(i_source) == 1) then
              if(SIMULATION_TYPE == 1) then  ! forward wavefield
                if(p_sv) then ! P-SV calculation
                  do j = 1,NGLLZ
                    do i = 1,NGLLX
                      iglob = ibool(i,j,ispec_selected_source(i_source))
                      hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
                      accel_elastic(1,iglob) = accel_elastic(1,iglob) &
                        - sin(anglesource(i_source))*source_time_function(i_source,it,i_stage)*hlagrange
                      accel_elastic(3,iglob) = accel_elastic(3,iglob) &
                        + cos(anglesource(i_source))*source_time_function(i_source,it,i_stage)*hlagrange
                    enddo
                  enddo
                else    ! SH (membrane) calculation
                  do j = 1,NGLLZ
                    do i = 1,NGLLX
                      iglob = ibool(i,j,ispec_selected_source(i_source))
                      hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
                      accel_elastic(2,iglob) = accel_elastic(2,iglob) &
                            + source_time_function(i_source,it,i_stage)*hlagrange
                    enddo
                  enddo
                endif
              else                   ! backward wavefield
                if(p_sv) then ! P-SV calculation
                  do j = 1,NGLLZ
                    do i = 1,NGLLX
                      iglob = ibool(i,j,ispec_selected_source(i_source))
                      hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
                      b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) &
                        - sin(anglesource(i_source))*source_time_function(i_source,NSTEP-it+1,stage_time_scheme-i_stage+1) &
                          *hlagrange
                      b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) &
                        + cos(anglesource(i_source))*source_time_function(i_source,NSTEP-it+1,stage_time_scheme-i_stage+1) &
                          *hlagrange
                    enddo
                  enddo
                else    ! SH (membrane) calculation
                  do j = 1,NGLLZ
                    do i = 1,NGLLX
                      iglob = ibool(i,j,ispec_selected_source(i_source))
                      hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
                      b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) &
                            + source_time_function(i_source,NSTEP-it+1,stage_time_scheme-i_stage+1)*hlagrange
                    enddo
                  enddo
                endif

              endif  !endif SIMULATION_TYPE == 1
            endif

          endif ! if this processor core carries the source and the source element is elastic
        enddo ! do i_source=1,NSOURCES

!<NOISE_TOMOGRAPHY

        ! inject wavefield sources for noise simulations

        if (NOISE_TOMOGRAPHY == 1) then
          call  add_point_source_noise()

        else if (NOISE_TOMOGRAPHY == 2) then
          call add_surface_movie_noise(accel_elastic)

        else if (NOISE_TOMOGRAPHY == 3) then
          if (.not. save_everywhere) then
            call add_surface_movie_noise(b_accel_elastic)
          endif
        endif

!>NOISE_TOMOGRAPHY


      endif ! if not using an initial field
    endif !if(any_elastic)

! assembling accel_elastic for elastic elements
#ifdef USE_MPI

    if(time_stepping_scheme == 2)then
    if(i_stage==1 .and. it == 1 .and. (.not. initialfield))then
    veloc_elastic_LDDRK_temp = veloc_elastic
      if (nproc > 1 .and. any_elastic .and. ninterface_elastic > 0) then
       call assemble_MPI_vector_el(veloc_elastic)
       endif
    endif
    endif

    call MPI_BARRIER(MPI_COMM_WORLD,ier)

    if (nproc > 1 .and. any_elastic .and. ninterface_elastic > 0) then
      call assemble_MPI_vector_el(accel_elastic)
    endif


    if (nproc > 1 .and. any_elastic .and. ninterface_elastic > 0 .and. SIMULATION_TYPE == 3) then
      call assemble_MPI_vector_el(b_accel_elastic)
    endif
#endif

      if(PML_BOUNDARY_CONDITIONS .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1)then
       if(any_elastic .and. nglob_interface > 0)then
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

      if(SIMULATION_TYPE == 3)then
        if(PML_BOUNDARY_CONDITIONS)then
          if(any_elastic .and. nglob_interface > 0)then
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

    if(any_elastic) then

!! DK DK this should be vectorized
      accel_elastic(1,:) = accel_elastic(1,:) * rmass_inverse_elastic_one
      accel_elastic(2,:) = accel_elastic(2,:) * rmass_inverse_elastic_one
      accel_elastic(3,:) = accel_elastic(3,:) * rmass_inverse_elastic_three

      if(time_stepping_scheme == 1)then
!! DK DK this should be vectorized
        veloc_elastic = veloc_elastic + deltatover2 * accel_elastic
      endif

      if(time_stepping_scheme == 2)then

!! DK DK this should be vectorized
        veloc_elastic_LDDRK = alpha_LDDRK(i_stage) * veloc_elastic_LDDRK + deltat * accel_elastic
        displ_elastic_LDDRK = alpha_LDDRK(i_stage) * displ_elastic_LDDRK + deltat * veloc_elastic
        if(i_stage==1 .and. it == 1 .and. (.not. initialfield))then
        veloc_elastic_LDDRK_temp = veloc_elastic_LDDRK_temp + beta_LDDRK(i_stage) * veloc_elastic_LDDRK
        veloc_elastic = veloc_elastic_LDDRK_temp
        else
        veloc_elastic = veloc_elastic + beta_LDDRK(i_stage) * veloc_elastic_LDDRK
        endif
        displ_elastic = displ_elastic + beta_LDDRK(i_stage) * displ_elastic_LDDRK

      endif

      if(time_stepping_scheme == 3)then

!! DK DK this should be vectorized
        accel_elastic_rk(1,:,i_stage) = deltat * accel_elastic(1,:)
        accel_elastic_rk(2,:,i_stage) = deltat * accel_elastic(2,:)
        accel_elastic_rk(3,:,i_stage) = deltat * accel_elastic(3,:)

        veloc_elastic_rk(1,:,i_stage) = deltat * veloc_elastic(1,:)
        veloc_elastic_rk(2,:,i_stage) = deltat * veloc_elastic(2,:)
        veloc_elastic_rk(3,:,i_stage) = deltat * veloc_elastic(3,:)

        if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then

        if(i_stage == 1)weight_rk = 0.5d0
        if(i_stage == 2)weight_rk = 0.5d0
        if(i_stage == 3)weight_rk = 1.0d0

        if(i_stage==1)then

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

        else if(i_stage==4)then

!! DK DK this should be vectorized
        veloc_elastic(1,:) = veloc_elastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
        (accel_elastic_rk(1,:,1) + 2.0d0 * accel_elastic_rk(1,:,2) + &
         2.0d0 * accel_elastic_rk(1,:,3) + accel_elastic_rk(1,:,4))

        veloc_elastic(2,:) = veloc_elastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
        (accel_elastic_rk(2,:,1) + 2.0d0 * accel_elastic_rk(2,:,2) + &
         2.0d0 * accel_elastic_rk(2,:,3) + accel_elastic_rk(2,:,4))

         veloc_elastic(3,:) = veloc_elastic_initial_rk(3,:) + 1.0d0 / 6.0d0 * &
        (accel_elastic_rk(3,:,1) + 2.0d0 * accel_elastic_rk(3,:,2) + &
         2.0d0 * accel_elastic_rk(3,:,3) + accel_elastic_rk(3,:,4))

        displ_elastic(1,:) = displ_elastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
        (veloc_elastic_rk(1,:,1) + 2.0d0 * veloc_elastic_rk(1,:,2) + &
         2.0d0 * veloc_elastic_rk(1,:,3) + veloc_elastic_rk(1,:,4))

        displ_elastic(2,:) = displ_elastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
        (veloc_elastic_rk(2,:,1) + 2.0d0 * veloc_elastic_rk(2,:,2) + &
         2.0d0 * veloc_elastic_rk(2,:,3) + veloc_elastic_rk(2,:,4))

        displ_elastic(3,:) = displ_elastic_initial_rk(3,:) + 1.0d0 / 6.0d0 * &
        (veloc_elastic_rk(3,:,1) + 2.0d0 * veloc_elastic_rk(3,:,2) + &
         2.0d0 * veloc_elastic_rk(3,:,3) + veloc_elastic_rk(3,:,4))

        endif

      endif

      if(SIMULATION_TYPE == 3) then
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

  if ( .NOT. GPU_MODE) then

    if(any_poroelastic) then

!--------------------------------------------------------------------------------------------
! implement viscous attenuation for poroelastic media
!
    if(ATTENUATION_PORO_FLUID_PART) then
      ! update memory variables with fourth-order Runge-Kutta time scheme for attenuation
      ! loop over spectral elements
      do ispec = 1,nspec

       if(poroelastic(ispec) .and. poroelastcoef(2,2,kmato(ispec)) >0.d0) then

        etal_f = poroelastcoef(2,2,kmato(ispec))
        permlxx = permeability(1,kmato(ispec))
        permlxz = permeability(2,kmato(ispec))
        permlzz = permeability(3,kmato(ispec))

        ! calcul of the inverse of k

        detk = permlxx*permlzz - permlxz*permlxz

        if(detk /= ZERO) then
          invpermlxx = permlzz/detk
          invpermlxz = -permlxz/detk
          invpermlzz = permlxx/detk
        else
          stop 'Permeability matrix is not invertible'
        endif

        ! relaxed viscous coef
        bl_unrelaxed_elastic(1) = etal_f*invpermlxx
        bl_unrelaxed_elastic(2) = etal_f*invpermlxz
        bl_unrelaxed_elastic(3) = etal_f*invpermlzz

        do j=1,NGLLZ
          do i=1,NGLLX

            iglob = ibool(i,j,ispec)
            viscox_loc(i,j) = velocw_poroelastic(1,iglob)*bl_unrelaxed_elastic(1) + &
                               velocw_poroelastic(2,iglob) * bl_unrelaxed_elastic(2)
            viscoz_loc(i,j) = velocw_poroelastic(1,iglob)*bl_unrelaxed_elastic(2) + &
                               velocw_poroelastic(2,iglob)*bl_unrelaxed_elastic(3)

            if(time_stepping_scheme == 1) then
              ! evolution rx_viscous
              Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscox(i,j,ispec)
              Snp1 = - (1.d0 - theta_e/theta_s)/theta_s*viscox_loc(i,j)
              rx_viscous(i,j,ispec) = alphaval * rx_viscous(i,j,ispec) &
                                    + betaval * Sn + gammaval * Snp1

              ! evolution rz_viscous
              Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscoz(i,j,ispec)
              Snp1 = - (1.d0 - theta_e/theta_s)/theta_s*viscoz_loc(i,j)
              rz_viscous(i,j,ispec) = alphaval * rz_viscous(i,j,ispec) &
                                    + betaval * Sn + gammaval * Snp1
            endif

            if(time_stepping_scheme == 2) then
              Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscox(i,j,ispec)
              rx_viscous_LDDRK(i,j,ispec) = alpha_LDDRK(i_stage) * rx_viscous_LDDRK(i,j,ispec) + &
                                            deltat * (Sn + thetainv * rx_viscous(i,j,ispec))
              rx_viscous(i,j,ispec)= rx_viscous(i,j,ispec)+beta_LDDRK(i_stage) * rx_viscous_LDDRK(i,j,ispec)

              Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscoz(i,j,ispec)
              rz_viscous_LDDRK(i,j,ispec)= alpha_LDDRK(i_stage) * rz_viscous_LDDRK(i,j,ispec)+&
                                           deltat * (Sn + thetainv * rz_viscous(i,j,ispec))
              rz_viscous(i,j,ispec)= rz_viscous(i,j,ispec)+beta_LDDRK(i_stage) * rz_viscous_LDDRK(i,j,ispec)
            endif

            if(time_stepping_scheme == 3) then

              Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscox(i,j,ispec)
              rx_viscous_force_RK(i,j,ispec,i_stage) = deltat * (Sn + thetainv * rx_viscous(i,j,ispec))

              if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then
                if(i_stage == 1)weight_rk = 0.5d0
                if(i_stage == 2)weight_rk = 0.5d0
                if(i_stage == 3)weight_rk = 1.0d0

                if(i_stage==1)then
                  rx_viscous_initial_rk(i,j,ispec) = rx_viscous(i,j,ispec)
                endif
                  rx_viscous(i,j,ispec) = rx_viscous_initial_rk(i,j,ispec) + &
                                          weight_rk * rx_viscous_force_RK(i,j,ispec,i_stage)
              else if(i_stage==4)then

                rx_viscous(i,j,ispec) = rx_viscous_initial_rk(i,j,ispec) + &
                                        1.0d0 / 6.0d0 * (rx_viscous_force_RK(i,j,ispec,i_stage) + &
                                        2.0d0 * rx_viscous_force_RK(i,j,ispec,i_stage) + &
                                        2.0d0 * rx_viscous_force_RK(i,j,ispec,i_stage) + &
                                        rx_viscous_force_RK(i,j,ispec,i_stage))
              endif

              Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscoz(i,j,ispec)
              rz_viscous_force_RK(i,j,ispec,i_stage) = deltat * (Sn + thetainv * rz_viscous(i,j,ispec))

              if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then
                if(i_stage == 1)weight_rk = 0.5d0
                if(i_stage == 2)weight_rk = 0.5d0
                if(i_stage == 3)weight_rk = 1.0d0
                if(i_stage==1)then
                  rz_viscous_initial_rk(i,j,ispec) = rz_viscous(i,j,ispec)
                endif
                rz_viscous(i,j,ispec) = rz_viscous_initial_rk(i,j,ispec) + &
                                        weight_rk * rz_viscous_force_RK(i,j,ispec,i_stage)
              else if(i_stage==4)then
                rz_viscous(i,j,ispec) = rz_viscous_initial_rk(i,j,ispec) + &
                                        1.0d0 / 6.0d0 * (rz_viscous_force_RK(i,j,ispec,i_stage) + &
                                        2.0d0 * rz_viscous_force_RK(i,j,ispec,i_stage) + &
                                        2.0d0 * rz_viscous_force_RK(i,j,ispec,i_stage) + &
                                        rz_viscous_force_RK(i,j,ispec,i_stage))
              endif
            endif
          enddo
        enddo

        if(stage_time_scheme == 1) then
        ! save visco for Runge-Kutta scheme when used together with Newmark
        viscox(:,:,ispec) = viscox_loc(:,:)
        viscoz(:,:,ispec) = viscoz_loc(:,:)
        endif

       endif  ! end of poroelastic element loop

      enddo   ! end of spectral element loop
    endif ! end of viscous attenuation for porous media

!-----------------------------------------

      if(SIMULATION_TYPE == 3) then
        ! if inviscid fluid, comment the reading and uncomment the zeroing
        !     read(23,rec=NSTEP-it+1) b_viscodampx
        !     read(24,rec=NSTEP-it+1) b_viscodampz
        b_viscodampx(:) = ZERO
        b_viscodampz(:) = ZERO
      endif

      call compute_forces_poro_solid(f0(1))

      call compute_forces_poro_fluid(f0(1))


      if(SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
        ! if inviscid fluid, comment
        !     write(23,rec=it) b_viscodampx
        !     write(24,rec=it) b_viscodampz
      endif

      if(anyabs .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1) then

        !--- left absorbing boundary
        if(nspec_left >0) then
          do ispec = 1,nspec_left
            do id =1,2
              do i=1,NGLLZ
                write(45) b_absorb_poro_s_left(id,i,ispec,it)
                write(25) b_absorb_poro_w_left(id,i,ispec,it)
              enddo
            enddo
          enddo
        endif

        !--- right absorbing boundary
        if(nspec_right >0) then
          do ispec = 1,nspec_right
            do id =1,2
              do i=1,NGLLZ
                write(46) b_absorb_poro_s_right(id,i,ispec,it)
                write(26) b_absorb_poro_w_right(id,i,ispec,it)
              enddo
            enddo
          enddo
        endif

        !--- bottom absorbing boundary
        if(nspec_bottom >0) then
          do ispec = 1,nspec_bottom
            do id =1,2
              do i=1,NGLLX
                write(47) b_absorb_poro_s_bottom(id,i,ispec,it)
                write(29) b_absorb_poro_w_bottom(id,i,ispec,it)
              enddo
            enddo
          enddo
        endif

        !--- top absorbing boundary
        if(nspec_top >0) then
          do ispec = 1,nspec_top
            do id =1,2
              do i=1,NGLLX
                write(48) b_absorb_poro_s_top(id,i,ispec,it)
                write(28) b_absorb_poro_w_top(id,i,ispec,it)
              enddo
            enddo
          enddo
        endif

      endif ! if(anyabs .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1)

    endif !if(any_poroelastic) then

! *********************************************************
! ************* add coupling with the acoustic side
! *********************************************************

    if(coupled_acoustic_poro) then

      ! loop on all the coupling edges
      do inum = 1,num_fluid_poro_edges

        ! get the edge of the acoustic element
        ispec_acoustic = fluid_poro_acoustic_ispec(inum)
        iedge_acoustic = fluid_poro_acoustic_iedge(inum)

        ! get the corresponding edge of the poroelastic element
        ispec_poroelastic = fluid_poro_poroelastic_ispec(inum)
        iedge_poroelastic = fluid_poro_poroelastic_iedge(inum)

        ! implement 1D coupling along the edge
        do ipoin1D = 1,NGLLX

          ! get point values for the acoustic side, which matches our side in the inverse direction
          i = ivalue_inverse(ipoin1D,iedge_acoustic)
          j = jvalue_inverse(ipoin1D,iedge_acoustic)
          iglob = ibool(i,j,ispec_acoustic)

          ! get poroelastic parameters
          phil = porosity(kmato(ispec_poroelastic))
          tortl = tortuosity(kmato(ispec_poroelastic))
          rhol_f = density(2,kmato(ispec_poroelastic))
          rhol_s = density(1,kmato(ispec_poroelastic))
          rhol_bar = (1._CUSTOM_REAL-phil)*rhol_s + phil*rhol_f

          ! compute pressure on the fluid/porous medium edge
          pressure = - potential_dot_dot_acoustic(iglob)
          if(SIMULATION_TYPE == 3) then
            b_pressure = - b_potential_dot_dot_acoustic(iglob)
            ! new definition of adjoint displacement and adjoint potential
            pressure = potential_acoustic_adj_coupling(iglob)
          endif

          ! get point values for the poroelastic side
          ii2 = ivalue(ipoin1D,iedge_poroelastic)
          jj2 = jvalue(ipoin1D,iedge_poroelastic)
          iglob = ibool(ii2,jj2,ispec_poroelastic)

          ! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
          ! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
          ! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
          ! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
          ! Blackwell Science, page 110, equation (4.60).
          if(iedge_acoustic == ITOP)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = - zxi / jacobian1D
            nz = + xxi / jacobian1D
            weight = jacobian1D * wxgll(i)
          else if(iedge_acoustic == IBOTTOM)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D
            weight = jacobian1D * wxgll(i)
          else if(iedge_acoustic ==ILEFT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = - zgamma / jacobian1D
            nz = + xgamma / jacobian1D
            weight = jacobian1D * wzgll(j)
          else if(iedge_acoustic ==IRIGHT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = + zgamma / jacobian1D
            nz = - xgamma / jacobian1D
            weight = jacobian1D * wzgll(j)
          endif

          ! contribution to the solid phase
          accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) &
            + weight*nx*pressure*(1._CUSTOM_REAL-phil/tortl)
          accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) &
            + weight*nz*pressure*(1._CUSTOM_REAL-phil/tortl)

          ! contribution to the fluid phase
          accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) &
            + weight*nx*pressure*(1._CUSTOM_REAL-rhol_f/rhol_bar)
          accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) &
            + weight*nz*pressure*(1._CUSTOM_REAL-rhol_f/rhol_bar)

          if(SIMULATION_TYPE == 3) then
            ! contribution to the solid phase
            b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) &
              + weight*nx*b_pressure*(1._CUSTOM_REAL-phil/tortl)
            b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) &
              + weight*nz*b_pressure*(1._CUSTOM_REAL-phil/tortl)

            ! contribution to the fluid phase
            b_accelw_poroelastic(1,iglob) = b_accelw_poroelastic(1,iglob) &
              + weight*nx*b_pressure*(1._CUSTOM_REAL-rhol_f/rhol_bar)
            b_accelw_poroelastic(2,iglob) = b_accelw_poroelastic(2,iglob) &
              + weight*nz*b_pressure*(1._CUSTOM_REAL-rhol_f/rhol_bar)
          endif !if(SIMULATION_TYPE == 3) then

        enddo ! do ipoin1D = 1,NGLLX

      enddo ! do inum = 1,num_fluid_poro_edges

    endif ! if(coupled_acoustic_poro)

! ****************************************************************************
! ************* add coupling with the elastic side
! ****************************************************************************

    if(coupled_elastic_poro) then

      ! loop on all the coupling edges
      do inum = 1,num_solid_poro_edges

        ! get the edge of the elastic element
        ispec_elastic = solid_poro_elastic_ispec(inum)
        iedge_elastic = solid_poro_elastic_iedge(inum)

        ! get the corresponding edge of the poroelastic element
        ispec_poroelastic = solid_poro_poroelastic_ispec(inum)
        iedge_poroelastic = solid_poro_poroelastic_iedge(inum)

        ! implement 1D coupling along the edge
        do ipoin1D = 1,NGLLX

          ! get point values for the elastic side, which matches our side in the inverse direction
          i = ivalue_inverse(ipoin1D,iedge_elastic)
          j = jvalue_inverse(ipoin1D,iedge_elastic)
          iglob = ibool(i,j,ispec_elastic)

          ! get elastic properties
          lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec_elastic))
          mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec_elastic))
          lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec_elastic))

          ! derivative along x and along z for u_s and w
          dux_dxi = ZERO
          duz_dxi = ZERO

          dux_dgamma = ZERO
          duz_dgamma = ZERO

          if(SIMULATION_TYPE == 3) then
            b_dux_dxi = ZERO
            b_duz_dxi = ZERO

            b_dux_dgamma = ZERO
            b_duz_dgamma = ZERO
          endif

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec_elastic))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j,ispec_elastic))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec_elastic))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i,k,ispec_elastic))*hprime_zz(j,k)

            if(SIMULATION_TYPE == 3) then
              b_dux_dxi = b_dux_dxi + b_displ_elastic(1,ibool(k,j,ispec_elastic))*hprime_xx(i,k)
              b_duz_dxi = b_duz_dxi + b_displ_elastic(3,ibool(k,j,ispec_elastic))*hprime_xx(i,k)
              b_dux_dgamma = b_dux_dgamma + b_displ_elastic(1,ibool(i,k,ispec_elastic))*hprime_zz(j,k)
              b_duz_dgamma = b_duz_dgamma + b_displ_elastic(3,ibool(i,k,ispec_elastic))*hprime_zz(j,k)
            endif
          enddo

          xixl = xix(i,j,ispec_elastic)
          xizl = xiz(i,j,ispec_elastic)
          gammaxl = gammax(i,j,ispec_elastic)
          gammazl = gammaz(i,j,ispec_elastic)

          ! derivatives of displacement
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

          if(SIMULATION_TYPE == 3) then
            b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
            b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl

            b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
            b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl
          endif
          ! compute stress tensor
          ! full anisotropy
          if(kmato(ispec_elastic) == 2) then
            ! implement anisotropy in 2D
            if(assign_external_model) then
              c11 = c11ext(i,j,ispec_elastic)
              c13 = c13ext(i,j,ispec_elastic)
              c15 = c15ext(i,j,ispec_elastic)
              c33 = c33ext(i,j,ispec_elastic)
              c35 = c35ext(i,j,ispec_elastic)
              c55 = c55ext(i,j,ispec_elastic)
              c12 = c12ext(i,j,ispec_elastic)
              c23 = c23ext(i,j,ispec_elastic)
              c25 = c25ext(i,j,ispec_elastic)
            else
              c11 = anisotropy(1,kmato(ispec_elastic))
              c13 = anisotropy(2,kmato(ispec_elastic))
              c15 = anisotropy(3,kmato(ispec_elastic))
              c33 = anisotropy(4,kmato(ispec_elastic))
              c35 = anisotropy(5,kmato(ispec_elastic))
              c55 = anisotropy(6,kmato(ispec_elastic))
              c12 = anisotropy(7,kmato(ispec_elastic))
              c23 = anisotropy(8,kmato(ispec_elastic))
              c25 = anisotropy(9,kmato(ispec_elastic))
            endif
            sigma_xx = c11*dux_dxl + c15*(duz_dxl + dux_dzl) + c13*duz_dzl
            sigma_zz = c13*dux_dxl + c35*(duz_dxl + dux_dzl) + c33*duz_dzl
            sigma_xz = c15*dux_dxl + c55*(duz_dxl + dux_dzl) + c35*duz_dzl
          else
            ! no attenuation
            sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
            sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
            sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl
          endif

          if(SIMULATION_TYPE == 3) then
            b_sigma_xx = lambdaplus2mu_unrelaxed_elastic*b_dux_dxl + lambdal_unrelaxed_elastic*b_duz_dzl
            b_sigma_xz = mul_unrelaxed_elastic*(b_duz_dxl + b_dux_dzl)
            b_sigma_zz = lambdaplus2mu_unrelaxed_elastic*b_duz_dzl + lambdal_unrelaxed_elastic*b_dux_dxl
          endif ! if(SIMULATION_TYPE == 3)

          ! get point values for the poroelastic side
          i = ivalue(ipoin1D,iedge_poroelastic)
          j = jvalue(ipoin1D,iedge_poroelastic)
          iglob = ibool(i,j,ispec_poroelastic)

          ! get poroelastic domain paramters
          phil = porosity(kmato(ispec_poroelastic))
          tortl = tortuosity(kmato(ispec_poroelastic))
          !solid properties
          mul_s = poroelastcoef(2,1,kmato(ispec_poroelastic))
          kappal_s = poroelastcoef(3,1,kmato(ispec_poroelastic)) - 4._CUSTOM_REAL*mul_s/3._CUSTOM_REAL
          rhol_s = density(1,kmato(ispec_poroelastic))
          !fluid properties
          kappal_f = poroelastcoef(1,2,kmato(ispec_poroelastic))
          rhol_f = density(2,kmato(ispec_poroelastic))
          !frame properties
          mul_fr = poroelastcoef(2,3,kmato(ispec_poroelastic))
          kappal_fr = poroelastcoef(3,3,kmato(ispec_poroelastic)) - 4._CUSTOM_REAL*mul_fr/3._CUSTOM_REAL
          rhol_bar =  (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f
          !Biot coefficients for the input phi
          D_biot = kappal_s*(1._CUSTOM_REAL + phil*(kappal_s/kappal_f - 1._CUSTOM_REAL))
          H_biot = (kappal_s - kappal_fr)*(kappal_s - kappal_fr)/(D_biot - kappal_fr) + &
                kappal_fr + 4._CUSTOM_REAL*mul_fr/3._CUSTOM_REAL
          C_biot = kappal_s*(kappal_s - kappal_fr)/(D_biot - kappal_fr)
          M_biot = kappal_s*kappal_s/(D_biot - kappal_fr)
          mul_G = mul_fr
          lambdal_G = H_biot - 2._CUSTOM_REAL*mul_fr
          lambdalplus2mul_G = lambdal_G + TWO*mul_G

          ! derivative along x and along z for u_s and w
          dux_dxi = ZERO
          duz_dxi = ZERO

          dux_dgamma = ZERO
          duz_dgamma = ZERO

          dwx_dxi = ZERO
          dwz_dxi = ZERO

          dwx_dgamma = ZERO
          dwz_dgamma = ZERO

          if(SIMULATION_TYPE == 3) then
            b_dux_dxi = ZERO
            b_duz_dxi = ZERO

            b_dux_dgamma = ZERO
            b_duz_dgamma = ZERO

            b_dwx_dxi = ZERO
            b_dwz_dxi = ZERO

            b_dwx_dgamma = ZERO
            b_dwz_dgamma = ZERO
          endif

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displs_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displs_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displs_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displs_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)

            dwx_dxi = dwx_dxi + displw_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
            dwz_dxi = dwz_dxi + displw_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
            dwx_dgamma = dwx_dgamma + displw_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
            dwz_dgamma = dwz_dgamma + displw_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
            if(SIMULATION_TYPE == 3) then
              b_dux_dxi = b_dux_dxi + b_displs_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
              b_duz_dxi = b_duz_dxi + b_displs_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
              b_dux_dgamma = b_dux_dgamma + b_displs_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
              b_duz_dgamma = b_duz_dgamma + b_displs_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)

              b_dwx_dxi = b_dwx_dxi + b_displw_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
              b_dwz_dxi = b_dwz_dxi + b_displw_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
              b_dwx_dgamma = b_dwx_dgamma + b_displw_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
              b_dwz_dgamma = b_dwz_dgamma + b_displw_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
            endif
          enddo

          xixl = xix(i,j,ispec_poroelastic)
          xizl = xiz(i,j,ispec_poroelastic)
          gammaxl = gammax(i,j,ispec_poroelastic)
          gammazl = gammaz(i,j,ispec_poroelastic)

          ! derivatives of displacement
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

          dwx_dxl = dwx_dxi*xixl + dwx_dgamma*gammaxl
          dwx_dzl = dwx_dxi*xizl + dwx_dgamma*gammazl

          dwz_dxl = dwz_dxi*xixl + dwz_dgamma*gammaxl
          dwz_dzl = dwz_dxi*xizl + dwz_dgamma*gammazl

          if(SIMULATION_TYPE == 3) then
            b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
            b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl

            b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
            b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl

            b_dwx_dxl = b_dwx_dxi*xixl + b_dwx_dgamma*gammaxl
            b_dwx_dzl = b_dwx_dxi*xizl + b_dwx_dgamma*gammazl

            b_dwz_dxl = b_dwz_dxi*xixl + b_dwz_dgamma*gammaxl
            b_dwz_dzl = b_dwz_dxi*xizl + b_dwz_dgamma*gammazl
          endif
          ! compute stress tensor

          ! no attenuation
          sigma_xx = sigma_xx + lambdalplus2mul_G*dux_dxl + lambdal_G*duz_dzl + C_biot*(dwx_dxl + dwz_dzl)
          sigma_xz = sigma_xz + mul_G*(duz_dxl + dux_dzl)
          sigma_zz = sigma_zz + lambdalplus2mul_G*duz_dzl + lambdal_G*dux_dxl + C_biot*(dwx_dxl + dwz_dzl)

          sigmap = C_biot*(dux_dxl + duz_dzl) + M_biot*(dwx_dxl + dwz_dzl)

          if(SIMULATION_TYPE == 3) then
            b_sigma_xx = b_sigma_xx + lambdalplus2mul_G*b_dux_dxl + lambdal_G*b_duz_dzl + C_biot*(b_dwx_dxl + b_dwz_dzl)
            b_sigma_xz = b_sigma_xz + mul_G*(b_duz_dxl + b_dux_dzl)
            b_sigma_zz = b_sigma_zz + lambdalplus2mul_G*b_duz_dzl + lambdal_G*b_dux_dxl + C_biot*(b_dwx_dxl + b_dwz_dzl)
            b_sigmap = C_biot*(b_dux_dxl + b_duz_dzl) + M_biot*(b_dwx_dxl + b_dwz_dzl)
          endif

          ! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
          ! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
          ! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
          ! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
          ! Blackwell Science, page 110, equation (4.60).
          if(iedge_poroelastic == ITOP)then
            xxi = + gammaz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zxi = - gammax(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = - zxi / jacobian1D
            nz = + xxi / jacobian1D
            weight = jacobian1D * wxgll(i)
          else if(iedge_poroelastic == IBOTTOM)then
            xxi = + gammaz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zxi = - gammax(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D
            weight = jacobian1D * wxgll(i)
          else if(iedge_poroelastic ==ILEFT)then
            xgamma = - xiz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zgamma = + xix(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = - zgamma / jacobian1D
            nz = + xgamma / jacobian1D
            weight = jacobian1D * wzgll(j)
          else if(iedge_poroelastic ==IRIGHT)then
            xgamma = - xiz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zgamma = + xix(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = + zgamma / jacobian1D
            nz = - xgamma / jacobian1D
            weight = jacobian1D * wzgll(j)
          endif

          ! contribution to the solid phase
          accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) + &
                weight*((sigma_xx*nx + sigma_xz*nz)/2.d0 -phil/tortl*sigmap*nx)

          accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) + &
                weight*((sigma_xz*nx + sigma_zz*nz)/2.d0 -phil/tortl*sigmap*nz)

          ! contribution to the fluid phase
          ! w = 0

          if(SIMULATION_TYPE == 3) then
            ! contribution to the solid phase
            b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) + &
                weight*((b_sigma_xx*nx + b_sigma_xz*nz)/2.d0 -phil/tortl*b_sigmap*nx)

            b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) + &
                weight*((b_sigma_xz*nx + b_sigma_zz*nz)/2.d0 -phil/tortl*b_sigmap*nz)

            ! contribution to the fluid phase
            ! w = 0
          endif !if(SIMULATION_TYPE == 3) then

        enddo

      enddo

    endif ! if(coupled_elastic_poro)


! ************************************************************************************
! ******************************** add force source
! ************************************************************************************

    if(any_poroelastic) then


      ! --- add the source if it is a collocated force
      if(.not. initialfield) then

        do i_source=1,NSOURCES
          ! if this processor core carries the source and the source element is elastic
          if (is_proc_source(i_source) == 1 .and. poroelastic(ispec_selected_source(i_source))) then

            phil = porosity(kmato(ispec_selected_source(i_source)))
            tortl = tortuosity(kmato(ispec_selected_source(i_source)))
            rhol_s = density(1,kmato(ispec_selected_source(i_source)))
            rhol_f = density(2,kmato(ispec_selected_source(i_source)))
            rhol_bar = (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f

            ! collocated force
            if(source_type(i_source) == 1) then
              if(SIMULATION_TYPE == 1) then  ! forward wavefield
                do j = 1,NGLLZ
                  do i = 1,NGLLX
                    iglob = ibool(i,j,ispec_selected_source(i_source))
                    hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
                    ! s
                    accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) - hlagrange * &
                      (1._CUSTOM_REAL - phil/tortl)*sin(anglesource(i_source))*source_time_function(i_source,it,i_stage)
                    accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) + hlagrange * &
                      (1._CUSTOM_REAL - phil/tortl)*cos(anglesource(i_source))*source_time_function(i_source,it,i_stage)
                    ! w
                    accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) - hlagrange * &
                      (1._CUSTOM_REAL - rhol_f/rhol_bar)*sin(anglesource(i_source))*source_time_function(i_source,it,i_stage)
                    accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) + hlagrange * &
                      (1._CUSTOM_REAL - rhol_f/rhol_bar)*cos(anglesource(i_source))*source_time_function(i_source,it,i_stage)
                  enddo
                enddo
              else                   ! backward wavefield
                do j = 1,NGLLZ
                  do i = 1,NGLLX
                    iglob = ibool(i,j,ispec_selected_source(i_source))
                    hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
                    ! b_s
                    b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) - hlagrange * &
                      (1._CUSTOM_REAL - phil/tortl)*sin(anglesource(i_source))* &
                      source_time_function(i_source,NSTEP-it+1,stage_time_scheme-i_stage+1)
                    b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) + hlagrange * &
                      (1._CUSTOM_REAL - phil/tortl)*cos(anglesource(i_source))* &
                      source_time_function(i_source,NSTEP-it+1,stage_time_scheme-i_stage+1)
                    !b_w
                    b_accelw_poroelastic(1,iglob) = b_accelw_poroelastic(1,iglob) - hlagrange * &
                      (1._CUSTOM_REAL - rhol_f/rhol_bar)*sin(anglesource(i_source))* &
                      source_time_function(i_source,NSTEP-it+1,stage_time_scheme-i_stage+1)
                    b_accelw_poroelastic(2,iglob) = b_accelw_poroelastic(2,iglob) + hlagrange * &
                      (1._CUSTOM_REAL - rhol_f/rhol_bar)*cos(anglesource(i_source))* &
                      source_time_function(i_source,NSTEP-it+1,stage_time_scheme-i_stage+1)
                  enddo
                enddo
              endif !endif SIMULATION_TYPE == 1
            endif

          endif ! if this processor core carries the source and the source element is elastic
        enddo ! do i_source=1,NSOURCES

      endif ! if not using an initial field
    endif !if(any_poroelastic)

! assembling accels_proelastic & accelw_poroelastic for poroelastic elements
#ifdef USE_MPI
    if (nproc > 1 .and. any_poroelastic .and. ninterface_poroelastic > 0) then
      call assemble_MPI_vector_po(accels_poroelastic,accelw_poroelastic)
    endif

    if (nproc > 1 .and. any_poroelastic .and. ninterface_poroelastic > 0 .and. SIMULATION_TYPE == 3) then
      call assemble_MPI_vector_po(b_accels_poroelastic,b_accelw_poroelastic)
    endif
#endif


! ************************************************************************************
! ************* multiply by the inverse of the mass matrix and update velocity
! ************************************************************************************

    if(any_poroelastic) then

      if(time_stepping_scheme == 1)then

      accels_poroelastic(1,:) = accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
      accels_poroelastic(2,:) = accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)
      velocs_poroelastic = velocs_poroelastic + deltatover2*accels_poroelastic

      accelw_poroelastic(1,:) = accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
      accelw_poroelastic(2,:) = accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)
      velocw_poroelastic = velocw_poroelastic + deltatover2*accelw_poroelastic

      endif

      if(time_stepping_scheme == 2)then

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(time_stepping_scheme == 3)then

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

        if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then

        if(i_stage == 1)weight_rk = 0.5d0
        if(i_stage == 2)weight_rk = 0.5d0
        if(i_stage == 3)weight_rk = 1.0d0

        if(i_stage==1)then

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


        else if(i_stage==4)then

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(SIMULATION_TYPE == 3) then
        b_accels_poroelastic(1,:) = b_accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
        b_accels_poroelastic(2,:) = b_accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)
        b_velocs_poroelastic = b_velocs_poroelastic + b_deltatover2*b_accels_poroelastic

        b_accelw_poroelastic(1,:) = b_accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
        b_accelw_poroelastic(2,:) = b_accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)
        b_velocw_poroelastic = b_velocw_poroelastic + b_deltatover2*b_accelw_poroelastic
      endif

    endif !if(any_poroelastic)

!*******************************************************************************
!         assembling the displacements on the elastic-poro boundaries
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

    if(coupled_elastic_poro) then
      icount(:)=ZERO

      ! loop on all the coupling edges
      do inum = 1,num_solid_poro_edges
        ! get the edge of the elastic element
        ispec_elastic = solid_poro_elastic_ispec(inum)
        iedge_elastic = solid_poro_elastic_iedge(inum)
        ! get the corresponding edge of the poroelastic element
        ispec_poroelastic = solid_poro_poroelastic_ispec(inum)
        iedge_poroelastic = solid_poro_poroelastic_iedge(inum)

        do ipoin1D = 1,NGLLX
          ! recovering original velocities and accelerations on boundaries (elastic side)
          i = ivalue(ipoin1D,iedge_poroelastic)
          j = jvalue(ipoin1D,iedge_poroelastic)
          iglob = ibool(i,j,ispec_poroelastic)
          icount(iglob) = icount(iglob) + 1

          if(icount(iglob) ==1)then

            if(time_stepping_scheme == 1)then

            veloc_elastic(1,iglob) = veloc_elastic(1,iglob) - deltatover2*accel_elastic(1,iglob)
            veloc_elastic(3,iglob) = veloc_elastic(3,iglob) - deltatover2*accel_elastic(3,iglob)
            accel_elastic(1,iglob) = accel_elastic(1,iglob) / rmass_inverse_elastic_one(iglob)
            accel_elastic(3,iglob) = accel_elastic(3,iglob) / rmass_inverse_elastic_three(iglob)
            ! recovering original velocities and accelerations on boundaries (poro side)
            velocs_poroelastic(1,iglob) = velocs_poroelastic(1,iglob) - deltatover2*accels_poroelastic(1,iglob)
            velocs_poroelastic(2,iglob) = velocs_poroelastic(2,iglob) - deltatover2*accels_poroelastic(2,iglob)
            accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) / rmass_s_inverse_poroelastic(iglob)
            accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) / rmass_s_inverse_poroelastic(iglob)
            ! assembling accelerations
            accel_elastic(1,iglob) = ( accel_elastic(1,iglob) + accels_poroelastic(1,iglob) ) / &
                                   ( 1.0/rmass_inverse_elastic_one(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
            accel_elastic(3,iglob) = ( accel_elastic(3,iglob) + accels_poroelastic(2,iglob) ) / &
                                   ( 1.0/rmass_inverse_elastic_three(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
            accels_poroelastic(1,iglob) = accel_elastic(1,iglob)
            accels_poroelastic(2,iglob) = accel_elastic(3,iglob)
            ! updating velocities
            velocs_poroelastic(1,iglob) = velocs_poroelastic(1,iglob) + deltatover2*accels_poroelastic(1,iglob)
            velocs_poroelastic(2,iglob) = velocs_poroelastic(2,iglob) + deltatover2*accels_poroelastic(2,iglob)
            veloc_elastic(1,iglob) = veloc_elastic(1,iglob) + deltatover2*accel_elastic(1,iglob)
            veloc_elastic(3,iglob) = veloc_elastic(3,iglob) + deltatover2*accel_elastic(3,iglob)
            ! zeros w
            accelw_poroelastic(1,iglob) = ZERO
            accelw_poroelastic(2,iglob) = ZERO
            velocw_poroelastic(1,iglob) = ZERO
            velocw_poroelastic(2,iglob) = ZERO

            endif

!            if(time_stepping_scheme == 2)then
            ! recovering original velocities and accelerations on boundaries (elastic side)
!      veloc_elastic = veloc_elastic - beta_LDDRK(i_stage) * veloc_elastic_LDDRK
!      displ_elastic = displ_elastic - beta_LDDRK(i_stage) * displ_elastic_LDDRK
!      veloc_elastic_LDDRK = (veloc_elastic_LDDRK - deltat * accel_elastic) / alpha_LDDRK(i_stage)
!      displ_elastic_LDDRK = (displ_elastic_LDDRK - deltat * veloc_elastic) / alpha_LDDRK(i_stage)
!            accel_elastic(1,iglob) = accel_elastic(1,iglob) / rmass_inverse_elastic(iglob)
!            accel_elastic(3,iglob) = accel_elastic(3,iglob) / rmass_inverse_elastic(iglob)

            ! recovering original velocities and accelerations on boundaries (poro side)
!      velocs_poroelastic = velocs_poroelastic - beta_LDDRK(i_stage) * velocs_poroelastic_LDDRK
!      displs_poroelastic = displs_poroelastic - beta_LDDRK(i_stage) * displs_poroelastic_LDDRK
!      velocs_poroelastic_LDDRK = (velocs_poroelastic_LDDRK - deltat * accels_poroelastic) / alpha_LDDRK(i_stage)
!      displs_poroelastic_LDDRK = (velocs_poroelastic_LDDRK - deltat * velocs_poroelastic) / alpha_LDDRK(i_stage)
!            accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) / rmass_s_inverse_poroelastic(iglob)
!            accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) / rmass_s_inverse_poroelastic(iglob)

            ! assembling accelerations
!            accel_elastic(1,iglob) = ( accel_elastic(1,iglob) + accels_poroelastic(1,iglob) ) / &
!                                   ( 1.0/rmass_inverse_elastic(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
!            accel_elastic(3,iglob) = ( accel_elastic(3,iglob) + accels_poroelastic(2,iglob) ) / &
!                                   ( 1.0/rmass_inverse_elastic(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
!            accels_poroelastic(1,iglob) = accel_elastic(1,iglob)
!            accels_poroelastic(2,iglob) = accel_elastic(3,iglob)

      ! updating velocities
            ! updating velocities(elastic side)
!      veloc_elastic_LDDRK = alpha_LDDRK(i_stage) * veloc_elastic_LDDRK + deltat * accel_elastic
!      displ_elastic_LDDRK = alpha_LDDRK(i_stage) * displ_elastic_LDDRK + deltat * veloc_elastic
!      veloc_elastic = veloc_elastic + beta_LDDRK(i_stage) * veloc_elastic_LDDRK
!      displ_elastic = displ_elastic + beta_LDDRK(i_stage) * displ_elastic_LDDRK
            ! updating velocities(poro side)
!      velocs_poroelastic_LDDRK = alpha_LDDRK(i_stage) * velocs_poroelastic_LDDRK + deltat * accels_poroelastic
!      displs_poroelastic_LDDRK = alpha_LDDRK(i_stage) * displs_poroelastic_LDDRK + deltat * velocs_poroelastic
!      velocs_poroelastic = velocs_poroelastic + beta_LDDRK(i_stage) * velocs_poroelastic_LDDRK
!      displs_poroelastic = displs_poroelastic + beta_LDDRK(i_stage) * displs_poroelastic_LDDRK

            ! zeros w
!            accelw_poroelastic(1,iglob) = ZERO
!            accelw_poroelastic(2,iglob) = ZERO
!            velocw_poroelastic(1,iglob) = ZERO
!            velocw_poroelastic(2,iglob) = ZERO
!            endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      if(time_stepping_scheme == 3)then

        ! recovering original velocities and accelerations on boundaries (elastic side)
!        if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then

!        if(i_stage == 1)weight_rk = 0.5d0
!        if(i_stage == 2)weight_rk = 0.5d0
!        if(i_stage == 3)weight_rk = 1.0d0

!        veloc_elastic(1,iglob) = veloc_elastic_initial_rk(1,iglob) - weight_rk * accel_elastic_rk(1,iglob,i_stage)
!  veloc_elastic(3,iglob) = veloc_elastic_initial_rk(3,iglob) - weight_rk * accel_elastic_rk(3,iglob,i_stage)
!        displ_elastic(1,iglob) = displ_elastic_initial_rk(1,iglob) - weight_rk * veloc_elastic_rk(1,iglob,i_stage)
!  displ_elastic(3,iglob) = displ_elastic_initial_rk(3,iglob) - weight_rk * veloc_elastic_rk(3,iglob,i_stage)


!        else if(i_stage==4)then

!        veloc_elastic(1,iglob) = veloc_elastic_initial_rk(1,iglob) - 1.0d0 / 6.0d0 * &
!        (accel_elastic_rk(1,iglob,1) + 2.0d0 * accel_elastic_rk(1,iglob,2) + &
!         2.0d0 * accel_elastic_rk(1,iglob,3) + accel_elastic_rk(1,iglob,4))

!        veloc_elastic(3,iglob) = veloc_elastic_initial_rk(3,iglob) - 1.0d0 / 6.0d0 * &
!        (accel_elastic_rk(3,iglob,1) + 2.0d0 * accel_elastic_rk(3,iglob,2) + &
!         2.0d0 * accel_elastic_rk(3,iglob,3) + accel_elastic_rk(3,iglob,4))

!        displ_elastic(1,iglob) = displ_elastic_initial_rk(1,iglob) - 1.0d0 / 6.0d0 * &
!        (veloc_elastic_rk(1,iglob,1) + 2.0d0 * veloc_elastic_rk(1,iglob,2) + &
!         2.0d0 * veloc_elastic_rk(1,iglob,3) + veloc_elastic_rk(1,iglob,4))

!        displ_elastic(3,iglob) = displ_elastic_initial_rk(3,iglob) - 1.0d0 / 6.0d0 * &
!        (veloc_elastic_rk(3,iglob,1) + 2.0d0 * veloc_elastic_rk(3,iglob,2) + &
!         2.0d0 * veloc_elastic_rk(3,iglob,3) + veloc_elastic_rk(3,iglob,4))

!        endif

!        accel_elastic(1,iglob) = accel_elastic(1,iglob) / rmass_inverse_elastic(iglob)
!        accel_elastic(3,iglob) = accel_elastic(3,iglob) / rmass_inverse_elastic(iglob)

!        accel_elastic_rk(1,iglob,i_stage) = accel_elastic(1,iglob) / deltat
!        accel_elastic_rk(3,iglob,i_stage) = accel_elastic(3,iglob) / deltat
!        veloc_elastic_rk(1,iglob,i_stage) = veloc_elastic(1,iglob) / deltat
!        veloc_elastic_rk(3,iglob,i_stage) = veloc_elastic(3,iglob) / deltat


        ! recovering original velocities and accelerations on boundaries (poro side)
!        if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then

!        if(i_stage == 1)weight_rk = 0.5d0
!        if(i_stage == 2)weight_rk = 0.5d0
!        if(i_stage == 3)weight_rk = 1.0d0

!        velocs_poroelastic(1,iglob) = velocs_poroelastic_initial_rk(1,iglob) - weight_rk * accels_poroelastic_rk(1,iglob,i_stage)
!  velocs_poroelastic(2,iglob) = velocs_poroelastic_initial_rk(2,iglob) - weight_rk * accels_poroelastic_rk(2,iglob,i_stage)
!        displs_poroelastic(1,iglob) = displs_poroelastic_initial_rk(1,iglob) - weight_rk * velocs_poroelastic_rk(1,iglob,i_stage)
!  displs_poroelastic(2,iglob) = displs_poroelastic_initial_rk(2,iglob) - weight_rk * velocs_poroelastic_rk(2,iglob,i_stage)


!        else if(i_stage==4)then

!        velocs_poroelastic(1,iglob) = velocs_poroelastic_initial_rk(1,iglob) - 1.0d0 / 6.0d0 * &
!        (accels_poroelastic_rk(1,iglob,1) + 2.0d0 * accels_poroelastic_rk(1,iglob,2) + &
!         2.0d0 * accels_poroelastic_rk(1,iglob,3) + accels_poroelastic_rk(1,iglob,4))

!        velocs_poroelastic(2,iglob) = velocs_poroelastic_initial_rk(2,iglob) - 1.0d0 / 6.0d0 * &
!        (accels_poroelastic_rk(2,iglob,1) + 2.0d0 * accels_poroelastic_rk(2,iglob,2) + &
!         2.0d0 * accels_poroelastic_rk(2,iglob,3) + accels_poroelastic_rk(2,iglob,4))

!        displs_poroelastic(1,iglob) = displs_poroelastic_initial_rk(1,iglob) - 1.0d0 / 6.0d0 * &
!        (velocs_poroelastic_rk(1,iglob,1) + 2.0d0 * velocs_poroelastic_rk(1,iglob,2) + &
!         2.0d0 * velocs_poroelastic_rk(1,iglob,3) + velocs_poroelastic_rk(1,iglob,4))

!        displs_poroelastic(2,iglob) = displs_poroelastic_initial_rk(2,iglob) - 1.0d0 / 6.0d0 * &
!        (velocs_poroelastic_rk(2,iglob,1) + 2.0d0 * velocs_poroelastic_rk(2,iglob,2) + &
!         2.0d0 * velocs_poroelastic_rk(2,iglob,3) + velocs_poroelastic_rk(2,iglob,4))

!        endif

!        accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) / rmass_s_inverse_poroelastic(iglob)
!        accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) / rmass_s_inverse_poroelastic(iglob)

!        accels_poroelastic_rk(1,iglob,i_stage) = accels_poroelastic(1,iglob) / deltat
!        accels_poroelastic_rk(2,iglob,i_stage) = accels_poroelastic(2,iglob) / deltat
!        velocs_poroelastic_rk(1,iglob,i_stage) = velocs_poroelastic(1,iglob) / deltat
!        velocs_poroelastic_rk(2,iglob,i_stage) = velocs_poroelastic(2,iglob) / deltat


        ! assembling accelerations
!            accel_elastic(1,iglob) = ( accel_elastic(1,iglob) + accels_poroelastic(1,iglob) ) / &
!                                   ( 1.0/rmass_inverse_elastic(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
!            accel_elastic(3,iglob) = ( accel_elastic(3,iglob) + accels_poroelastic(2,iglob) ) / &
!                                   ( 1.0/rmass_inverse_elastic(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
!            accels_poroelastic(1,iglob) = accel_elastic(1,iglob)
!            accels_poroelastic(2,iglob) = accel_elastic(3,iglob)

   ! updating velocities
        ! updating velocities(elastic side)

 !       accel_elastic_rk(1,iglob,i_stage) = accel_elastic(1,iglob) * deltat
 !       accel_elastic_rk(3,iglob,i_stage) = accel_elastic(3,iglob) * deltat

 !       if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then

 !       if(i_stage == 1)weight_rk = 0.5d0
 !       if(i_stage == 2)weight_rk = 0.5d0
 !       if(i_stage == 3)weight_rk = 1.0d0

 !       veloc_elastic(1,iglob) = veloc_elastic_initial_rk(1,iglob) + weight_rk * accel_elastic_rk(1,iglob,i_stage)
 ! veloc_elastic(3,iglob) = veloc_elastic_initial_rk(3,iglob) + weight_rk * accel_elastic_rk(3,iglob,i_stage)
 !       displ_elastic(1,iglob) = displ_elastic_initial_rk(1,iglob) + weight_rk * veloc_elastic_rk(1,iglob,i_stage)
 ! displ_elastic(3,iglob) = displ_elastic_initial_rk(3,iglob) + weight_rk * veloc_elastic_rk(3,iglob,i_stage)


 !       else if(i_stage==4)then

 !       veloc_elastic(1,iglob) = veloc_elastic_initial_rk(1,iglob) + 1.0d0 / 6.0d0 * &
 !       (accel_elastic_rk(1,iglob,1) + 2.0d0 * accel_elastic_rk(1,iglob,2) + &
 !        2.0d0 * accel_elastic_rk(1,iglob,3) + accel_elastic_rk(1,iglob,4))
!
 !       veloc_elastic(3,iglob) = veloc_elastic_initial_rk(3,iglob) + 1.0d0 / 6.0d0 * &
 !       (accel_elastic_rk(3,iglob,1) + 2.0d0 * accel_elastic_rk(3,iglob,2) + &
 !        2.0d0 * accel_elastic_rk(3,iglob,3) + accel_elastic_rk(3,iglob,4))

 !       displ_elastic(1,iglob) = displ_elastic_initial_rk(1,iglob) + 1.0d0 / 6.0d0 * &
 !       (veloc_elastic_rk(1,iglob,1) + 2.0d0 * veloc_elastic_rk(1,iglob,2) + &
 !        2.0d0 * veloc_elastic_rk(1,iglob,3) + veloc_elastic_rk(1,iglob,4))

 !       displ_elastic(3,iglob) = displ_elastic_initial_rk(3,iglob) + 1.0d0 / 6.0d0 * &
 !       (veloc_elastic_rk(3,iglob,1) + 2.0d0 * veloc_elastic_rk(3,iglob,2) + &
 !        2.0d0 * veloc_elastic_rk(3,iglob,3) + veloc_elastic_rk(3,iglob,4))

 !       endif
        ! updating velocities(poro side)

 !       accels_poroelastic_rk(1,iglob,i_stage) = deltat * accels_poroelastic(1,iglob)
 !       accels_poroelastic_rk(2,iglob,i_stage) = deltat * accels_poroelastic(2,iglob)
 !       velocs_poroelastic_rk(1,iglob,i_stage) = deltat * velocs_poroelastic(1,iglob)
 !       velocs_poroelastic_rk(2,iglob,i_stage) = deltat * velocs_poroelastic(2,iglob)


 !       if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then

 !       if(i_stage == 1)weight_rk = 0.5d0
 !       if(i_stage == 2)weight_rk = 0.5d0
 !       if(i_stage == 3)weight_rk = 1.0d0

 !       velocs_poroelastic(1,iglob) = velocs_poroelastic_initial_rk(1,iglob) + weight_rk * accels_poroelastic_rk(1,iglob,i_stage)
 ! velocs_poroelastic(2,iglob) = velocs_poroelastic_initial_rk(2,iglob) + weight_rk * accels_poroelastic_rk(2,iglob,i_stage)
 !       displs_poroelastic(1,iglob) = displs_poroelastic_initial_rk(1,iglob) + weight_rk * velocs_poroelastic_rk(1,iglob,i_stage)
 ! displs_poroelastic(2,iglob) = displs_poroelastic_initial_rk(2,iglob) + weight_rk * velocs_poroelastic_rk(2,iglob,i_stage)


 !       else if(i_stage==4)then

 !       velocs_poroelastic(1,iglob) = velocs_poroelastic_initial_rk(1,iglob) + 1.0d0 / 6.0d0 * &
 !       (accels_poroelastic_rk(1,iglob,1) + 2.0d0 * accels_poroelastic_rk(1,iglob,2) + &
 !        2.0d0 * accels_poroelastic_rk(1,iglob,3) + accels_poroelastic_rk(1,iglob,4))

 !       velocs_poroelastic(2,iglob) = velocs_poroelastic_initial_rk(2,iglob) + 1.0d0 / 6.0d0 * &
 !       (accels_poroelastic_rk(2,iglob,1) + 2.0d0 * accels_poroelastic_rk(2,iglob,2) + &
 !        2.0d0 * accels_poroelastic_rk(2,iglob,3) + accels_poroelastic_rk(2,iglob,4))
!
 !       displs_poroelastic(1,iglob) = displs_poroelastic_initial_rk(1,iglob) + 1.0d0 / 6.0d0 * &
 !       (velocs_poroelastic_rk(1,iglob,1) + 2.0d0 * velocs_poroelastic_rk(1,iglob,2) + &
 !        2.0d0 * velocs_poroelastic_rk(1,iglob,3) + velocs_poroelastic_rk(1,iglob,4))
!
 !       displs_poroelastic(2,iglob) = displs_poroelastic_initial_rk(2,iglob) + 1.0d0 / 6.0d0 * &
 !       (velocs_poroelastic_rk(2,iglob,1) + 2.0d0 * velocs_poroelastic_rk(2,iglob,2) + &
 !        2.0d0 * velocs_poroelastic_rk(2,iglob,3) + velocs_poroelastic_rk(2,iglob,4))

 !       endif

 !     endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if(SIMULATION_TYPE == 3) then
              b_veloc_elastic(1,iglob) = b_veloc_elastic(1,iglob) - b_deltatover2*b_accel_elastic(1,iglob)
              b_veloc_elastic(3,iglob) = b_veloc_elastic(3,iglob) - b_deltatover2*b_accel_elastic(3,iglob)
              b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) / rmass_inverse_elastic_one(iglob)
              b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) / rmass_inverse_elastic_three(iglob)
              ! recovering original velocities and accelerations on boundaries (poro side)
              b_velocs_poroelastic(1,iglob) = b_velocs_poroelastic(1,iglob) - b_deltatover2*b_accels_poroelastic(1,iglob)
              b_velocs_poroelastic(2,iglob) = b_velocs_poroelastic(2,iglob) - b_deltatover2*b_accels_poroelastic(2,iglob)
              b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) / rmass_s_inverse_poroelastic(iglob)
              b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) / rmass_s_inverse_poroelastic(iglob)
              ! assembling accelerations
              b_accel_elastic(1,iglob) = ( b_accel_elastic(1,iglob) + b_accels_poroelastic(1,iglob) ) / &
                                   ( 1.0/rmass_inverse_elastic_one(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
              b_accel_elastic(3,iglob) = ( b_accel_elastic(3,iglob) + b_accels_poroelastic(2,iglob) ) / &
                                   ( 1.0/rmass_inverse_elastic_three(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
              b_accels_poroelastic(1,iglob) = b_accel_elastic(1,iglob)
              b_accels_poroelastic(2,iglob) = b_accel_elastic(3,iglob)
              ! updating velocities
              b_velocs_poroelastic(1,iglob) = b_velocs_poroelastic(1,iglob) + b_deltatover2*b_accels_poroelastic(1,iglob)
              b_velocs_poroelastic(2,iglob) = b_velocs_poroelastic(2,iglob) + b_deltatover2*b_accels_poroelastic(2,iglob)
              b_veloc_elastic(1,iglob) = b_veloc_elastic(1,iglob) + b_deltatover2*b_accel_elastic(1,iglob)
              b_veloc_elastic(3,iglob) = b_veloc_elastic(3,iglob) + b_deltatover2*b_accel_elastic(3,iglob)
              ! zeros w
              b_accelw_poroelastic(1,iglob) = ZERO
              b_accelw_poroelastic(2,iglob) = ZERO
              b_velocw_poroelastic(1,iglob) = ZERO
              b_velocw_poroelastic(2,iglob) = ZERO
            endif !if(SIMULATION_TYPE == 3)

          endif !if(icount(iglob) ==1)

        enddo

      enddo
    endif

  else !GPU_MODE

    if(any_poroelastic)   call exit_mpi('poroelastic not implemented in GPU MODE yet')

  endif

   enddo !LDDRK or RK

! ********************************************************************************************
!                       reading lastframe for adjoint/kernels calculation
! ********************************************************************************************
    if(it == 1 .and. SIMULATION_TYPE == 3) then

      ! acoustic medium
      if(any_acoustic) then
        write(outputname,'(a,i6.6,a)') 'lastframe_acoustic',myrank,'.bin'
        open(unit=55,file='OUTPUT_FILES/'//outputname,status='old',action='read',form='unformatted')
        do j=1,nglob
          read(55) b_potential_acoustic(j),&
                  b_potential_dot_acoustic(j),&
                  b_potential_dot_dot_acoustic(j)
          enddo
        close(55)



        if(GPU_MODE) then
        ! transfers fields onto GPU
        call transfer_b_fields_ac_to_device(NGLOB_AB,b_potential_acoustic, &
                                            b_potential_dot_acoustic,      &
                                            b_potential_dot_dot_acoustic,  &
                                            Mesh_pointer)

        else
          ! free surface for an acoustic medium
          if ( nelem_acoustic_surface > 0 ) then
            call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                            b_potential_acoustic)
          endif

        endif


       endif




      ! elastic medium
      if(any_elastic) then
        write(outputname,'(a,i6.6,a)') 'lastframe_elastic',myrank,'.bin'
        open(unit=55,file='OUTPUT_FILES/'//outputname,status='old',action='read',form='unformatted')
        if(p_sv)then !P-SV waves
          do j=1,nglob
            read(55) (b_displ_elastic(i,j), i=1,NDIM), &
                      (b_veloc_elastic(i,j), i=1,NDIM), &
                      (b_accel_elastic(i,j), i=1,NDIM)
          enddo

          if(GPU_MODE) then
            b_displ_2D(1,:) = b_displ_elastic(1,:)
            b_displ_2D(2,:) = b_displ_elastic(2,:)
            b_veloc_2D(1,:) = b_veloc_elastic(1,:)
            b_veloc_2D(2,:) = b_veloc_elastic(2,:)
            b_accel_2D(1,:) = b_accel_elastic(1,:)
            b_accel_2D(2,:) = b_accel_elastic(2,:)
            call transfer_b_fields_to_device(NDIM*NGLOB_AB,b_displ_2D,b_veloc_2D,b_accel_2D,Mesh_pointer)
          else
            b_displ_elastic(3,:) = b_displ_elastic(2,:)
            b_displ_elastic(2,:) = 0._CUSTOM_REAL
            b_veloc_elastic(3,:) = b_veloc_elastic(2,:)
            b_veloc_elastic(2,:) = 0._CUSTOM_REAL
            b_accel_elastic(3,:) = b_accel_elastic(2,:)
            b_accel_elastic(2,:) = 0._CUSTOM_REAL
          endif

        else !SH (membrane) waves
          do j=1,nglob
            read(55) b_displ_elastic(2,j), &
                      b_veloc_elastic(2,j), &
                      b_accel_elastic(2,j)
          enddo
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

  if ( NOISE_TOMOGRAPHY == 1 ) then
    call save_surface_movie_noise()

  else if ( NOISE_TOMOGRAPHY == 2 .and. save_everywhere ) then
    call save_surface_movie_noise()

  else if ( NOISE_TOMOGRAPHY == 3 .and. save_everywhere ) then
    if (it==1) &
      open(unit=500,file='OUTPUT_FILES/NOISE_TOMOGRAPHY/phi',access='direct', &
      recl=nglob*CUSTOM_REAL,action='write',iostat=ios)
    if( ios /= 0) write(*,*) 'Error retrieving ensemble forward wavefield.'
    if(p_sv) then
      call exit_mpi('P-SV case not yet implemented.')
    else
      read(unit=500,rec=NSTEP-it+1) b_displ_elastic(2,:)
    endif

  endif



!>NOISE_TOMOGRAPHY


if (GPU_MODE) then


! Kernel calculation
if (SIMULATION_TYPE == 3) then

  if (any_acoustic) call compute_kernels_acoustic_cuda(Mesh_pointer,deltatf)

  if (any_elastic) call compute_kernels_elastic_cuda(Mesh_pointer,deltatf)

  if ( APPROXIMATE_HESS_KL ) then
     ! computes contribution to density and bulk modulus kernel
    call compute_kernels_hess_cuda(Mesh_pointer,any_elastic,any_acoustic)
  endif


! Kernel transfer

   if(it == NSTEP) then

     if( any_acoustic ) then
       call transfer_kernels_ac_to_host(Mesh_pointer,rho_ac_kl,kappa_ac_kl,NSPEC_AB)
     endif

     if( any_elastic ) then
       call transfer_kernels_el_to_host(Mesh_pointer,rho_kl,mu_kl,kappa_kl,NSPEC_AB)


     ! Multiply each kernel point with the local coefficient
        do ispec = 1, nspec
          if(elastic(ispec)) then
            do j = 1, NGLLZ
              do i = 1, NGLLX
                iglob = ibool(i,j,ispec)
                if (.not. assign_external_model) then
                   mul_global(iglob) = poroelastcoef(2,1,kmato(ispec))
                   kappal_global(iglob) = poroelastcoef(3,1,kmato(ispec)) &
                                       - 4._CUSTOM_REAL*mul_global(iglob)/3._CUSTOM_REAL
                   rhol_global(iglob) = density(1,kmato(ispec))
                else
                   rhol_global(iglob)   = rhoext(i,j,ispec)
                   mul_global(iglob)    = rhoext(i,j,ispec)*vsext(i,j,ispec)*vsext(i,j,ispec)
                   kappal_global(iglob) = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec) &
                                       -4._CUSTOM_REAL*mul_global(iglob)/3._CUSTOM_REAL
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


! Simulating seismograms

   if(mod(it-1,subsamp_seismos) == 0 .and. SIMULATION_TYPE == 1) then


    seismo_current = seismo_current + 1

    if ( nrecloc > 0 ) call compute_seismograms_cuda(Mesh_pointer,seismotype,sisux,&
                       sisuz,seismo_current,NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,any_elastic_glob,any_acoustic_glob)

   endif


! Fields transfer for imaging

    if( (output_color_image .and. ( (mod(it,NSTEP_BETWEEN_OUTPUT_IMAGES) == 0 .or. it == 5)) .or. it == NSTEP) ) then


    if( any_acoustic ) &
      call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic, &
                                          potential_dot_acoustic, potential_dot_dot_acoustic, &
                                          Mesh_pointer)


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


if (.NOT. GPU_MODE) then

! ********************************************************************************************
!                                      kernels calculation
! ********************************************************************************************
    if(any_elastic .and. SIMULATION_TYPE == 3) then ! kernels calculation
      do iglob = 1,nglob
        rho_k(iglob) =  accel_elastic(1,iglob)*b_displ_elastic(1,iglob) +&
                            accel_elastic(2,iglob)*b_displ_elastic(2,iglob) +&
                            accel_elastic(3,iglob)*b_displ_elastic(3,iglob)
        rhorho_el_hessian_temp1(iglob) = accel_elastic(1,iglob)*accel_elastic(1,iglob) +&
                                            accel_elastic(2,iglob)*accel_elastic(2,iglob)  +&
                                            accel_elastic(3,iglob)*accel_elastic(3,iglob)
        rhorho_el_hessian_temp2(iglob) = accel_elastic(1,iglob)*b_accel_elastic(1,iglob) +&
                                            accel_elastic(2,iglob)*b_accel_elastic(2,iglob)  +&
                                            accel_elastic(3,iglob)*b_accel_elastic(3,iglob)
      enddo
    endif

    if(any_poroelastic .and. SIMULATION_TYPE == 3) then
      do iglob =1,nglob
        rhot_k(iglob) = accels_poroelastic(1,iglob) * b_displs_poroelastic(1,iglob) + &
                  accels_poroelastic(2,iglob) * b_displs_poroelastic(2,iglob)
        rhof_k(iglob) = accelw_poroelastic(1,iglob) * b_displs_poroelastic(1,iglob) + &
                  accelw_poroelastic(2,iglob) * b_displs_poroelastic(2,iglob) + &
                  accels_poroelastic(1,iglob) * b_displw_poroelastic(1,iglob) + &
                  accels_poroelastic(2,iglob) * b_displw_poroelastic(2,iglob)
        sm_k(iglob) =  accelw_poroelastic(1,iglob) * b_displw_poroelastic(1,iglob) + &
                  accelw_poroelastic(2,iglob) * b_displw_poroelastic(2,iglob)
        eta_k(iglob) = velocw_poroelastic(1,iglob) * b_displw_poroelastic(1,iglob) + &
                  velocw_poroelastic(2,iglob) * b_displw_poroelastic(2,iglob)
      enddo
    endif

!----  compute kinetic and potential energy
    if(output_energy) then

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
    if(mod(it,NSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then
      call check_stability()

    endif

!---- loop on all the receivers to compute and store the seismograms
    if(mod(it-1,subsamp_seismos) == 0) then
      call write_seismograms()
    endif

!----- writing the kernels
    ! kernels output
    if(SIMULATION_TYPE == 3) then

      if(any_acoustic) then

        do ispec = 1, nspec
          if(acoustic(ispec)) then
            do j = 1, NGLLZ
              do i = 1, NGLLX
                iglob = ibool(i,j,ispec)
                if (.not. assign_external_model) then
                   kappal_ac_global(iglob) = poroelastcoef(3,1,kmato(ispec))
                   rhol_ac_global(iglob) = density(1,kmato(ispec))
                else
                   rhol_ac_global(iglob)   = rhoext(i,j,ispec)
                   kappal_ac_global(iglob) = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec)
                endif

! calcul the displacement by computing the gradient of potential / rho
! and calcul the acceleration by computing the gradient of potential_dot_dot / rho
                tempx1l = ZERO
                tempx2l = ZERO
                b_tempx1l = ZERO
                b_tempx2l = ZERO
                bb_tempx1l = ZERO
                bb_tempx2l = ZERO
                do k = 1,NGLLX
                  ! derivative along x
                  !tempx1l = tempx1l + potential_dot_dot_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
                  tempx1l = tempx1l + potential_acoustic(ibool(k,j,ispec))*hprime_xx(i,k) !!! YANGL
                  b_tempx1l = b_tempx1l + b_potential_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
                  bb_tempx1l = bb_tempx1l + b_potential_dot_dot_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
                  ! derivative along z
                  !tempx2l = tempx2l + potential_dot_dot_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
                  tempx2l = tempx2l + potential_acoustic(ibool(i,k,ispec))*hprime_zz(j,k) !!! YANGL
                  b_tempx2l = b_tempx2l + b_potential_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
                  bb_tempx2l = bb_tempx2l + b_potential_dot_dot_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
                enddo

                xixl = xix(i,j,ispec)
                xizl = xiz(i,j,ispec)
                gammaxl = gammax(i,j,ispec)
                gammazl = gammaz(i,j,ispec)

                if(assign_external_model) rhol_ac_global(iglob) = rhoext(i,j,ispec)

                ! derivatives of potential
                accel_ac(1,iglob) = (tempx1l*xixl + tempx2l*gammaxl) / rhol_ac_global(iglob)
                accel_ac(2,iglob) = (tempx1l*xizl + tempx2l*gammazl) / rhol_ac_global(iglob)
                b_displ_ac(1,iglob) = (b_tempx1l*xixl + b_tempx2l*gammaxl) / rhol_ac_global(iglob)
                b_displ_ac(2,iglob) = (b_tempx1l*xizl + b_tempx2l*gammazl) / rhol_ac_global(iglob)
                b_accel_ac(1,iglob) = (bb_tempx1l*xixl + bb_tempx2l*gammaxl) / rhol_ac_global(iglob)
                b_accel_ac(2,iglob) = (bb_tempx1l*xizl + bb_tempx2l*gammazl) / rhol_ac_global(iglob)

              enddo !i = 1, NGLLX
            enddo !j = 1, NGLLZ
          endif
        enddo

        do ispec = 1,nspec
          if(acoustic(ispec)) then
            do j = 1, NGLLZ
              do i = 1, NGLLX
                iglob = ibool(i,j,ispec)
                !<YANGL
                !!!! old expression (from elastic kernels)
                !!!rho_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) - rhol_ac_global(iglob)  * &
                !!!           dot_product(accel_ac(:,iglob),b_displ_ac(:,iglob)) * deltat
                !!!kappa_ac_kl(i,j,ispec) = kappa_ac_kl(i,j,ispec) - kappal_ac_global(iglob) * &
                !!!           potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob) * &
                !!!           b_potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob)&
                !!!           * deltat
                !!!! new expression (from PDE-constrained optimization, coupling terms changed as well)
                rho_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) + rhol_ac_global(iglob)  * &
                           dot_product(accel_ac(:,iglob),b_displ_ac(:,iglob)) * deltat
                kappa_ac_kl(i,j,ispec) = kappa_ac_kl(i,j,ispec) + kappal_ac_global(iglob) * &
                           potential_acoustic(iglob)/kappal_ac_global(iglob) * &
                           b_potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob) * deltat
                !>YANGL
                !
                rhop_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) + kappa_ac_kl(i,j,ispec)
                alpha_ac_kl(i,j,ispec) = TWO *  kappa_ac_kl(i,j,ispec)
                rhorho_ac_hessian_final1(i,j,ispec) =  rhorho_ac_hessian_final1(i,j,ispec) + &
                             dot_product(accel_ac(:,iglob),accel_ac(:,iglob)) * deltat
                rhorho_ac_hessian_final2(i,j,ispec) =  rhorho_ac_hessian_final2(i,j,ispec) + &
                             dot_product(accel_ac(:,iglob),b_accel_ac(:,iglob)) * deltat
              enddo
            enddo
          endif
        enddo

      endif !if(any_acoustic)

      if(any_elastic) then

        do ispec = 1, nspec
          if(elastic(ispec)) then
            do j = 1, NGLLZ
              do i = 1, NGLLX
                iglob = ibool(i,j,ispec)
                if (.not. assign_external_model) then
                   mul_global(iglob) = poroelastcoef(2,1,kmato(ispec))
                   kappal_global(iglob) = poroelastcoef(3,1,kmato(ispec)) &
                                       - 4._CUSTOM_REAL*mul_global(iglob)/3._CUSTOM_REAL
                   rhol_global(iglob) = density(1,kmato(ispec))
                else
                   rhol_global(iglob)   = rhoext(i,j,ispec)
                   mul_global(iglob)    = rhoext(i,j,ispec)*vsext(i,j,ispec)*vsext(i,j,ispec)
                   kappal_global(iglob) = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec) &
                                       -4._CUSTOM_REAL*mul_global(iglob)/3._CUSTOM_REAL
                endif

                rho_kl(i,j,ispec) = rho_kl(i,j,ispec) - rhol_global(iglob)  * rho_k(iglob) * deltat
                mu_kl(i,j,ispec) =  mu_kl(i,j,ispec) - TWO * mul_global(iglob) * mu_k(iglob) * deltat
                kappa_kl(i,j,ispec) = kappa_kl(i,j,ispec) - kappal_global(iglob) * kappa_k(iglob) * deltat
                !
                rhop_kl(i,j,ispec) = rho_kl(i,j,ispec) + kappa_kl(i,j,ispec) + mu_kl(i,j,ispec)
                beta_kl(i,j,ispec) = TWO * (mu_kl(i,j,ispec) - 4._CUSTOM_REAL * mul_global(iglob) &
                    / (3._CUSTOM_REAL * kappal_global(iglob)) * kappa_kl(i,j,ispec))
                alpha_kl(i,j,ispec) = TWO * (1._CUSTOM_REAL + 4._CUSTOM_REAL * mul_global(iglob)/&
                     (3._CUSTOM_REAL * kappal_global(iglob))) * kappa_kl(i,j,ispec)
                rhorho_el_hessian_final1(i,j,ispec) = rhorho_el_hessian_final1(i,j,ispec) &
                                                  + rhorho_el_hessian_temp1(iglob) * deltat
                rhorho_el_hessian_final2(i,j,ispec) = rhorho_el_hessian_final2(i,j,ispec) &
                                                  + rhorho_el_hessian_temp2(iglob) * deltat

              enddo
            enddo
          endif
        enddo

      endif !if(any_elastic)

      if(any_poroelastic) then

        do ispec = 1, nspec
          if(poroelastic(ispec)) then
            do j = 1, NGLLZ
              do i = 1, NGLLX
                iglob = ibool(i,j,ispec)
                phil_global(iglob) = porosity(kmato(ispec))
                tortl_global(iglob) = tortuosity(kmato(ispec))
                rhol_s_global(iglob) = density(1,kmato(ispec))
                rhol_f_global(iglob) = density(2,kmato(ispec))
                rhol_bar_global(iglob) =  (1._CUSTOM_REAL - phil_global(iglob))*rhol_s_global(iglob) &
                  + phil_global(iglob)*rhol_f_global(iglob)
                etal_f_global(iglob) = poroelastcoef(2,2,kmato(ispec))
                permlxx_global(iglob) = permeability(1,kmato(ispec))
                permlxz_global(iglob) = permeability(2,kmato(ispec))
                permlzz_global(iglob) = permeability(3,kmato(ispec))
                mulfr_global(iglob) = poroelastcoef(2,3,kmato(ispec))

                rhot_kl(i,j,ispec) = rhot_kl(i,j,ispec) - deltat * rhol_bar_global(iglob) * rhot_k(iglob)
                rhof_kl(i,j,ispec) = rhof_kl(i,j,ispec) - deltat * rhol_f_global(iglob) * rhof_k(iglob)
                sm_kl(i,j,ispec) = sm_kl(i,j,ispec) - &
                        deltat * rhol_f_global(iglob)*tortl_global(iglob)/phil_global(iglob) * sm_k(iglob)
                !at the moment works with constant permeability
                eta_kl(i,j,ispec) = eta_kl(i,j,ispec) - deltat * etal_f_global(iglob)/permlxx_global(iglob) * eta_k(iglob)
                B_kl(i,j,ispec) = B_kl(i,j,ispec) - deltat * B_k(iglob)
                C_kl(i,j,ispec) = C_kl(i,j,ispec) - deltat * C_k(iglob)
                M_kl(i,j,ispec) = M_kl(i,j,ispec) - deltat * M_k(iglob)
                mufr_kl(i,j,ispec) = mufr_kl(i,j,ispec) - TWO * deltat * mufr_k(iglob)
                ! density kernels
                rholb = rhol_bar_global(iglob) - phil_global(iglob)*rhol_f_global(iglob)/tortl_global(iglob)
                rhob_kl(i,j,ispec) = rhot_kl(i,j,ispec) + B_kl(i,j,ispec) + mufr_kl(i,j,ispec)
                rhofb_kl(i,j,ispec) = rhof_kl(i,j,ispec) + C_kl(i,j,ispec) + M_kl(i,j,ispec) + sm_kl(i,j,ispec)
                Bb_kl(i,j,ispec) = B_kl(i,j,ispec)
                Cb_kl(i,j,ispec) = C_kl(i,j,ispec)
                Mb_kl(i,j,ispec) = M_kl(i,j,ispec)
                mufrb_kl(i,j,ispec) = mufr_kl(i,j,ispec)
                phi_kl(i,j,ispec) = - sm_kl(i,j,ispec) - M_kl(i,j,ispec)
                ! wave speed kernels
                dd1 = (1._CUSTOM_REAL+rholb/rhol_f_global(iglob))*ratio**2 &
                      + 2._CUSTOM_REAL*ratio &
                      + tortl_global(iglob)/phil_global(iglob)

                rhobb_kl(i,j,ispec) = rhob_kl(i,j,ispec) - &
                      phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*B_biot) * &
                      (cpIIsquare + (cpIsquare - cpIIsquare)*( (phil_global(iglob) / &
                      tortl_global(iglob)*ratio +1._CUSTOM_REAL)/dd1 + &
                      (rhol_bar_global(iglob)**2*ratio**2/rhol_f_global(iglob)**2*(phil_global(iglob) / &
                      tortl_global(iglob)*ratio+1)*(phil_global(iglob)/tortl_global(iglob)*ratio + &
                      phil_global(iglob)/tortl_global(iglob) * &
                      (1+rhol_f_global(iglob)/rhol_bar_global(iglob))-1) )/dd1**2 ) - &
                      FOUR_THIRDS*cssquare )*Bb_kl(i,j,ispec) - &
                      rhol_bar_global(iglob)*ratio**2/M_biot * (cpIsquare - cpIIsquare)* &
                      (phil_global(iglob)/tortl_global(iglob)*ratio + &
                      1._CUSTOM_REAL)**2/dd1**2*Mb_kl(i,j,ispec) + &
                      rhol_bar_global(iglob)*ratio/C_biot * (cpIsquare - cpIIsquare)* (&
                      (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                      phil_global(iglob)*ratio/tortl_global(iglob)*(phil_global(iglob) / &
                      tortl_global(iglob)*ratio+1._CUSTOM_REAL)*&
                      (1+rhol_bar_global(iglob)*ratio/rhol_f_global(iglob))/dd1**2)*Cb_kl(i,j,ispec)+ &
                      phil_global(iglob)*rhol_f_global(iglob)*cssquare / &
                      (tortl_global(iglob)*mulfr_global(iglob))*mufrb_kl(i,j,ispec)
                rhofbb_kl(i,j,ispec) = rhofb_kl(i,j,ispec) + &
                        phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*B_biot) * &
                       (cpIIsquare + (cpIsquare - cpIIsquare)*( (phil_global(iglob)/ &
                       tortl_global(iglob)*ratio +1._CUSTOM_REAL)/dd1+&
                       (rhol_bar_global(iglob)**2*ratio**2/rhol_f_global(iglob)**2*(phil_global(iglob)/ &
                       tortl_global(iglob)*ratio+1)*(phil_global(iglob)/tortl_global(iglob)*ratio+ &
                       phil_global(iglob)/tortl_global(iglob)*&
                       (1+rhol_f_global(iglob)/rhol_bar_global(iglob))-1) )/dd1**2 )- &
                       FOUR_THIRDS*cssquare )*Bb_kl(i,j,ispec) + &
                        rhol_bar_global(iglob)*ratio**2/M_biot * (cpIsquare - cpIIsquare)* &
                       (phil_global(iglob)/tortl_global(iglob)*ratio + &
                       1._CUSTOM_REAL)**2/dd1**2*Mb_kl(i,j,ispec) - &
                       rhol_bar_global(iglob)*ratio/C_biot * (cpIsquare - cpIIsquare)* (&
                       (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                       phil_global(iglob)*ratio/tortl_global(iglob)*(phil_global(iglob)/ &
                       tortl_global(iglob)*ratio+1._CUSTOM_REAL)*&
                       (1+rhol_bar_global(iglob)*ratio/rhol_f_global(iglob))/dd1**2)*Cb_kl(i,j,ispec)- &
                       phil_global(iglob)*rhol_f_global(iglob)*cssquare/ &
                       (tortl_global(iglob)*mulfr_global(iglob))*mufrb_kl(i,j,ispec)
                phib_kl(i,j,ispec) = phi_kl(i,j,ispec) - &
                       phil_global(iglob)*rhol_bar_global(iglob)/(tortl_global(iglob)*B_biot) &
                       * ( cpIsquare - rhol_f_global(iglob)/rhol_bar_global(iglob)*cpIIsquare- &
                       (cpIsquare-cpIIsquare)*( (TWO*ratio**2*phil_global(iglob)/ &
                       tortl_global(iglob) + (1._CUSTOM_REAL+&
                       rhol_f_global(iglob)/rhol_bar_global(iglob))* &
                       (TWO*ratio*phil_global(iglob)/tortl_global(iglob)+&
                       1._CUSTOM_REAL))/dd1 + (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                       1._CUSTOM_REAL)*(phil_global(iglob)*&
                       ratio/tortl_global(iglob)+phil_global(iglob)/tortl_global(iglob)* &
                       (1._CUSTOM_REAL+rhol_f_global(iglob)/&
                       rhol_bar_global(iglob))-1._CUSTOM_REAL)*((1._CUSTOM_REAL+ &
                       rhol_bar_global(iglob)/rhol_f_global(iglob)-&
                       TWO*phil_global(iglob)/tortl_global(iglob))*ratio**2+TWO*ratio)/dd1**2 ) - &
                       FOUR_THIRDS*rhol_f_global(iglob)*cssquare/rhol_bar_global(iglob) )*Bb_kl(i,j,ispec) + &
                       rhol_f_global(iglob)/M_biot * (cpIsquare-cpIIsquare)*(&
                       TWO*ratio*(phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                       (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2*( &
                       (1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                       rhol_f_global(iglob)-TWO*phil_global(iglob)/tortl_global(iglob))*ratio**2+TWO*ratio)/dd1**2 &
                       )*Mb_kl(i,j,ispec) + &
                       phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*C_biot)* &
                       (cpIsquare-cpIIsquare)*ratio* (&
                       (1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob)*ratio)/dd1 - &
                       (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)* &
                       (1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                       rhol_f_global(iglob)*ratio)*((1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)-TWO*&
                       phil_global(iglob)/tortl_global(iglob))*ratio+TWO)/dd1**2&
                        )*Cb_kl(i,j,ispec) -&
                       phil_global(iglob)*rhol_f_global(iglob)*cssquare &
                       /(tortl_global(iglob)*mulfr_global(iglob))*mufrb_kl(i,j,ispec)
                cpI_kl(i,j,ispec) = 2._CUSTOM_REAL*cpIsquare/B_biot*rhol_bar_global(iglob)*( &
                       1._CUSTOM_REAL-phil_global(iglob)/tortl_global(iglob) + &
                       (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                       1._CUSTOM_REAL)*(phil_global(iglob)/tortl_global(iglob)*&
                       ratio+phil_global(iglob)/tortl_global(iglob)* &
                       (1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-&
                       1._CUSTOM_REAL)/dd1 &
                        )* Bb_kl(i,j,ispec) +&
                       2._CUSTOM_REAL*cpIsquare*rhol_f_global(iglob)*tortl_global(iglob)/(phil_global(iglob)*M_biot) *&
                       (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2/dd1*Mb_kl(i,j,ispec)+&
                       2._CUSTOM_REAL*cpIsquare*rhol_f_global(iglob)/C_biot * &
                       (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)* &
                       (1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                       rhol_f_global(iglob)*ratio)/dd1*Cb_kl(i,j,ispec)
                cpII_kl(i,j,ispec) = 2._CUSTOM_REAL*cpIIsquare*rhol_bar_global(iglob)/B_biot * (&
                       phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*rhol_bar_global(iglob)) - &
                       (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                       1._CUSTOM_REAL)*(phil_global(iglob)/tortl_global(iglob)*&
                       ratio+phil_global(iglob)/tortl_global(iglob)* &
                       (1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-&
                       1._CUSTOM_REAL)/dd1  ) * Bb_kl(i,j,ispec) +&
                       2._CUSTOM_REAL*cpIIsquare*rhol_f_global(iglob)*tortl_global(iglob)/(phil_global(iglob)*M_biot) * (&
                       1._CUSTOM_REAL - (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                       1._CUSTOM_REAL)**2/dd1  )*Mb_kl(i,j,ispec) + &
                       2._CUSTOM_REAL*cpIIsquare*rhol_f_global(iglob)/C_biot * (&
                       1._CUSTOM_REAL - (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                       1._CUSTOM_REAL)*(1._CUSTOM_REAL+&
                       rhol_bar_global(iglob)/rhol_f_global(iglob)*ratio)/dd1  )*Cb_kl(i,j,ispec)
                cs_kl(i,j,ispec) = - 8._CUSTOM_REAL/3._CUSTOM_REAL*cssquare* &
                       rhol_bar_global(iglob)/B_biot*(1._CUSTOM_REAL-&
                       phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)* &
                       rhol_bar_global(iglob)))*Bb_kl(i,j,ispec) + &
                       2._CUSTOM_REAL*(rhol_bar_global(iglob)-rhol_f_global(iglob)*&
                       phil_global(iglob)/tortl_global(iglob))/&
                       mulfr_global(iglob)*cssquare*mufrb_kl(i,j,ispec)
                ratio_kl(i,j,ispec) = ratio*rhol_bar_global(iglob)*phil_global(iglob)/(tortl_global(iglob)*B_biot) * &
                       (cpIsquare-cpIIsquare) * ( &
                       phil_global(iglob)/tortl_global(iglob)*(2._CUSTOM_REAL*ratio+1._CUSTOM_REAL+rhol_f_global(iglob)/ &
                       rhol_bar_global(iglob))/dd1 - (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*&
                       (phil_global(iglob)/tortl_global(iglob)*ratio+phil_global(iglob)/tortl_global(iglob)*(&
                       1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-1._CUSTOM_REAL)*(2._CUSTOM_REAL*ratio*(&
                       1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)-phil_global(iglob)/tortl_global(iglob)) +&
                       2._CUSTOM_REAL)/dd1**2  )*Bb_kl(i,j,ispec) + &
                       ratio*rhol_f_global(iglob)*tortl_global(iglob)/(phil_global(iglob)*M_biot)*(cpIsquare-cpIIsquare) * &
                       2._CUSTOM_REAL*phil_global(iglob)/tortl_global(iglob) * (&
                       (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                       (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2*( &
                       (1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                       rhol_f_global(iglob)-phil_global(iglob)/tortl_global(iglob))*ratio+ &
                       1._CUSTOM_REAL)/dd1**2 )*Mb_kl(i,j,ispec) +&
                       ratio*rhol_f_global(iglob)/C_biot*(cpIsquare-cpIIsquare) * (&
                       (2._CUSTOM_REAL*phil_global(iglob)*rhol_bar_global(iglob)* &
                       ratio/(tortl_global(iglob)*rhol_f_global(iglob))+&
                       phil_global(iglob)/tortl_global(iglob)+rhol_bar_global(iglob)/rhol_f_global(iglob))/dd1 - &
                       2._CUSTOM_REAL*phil_global(iglob)/tortl_global(iglob)*(phil_global(iglob)/tortl_global(iglob)*ratio+&
                       1._CUSTOM_REAL)*(1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)*ratio)*((1._CUSTOM_REAL+&
                       rhol_bar_global(iglob)/rhol_f_global(iglob)- &
                       phil_global(iglob)/tortl_global(iglob))*ratio+1._CUSTOM_REAL)/&
                       dd1**2 )*Cb_kl(i,j,ispec)

              enddo
            enddo
          endif
        enddo

      endif ! if(any_poroelastic)

    endif ! if(SIMULATION_TYPE == 3)


  endif !Not GPU_MODE

!
!----  display results at given time steps
!
    if(mod(it,NSTEP_BETWEEN_OUTPUT_IMAGES) == 0 .or. it == 5 .or. it == NSTEP) then

!
! write kernel files
!

      if(SIMULATION_TYPE == 3 .and. it == NSTEP) then
          call save_adjoint_kernels()
      endif


!<NOISE_TOMOGRAPHY

if (.NOT. GPU_MODE ) then

      if (NOISE_TOMOGRAPHY == 3 .and. output_wavefields_noise) then

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




!
!----  PostScript display
!
      if(output_postscript_snapshot) then

        if (myrank == 0) then
          write(IOUT,*)
          write(IOUT,*) 'Writing PostScript vector plot for time step ',it
        endif

        if(imagetype_postscript == 1 .and. p_sv) then

          if (myrank == 0) write(IOUT,*) 'drawing displacement vector as small arrows...'

          call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                          potential_gravito,displ_elastic,displs_poroelastic)

          call plotpost()

        else if(imagetype_postscript == 2 .and. p_sv) then

          if (myrank == 0) write(IOUT,*) 'drawing velocity vector as small arrows...'

          call compute_vector_whole_medium(potential_dot_acoustic,potential_dot_gravitoacoustic, &
                          potential_dot_gravito,veloc_elastic,velocs_poroelastic)

          call plotpost()

        else if(imagetype_postscript == 3 .and. p_sv) then

          if (myrank == 0) write(IOUT,*) 'drawing acceleration vector as small arrows...'

          call compute_vector_whole_medium(potential_dot_dot_acoustic,potential_dot_dot_gravitoacoustic, &
                          potential_dot_dot_gravito,accel_elastic,accels_poroelastic)

          call plotpost()

        else if(.not. p_sv) then
          call exit_MPI('cannot draw a SH scalar field as a vector plot, turn PostScript plots off')

        else
          call exit_MPI('wrong type for PostScript snapshots')
        endif

        if (myrank == 0 .and. imagetype_postscript /= 4 .and. p_sv) write(IOUT,*) 'PostScript file written'

      endif

!
!----  display color image
!
      if(output_color_image) then

        if (myrank == 0) then
          write(IOUT,*)
          write(IOUT,*) 'Creating color image of size ',NX_IMAGE_color,' x ',NZ_IMAGE_color,' for time step ',it
        endif

        if(imagetype_JPEG >= 1 .and. imagetype_JPEG <= 3) then

          if (myrank == 0) write(IOUT,*) 'drawing scalar image of part of the displacement vector...'
          call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                          potential_gravito,displ_elastic,displs_poroelastic)

        else if(imagetype_JPEG >= 4 .and. imagetype_JPEG <= 6) then

          if (myrank == 0) write(IOUT,*) 'drawing scalar image of part of the velocity vector...'
          call compute_vector_whole_medium(potential_dot_acoustic,potential_dot_gravitoacoustic, &
                          potential_dot_gravito,veloc_elastic,velocs_poroelastic)

        else if(imagetype_JPEG >= 7 .and. imagetype_JPEG <= 9) then

          if (myrank == 0) write(IOUT,*) 'drawing scalar image of part of the acceleration vector...'
          call compute_vector_whole_medium(potential_dot_dot_acoustic,potential_dot_dot_gravitoacoustic, &
                          potential_dot_dot_gravito,accel_elastic,accels_poroelastic)

        else if(imagetype_JPEG >= 11 .and. imagetype_JPEG <= 13) then
! allocation for normalized representation in JPEG image
! for an atmosphere model

          if (myrank == 0) write(IOUT,*) 'drawing scalar image of part of normalized displacement vector...'

          call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                          potential_gravito,displ_elastic,displs_poroelastic)


          do ispec = 1,nspec
            do j = 1,NGLLZ
              do i = 1,NGLLX
                iglob = ibool(i,j,ispec)
                vector_field_display(1,iglob) = sqrt(rhoext(i,j,ispec)) * vector_field_display(1,iglob)
                vector_field_display(2,iglob) = sqrt(rhoext(i,j,ispec)) * vector_field_display(2,iglob)
                vector_field_display(3,iglob) = sqrt(rhoext(i,j,ispec)) * vector_field_display(3,iglob)
              enddo
            enddo
          enddo

        else if(imagetype_JPEG >= 14 .and. imagetype_JPEG <= 16) then
! allocation for normalized representation in JPEG image
! for an atmosphere model
          call compute_vector_whole_medium(potential_dot_acoustic,potential_dot_gravitoacoustic, &
                          potential_dot_gravito,veloc_elastic,velocs_poroelastic)

          do ispec = 1,nspec
            do j = 1,NGLLZ
              do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            vector_field_display(1,iglob) = sqrt(rhoext(i,j,ispec)) * vector_field_display(1,iglob)
            vector_field_display(2,iglob) = sqrt(rhoext(i,j,ispec)) * vector_field_display(2,iglob)
            vector_field_display(3,iglob) = sqrt(rhoext(i,j,ispec)) * vector_field_display(3,iglob)
              enddo
            enddo
          enddo

        else if(imagetype_JPEG == 10 .and. p_sv) then

          if (myrank == 0) write(IOUT,*) 'drawing image of pressure field...'
          call compute_pressure_whole_medium()

        else if(imagetype_JPEG == 10 .and. .not. p_sv) then
          call exit_MPI('cannot draw pressure field for SH (membrane) waves')

        else
          call exit_MPI('wrong type for JPEG snapshots')
        endif

!! DK DK quick hack to remove the PMLs from JPEG images if needed: set the vector field to zero there
        if(PML_BOUNDARY_CONDITIONS .and. REMOVE_PMLS_FROM_JPEG_IMAGES) then
          do ispec = 1,nspec
            if(is_PML(ispec)) then
              do j = 1,NGLLZ
                do i = 1,NGLLX
                  iglob = ibool(i,j,ispec)
                  vector_field_display(1,iglob) = 0.d0
                  vector_field_display(2,iglob) = 0.d0
                  vector_field_display(3,iglob) = 0.d0
                enddo
              enddo
            endif
          enddo
        endif
!! DK DK quick hack to remove the PMLs from JPEG images if needed

        image_color_data(:,:) = 0.d0

        do k = 1, nb_pixel_loc
          j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
          i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color

          ! avoid edge effects
          if(i < 1) i = 1
          if(j < 1) j = 1

          if(i > NX_IMAGE_color) i = NX_IMAGE_color
          if(j > NZ_IMAGE_color) j = NZ_IMAGE_color

          if(p_sv) then ! P-SH waves, plot a component of vector, its norm, or else pressure
            if(iglob_image_color(i,j) /= -1) then
              if(imagetype_JPEG == 1 .or. imagetype_JPEG == 4 .or. imagetype_JPEG == 7 .or. imagetype_JPEG == 11 &
                                     .or. imagetype_JPEG == 14) then
                image_color_data(i,j) = vector_field_display(1,iglob_image_color(i,j))  ! draw the X component of the vector

              else if(imagetype_JPEG == 2 .or. imagetype_JPEG == 5 .or. imagetype_JPEG == 8 .or. imagetype_JPEG == 12 &
                                          .or. imagetype_JPEG == 15) then
                image_color_data(i,j) = vector_field_display(3,iglob_image_color(i,j))  ! draw the Z component of the vector
              else if(imagetype_JPEG == 3 .or. imagetype_JPEG == 6 .or. imagetype_JPEG == 9 .or. imagetype_JPEG == 13 &
                                          .or. imagetype_JPEG == 16) then
                image_color_data(i,j) = sqrt(vector_field_display(1,iglob_image_color(i,j))**2 + &
                                             vector_field_display(3,iglob_image_color(i,j))**2)  ! draw the norm of the vector

              else if(imagetype_JPEG == 10) then
! by convention we have stored pressure in the third component of the array
                image_color_data(i,j) = vector_field_display(3,iglob_image_color(i,j))

              else
                call exit_MPI('wrong type for JPEG snapshots')
              endif
            endif

          else ! SH (membrane) waves, plot y-component
            if(iglob_image_color(i,j) /= -1) image_color_data(i,j) = vector_field_display(2,iglob_image_color(i,j))
          endif
        enddo

! assembling array image_color_data on process zero for color output
#ifdef USE_MPI

        if (nproc > 1) then
          if (myrank == 0) then
            do iproc = 1, nproc-1
              call MPI_RECV(data_pixel_recv(1),nb_pixel_per_proc(iproc+1), MPI_DOUBLE_PRECISION, &
                  iproc, 43, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

              do k = 1, nb_pixel_per_proc(iproc+1)
                j = ceiling(real(num_pixel_recv(k,iproc+1)) / real(NX_IMAGE_color))
                i = num_pixel_recv(k,iproc+1) - (j-1)*NX_IMAGE_color

                ! avoid edge effects
                if(i < 1) i = 1
                if(j < 1) j = 1

                if(i > NX_IMAGE_color) i = NX_IMAGE_color
                if(j > NZ_IMAGE_color) j = NZ_IMAGE_color

                image_color_data(i,j) = data_pixel_recv(k)
              enddo
            enddo
          else
            do k = 1, nb_pixel_loc
              j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
              i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color

              ! avoid edge effects
              if(i < 1) i = 1
              if(j < 1) j = 1

              if(i > NX_IMAGE_color) i = NX_IMAGE_color
              if(j > NZ_IMAGE_color) j = NZ_IMAGE_color

              if(p_sv) then ! P-SH waves, plot a component of vector, its norm, or else pressure

              if(imagetype_JPEG == 1 .or. imagetype_JPEG == 4 .or. imagetype_JPEG == 7 .or. imagetype_JPEG == 11 &
                                     .or. imagetype_JPEG == 14) then
                  data_pixel_send(k) = vector_field_display(1,iglob_image_color(i,j))  ! draw the X component of the vector

              else if(imagetype_JPEG == 2 .or. imagetype_JPEG == 5 .or. imagetype_JPEG == 8 .or. imagetype_JPEG == 12 &
                                          .or. imagetype_JPEG == 15) then
                  data_pixel_send(k) = vector_field_display(3,iglob_image_color(i,j))  ! draw the Z component of the vector

              else if(imagetype_JPEG == 3 .or. imagetype_JPEG == 6 .or. imagetype_JPEG == 9 .or. imagetype_JPEG == 13 &
                                          .or. imagetype_JPEG == 16) then
                  data_pixel_send(k) = sqrt(vector_field_display(1,iglob_image_color(i,j))**2 + &
                                            vector_field_display(3,iglob_image_color(i,j))**2)  ! draw the norm of the vector

                else if(imagetype_JPEG == 10) then
! by convention we have stored pressure in the third component of the array

                  data_pixel_send(k) = vector_field_display(3,iglob_image_color(i,j))

                else
                  call exit_MPI('wrong type for JPEG snapshots')
                endif

              else ! SH (membrane) waves, plot y-component
                if(iglob_image_color(i,j) /= -1) data_pixel_send(k) = vector_field_display(2,iglob_image_color(i,j))
              endif
            enddo

            call MPI_SEND(data_pixel_send(1),nb_pixel_loc,MPI_DOUBLE_PRECISION, 0, 43, MPI_COMM_WORLD, ier)
          endif
        endif
#endif

        if (myrank == 0) then
          call create_color_image()
          write(IOUT,*) 'Color image created'
        endif

      endif

    endif  ! of display images at a given time step

!----------------------------------------------

! dump the full (local) wavefield to a file
! note: in the case of MPI, in the future it would be more convenient to output a single file rather than one for each myrank

    if(output_wavefield_dumps .and. (mod(it,NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS) == 0 .or. it == 5 .or. it == NSTEP)) then

      if (myrank == 0) then
        write(IOUT,*)
        write(IOUT,*) 'Dumping the wave field to a file for time step ',it
      endif

      if(this_is_the_first_time_we_dump) then

        allocate(mask_ibool(nglob))

! save the grid separately once and for all
        if(use_binary_for_wavefield_dumps) then
          write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_for_dumps_',i3.3,'.bin')") myrank
          open(unit=27,file=wavefield_file,form='unformatted',access='direct',status='unknown', &
               action='write',recl=2*SIZE_REAL)
        else
          write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_for_dumps_',i3.3,'.txt')") myrank
          open(unit=27,file=wavefield_file,status='unknown',action='write')
        endif

        icounter = 0
        mask_ibool(:) = .false.
        do ispec = 1,nspec
          do j = 1,NGLLZ
            do i = 1,NGLLX
               iglob = ibool(i,j,ispec)
               if(.not. mask_ibool(iglob)) then
                 icounter = icounter + 1
                 mask_ibool(iglob) = .true.
                 if(use_binary_for_wavefield_dumps) then
                   write(27,rec=icounter) sngl(coord(1,iglob)),sngl(coord(2,iglob))
                 else
                   write(27,'(2e16.6)') coord(1,iglob),coord(2,iglob)
                 endif
               endif
            enddo
          enddo
        enddo

        close(27)

! save nglob to a file once and for all
        write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_value_of_nglob_',i3.3,'.txt')") myrank
        open(unit=27,file=wavefield_file,status='unknown',action='write')
        write(27,*) icounter
        close(27)
        if(icounter /= nglob) stop 'error: should have icounter == nglob in wavefield dumps'

        this_is_the_first_time_we_dump = .false.

      endif

        if(imagetype_wavefield_dumps == 1) then

          if (myrank == 0) write(IOUT,*) 'dumping the displacement vector...'
          call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                            potential_gravito,displ_elastic,displs_poroelastic)
        else if(imagetype_wavefield_dumps == 2) then

          if (myrank == 0) write(IOUT,*) 'dumping the velocity vector...'
          call compute_vector_whole_medium(potential_dot_acoustic,potential_gravitoacoustic, &
                            potential_gravito,veloc_elastic,velocs_poroelastic)

        else if(imagetype_wavefield_dumps == 3) then

          if (myrank == 0) write(IOUT,*) 'dumping the acceleration vector...'
          call compute_vector_whole_medium(potential_dot_dot_acoustic,potential_gravitoacoustic, &
                            potential_gravito,accel_elastic,accels_poroelastic)

        else if(imagetype_wavefield_dumps == 4 .and. p_sv) then

          if (myrank == 0) write(IOUT,*) 'dumping the pressure field...'
          call compute_pressure_whole_medium()

        else if(imagetype_wavefield_dumps == 4 .and. .not. p_sv) then
          call exit_MPI('cannot dump the pressure field for SH (membrane) waves')

        else
          call exit_MPI('wrong type of flag for wavefield dumping')
        endif

        if(use_binary_for_wavefield_dumps) then
          if(p_sv .and. .not. imagetype_wavefield_dumps == 4) then
            nb_of_values_to_save = 2
          else
            nb_of_values_to_save = 1
          endif
          write(wavefield_file,"('OUTPUT_FILES/wavefield',i7.7,'_',i2.2,'_',i3.3,'.bin')") it,SIMULATION_TYPE,myrank
          open(unit=27,file=wavefield_file,form='unformatted',access='direct',status='unknown', &
                       action='write',recl=nb_of_values_to_save*SIZE_REAL)
        else
          write(wavefield_file,"('OUTPUT_FILES/wavefield',i7.7,'_',i2.2,'_',i3.3,'.txt')") it,SIMULATION_TYPE,myrank
          open(unit=27,file=wavefield_file,status='unknown',action='write')
        endif

        icounter = 0
        mask_ibool(:) = .false.
        do ispec = 1,nspec
          do j = 1,NGLLZ
            do i = 1,NGLLX
               iglob = ibool(i,j,ispec)
               if(.not. mask_ibool(iglob)) then
                 icounter = icounter + 1
                 mask_ibool(iglob) = .true.
                 if(use_binary_for_wavefield_dumps) then

                   if(p_sv .and. .not. imagetype_wavefield_dumps == 4) then
                     write(27,rec=icounter) sngl(vector_field_display(1,iglob)),sngl(vector_field_display(3,iglob))
                   else if(p_sv .and. imagetype_wavefield_dumps == 4) then
! by convention we use the third component of the array to store the pressure above
                     write(27,rec=icounter) sngl(vector_field_display(3,iglob))
                   else ! SH case
                     write(27,rec=icounter) sngl(vector_field_display(2,iglob))
                   endif

                 else

                   if(p_sv .and. .not. imagetype_wavefield_dumps == 4) then
                     write(27,*) sngl(vector_field_display(1,iglob)),sngl(vector_field_display(3,iglob))
                   else if(p_sv .and. imagetype_wavefield_dumps == 4) then
! by convention we use the third component of the array to store the pressure above
                     write(27,*) sngl(vector_field_display(3,iglob))
                   else ! SH case
                     write(27,*) sngl(vector_field_display(2,iglob))
                   endif

                 endif
               endif
            enddo
          enddo
        enddo

        close(27)

        if(myrank ==0) write(IOUT,*) 'Wave field dumped'

    endif  ! of display wavefield dumps at a given time step

!----  save temporary or final seismograms
! suppress seismograms if we generate traces of the run for analysis with "ParaVer", because time consuming
    if(mod(it,NSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then

      call write_seismograms_to_file(x_source(1),z_source(1))

      seismo_offset = seismo_offset + seismo_current
      seismo_current = 0

    endif  ! of display images at a given time step

  enddo ! end of the main time loop

! *********************************************************
! ************* END MAIN LOOP OVER THE TIME STEPS *********
! *********************************************************



!
!----  formats
!

 400 format(/1x,41('=')/,' =  T i m e  e v o l u t i o n  l o o p  ='/1x,41('=')/)


end subroutine
