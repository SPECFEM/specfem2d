
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

! write seismograms to text files

  subroutine write_seismograms

  use specfem_par

  implicit none

  integer :: i, j, iglob, irecloc, irec, ispec



! update position in seismograms
    seismo_current = seismo_current + 1

    do irecloc = 1,nrecloc

      irec = recloc(irecloc)

      ispec = ispec_selected_rec(irec)

      ! compute pressure in this element if needed
      if(seismotype == 4) then

        call compute_pressure_one_element(ispec)

      else if(acoustic(ispec)) then

        ! for acoustic medium, compute vector field from gradient of potential for seismograms
        if(seismotype == 1 .or. seismotype == 7) then
          call compute_vector_one_element(potential_acoustic,potential_gravitoacoustic, &
                              potential_gravito,displ_elastic,displs_poroelastic,ispec)
        else if(seismotype == 2) then
          call compute_vector_one_element(potential_dot_acoustic,potential_dot_gravitoacoustic, &
                              potential_dot_gravito,veloc_elastic,velocs_poroelastic,ispec)
        else if(seismotype == 3) then
          call compute_vector_one_element(potential_dot_dot_acoustic,potential_dot_dot_gravitoacoustic, &
                              potential_dot_dot_gravito,accel_elastic,accels_poroelastic,ispec)
        endif

      else if(gravitoacoustic(ispec)) then

        ! for acoustic medium, compute vector field from gradient of potential for seismograms
        if(seismotype == 1 .or. seismotype == 7) then
          call compute_vector_one_element(potential_acoustic,potential_gravitoacoustic, &
                              potential_gravito,displ_elastic,displs_poroelastic,ispec)
        else if(seismotype == 2) then
          call compute_vector_one_element(potential_dot_acoustic,potential_dot_gravitoacoustic, &
                              potential_dot_gravito,veloc_elastic,velocs_poroelastic,ispec)
        else if(seismotype == 3) then
          call compute_vector_one_element(potential_dot_dot_acoustic,potential_dot_dot_gravitoacoustic, &
                              potential_dot_dot_gravito,accel_elastic,accels_poroelastic,ispec)
        endif

      else if(seismotype == 5) then
        call compute_curl_one_element(ispec)
      endif

      ! perform the general interpolation using Lagrange polynomials
      valux = ZERO
      valuy = ZERO
      valuz = ZERO

      valcurl = ZERO

      dxd = 0
      dyd = 0
      dzd = 0

      dcurld = 0
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          hlagrange = hxir_store(irec,i)*hgammar_store(irec,j)
          dcurld=ZERO

          if(seismotype == 4) then
            if (USE_TRICK_FOR_BETTER_PRESSURE .and. (acoustic(ispec) .or. gravitoacoustic(ispec))) then
              ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
              ! use the second derivative of the source for the source time function instead of the source itself,
              ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
              ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
              ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
              ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
              ! is accurate at second order and thus contains significantly less numerical noise.
              ! compute pressure on the fluid/porous medium edge
              if(gravitoacoustic(ispec)) then
                dxd = - potential_gravitoacoustic(iglob)
              else
                dxd = - potential_acoustic(iglob)
              endif
            else
              dxd = pressure_element(i,j) ! == potential_dot_dot_acoustic(iglob) (or potential_dot_dot_gravitoacoustic(iglob))
            endif
            dzd = ZERO

          else if((acoustic(ispec) .or. gravitoacoustic(ispec)) .and.  seismotype /= 6) then

            dxd = vector_field_element(1,i,j)
            dzd = vector_field_element(3,i,j)

          else if(acoustic(ispec) .and. seismotype == 6) then

            dxd = potential_acoustic(iglob)
            dzd = ZERO

          else if(gravitoacoustic(ispec) .and. seismotype == 6) then

            dxd = potential_acoustic(iglob)
            dzd = ZERO

          else if(seismotype == 1) then

            if(poroelastic(ispec)) then
              dxd = displs_poroelastic(1,iglob)
              dzd = displs_poroelastic(2,iglob)
            else if(elastic(ispec)) then
              dxd = displ_elastic(1,iglob)
              dyd = displ_elastic(2,iglob)
              dzd = displ_elastic(3,iglob)
            endif

          else if(seismotype == 2) then

            if(poroelastic(ispec)) then
              dxd = velocs_poroelastic(1,iglob)
              dzd = velocs_poroelastic(2,iglob)
            else if(elastic(ispec)) then
              dxd = veloc_elastic(1,iglob)
              dyd = veloc_elastic(2,iglob)
              dzd = veloc_elastic(3,iglob)
            endif

          else if(seismotype == 3) then

            if(poroelastic(ispec)) then
              dxd = accels_poroelastic(1,iglob)
              dzd = accels_poroelastic(2,iglob)
            else if(elastic(ispec)) then
              dxd = accel_elastic(1,iglob)
              dyd = accel_elastic(2,iglob)
              dzd = accel_elastic(3,iglob)
            endif

          else if(seismotype == 5) then

            if(poroelastic(ispec)) then
              dxd = displs_poroelastic(1,iglob)
              dzd = displs_poroelastic(2,iglob)
            else if(elastic(ispec)) then
              dxd = displ_elastic(1,iglob)
              dzd = displ_elastic(2,iglob)
            endif
            dcurld = curl_element(i,j)

          endif

          ! compute interpolated field
          valux = valux + dxd*hlagrange
          if(elastic(ispec))  valuy = valuy + dyd*hlagrange
          valuz = valuz + dzd*hlagrange
          valcurl = valcurl + dcurld*hlagrange

        enddo
      enddo

      ! check for edge effects
      if(seismo_current < 1 .or. seismo_current > NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos) &
        stop 'error: seismo_current out of bounds in recording of seismograms'

      ! rotate seismogram components if needed, except if recording pressure, which is a scalar
      if(seismotype /= 4 .and. seismotype /= 6) then
        if(p_sv) then
          sisux(seismo_current,irecloc) =   cosrot_irec(irecloc)*valux + sinrot_irec(irecloc)*valuz
          sisuz(seismo_current,irecloc) = - sinrot_irec(irecloc)*valux + cosrot_irec(irecloc)*valuz
        else
          sisux(seismo_current,irecloc) = valuy
          sisuz(seismo_current,irecloc) = ZERO
        endif
      else
        sisux(seismo_current,irecloc) = valux
        sisuz(seismo_current,irecloc) = ZERO
      endif
      siscurl(seismo_current,irecloc) = valcurl

    enddo

  end subroutine write_seismograms



!================================================================

  subroutine write_seismograms_to_file(x_source,z_source)

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par, only : sisux,sisuz,siscurl,station_name,network_name, &
                          NSTEP,which_proc_receiver,nrec,myrank,deltat,seismotype,t0, &
                          NSTEP_BETWEEN_OUTPUT_SEISMOS,seismo_offset,seismo_current,p_sv, &
                          SU_FORMAT,save_ASCII_seismograms, &
                          save_binary_seismograms_single,save_binary_seismograms_double,subsamp_seismos

! uncomment this to save the ASCII *.sem* seismograms in binary instead, to save disk space and/or writing time
! we could/should move this flag to DATA/Par_file one day.
!
! #define PAUL_SAVE_ASCII_IN_BINARY

  include "constants.h"

  logical :: save_binary_seismograms
  integer irec,length_station_name,length_network_name,iorientation,isample,number_of_components

  character(len=4) channel
  character(len=1) component
  character(len=4) suffix
  character(len=150) sisname

! to write seismograms in single precision SEP and double precision binary format
  double precision, dimension(:,:,:), allocatable :: buffer_binary

! scaling factor for Seismic Unix xsu dislay
  double precision, parameter :: FACTORXSU = 1.d0


  integer  :: irecloc

  double precision :: x_source,z_source

#ifdef USE_MPI
  integer  :: ierror
#endif

!----

! see if we need to save any seismogram in binary format
  save_binary_seismograms = save_binary_seismograms_single .or. save_binary_seismograms_double

  if(SU_FORMAT .and. .not. save_binary_seismograms_single) &
     stop 'error: SU_FORMAT seismograms are single precision and thus require save_binary_seismograms_single'

! write seismograms in ASCII format

! save displacement, velocity, acceleration or pressure
  if(seismotype == 1) then
    component = 'd'
  else if(seismotype == 2) then
    component = 'v'
  else if(seismotype == 3) then
    component = 'a'
  else if(seismotype == 4 .or. seismotype == 6) then
    component = 'p'
  else if(seismotype == 5) then
    component = 'c'
  else
    call exit_MPI('wrong component to save for seismograms')
  endif


! only one seismogram if pressures or SH (membrane) waves
  if(seismotype == 4 .or. seismotype == 6 .or. .not. p_sv) then
     number_of_components = 1
  else if(seismotype == 5) then
     number_of_components = NDIM+1
  else
     number_of_components = NDIM
  endif


! filename extension
  if (.not. SU_FORMAT) then
    suffix = '.bin'
  else
    suffix = '.su'
  endif

  allocate(buffer_binary(NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,nrec,number_of_components))

  if (save_binary_seismograms .and. myrank == 0 .and. seismo_offset == 0) then

! delete old files
     open(unit=12,file='OUTPUT_FILES/Ux_file_single'//suffix,status='unknown')
     close(12,status='delete')

     open(unit=12,file='OUTPUT_FILES/Ux_file_double'//suffix,status='unknown')
     close(12,status='delete')

     open(unit=12,file='OUTPUT_FILES/Uy_file_single'//suffix,status='unknown')
     close(12,status='delete')

     open(unit=12,file='OUTPUT_FILES/Uy_file_double'//suffix,status='unknown')
     close(12,status='delete')

     open(unit=12,file='OUTPUT_FILES/Uz_file_single'//suffix,status='unknown')
     close(12,status='delete')

     open(unit=12,file='OUTPUT_FILES/Uz_file_double'//suffix,status='unknown')
     close(12,status='delete')

     open(unit=12,file='OUTPUT_FILES/Up_file_single'//suffix,status='unknown')
     close(12,status='delete')

     open(unit=12,file='OUTPUT_FILES/Up_file_double'//suffix,status='unknown')
     close(12,status='delete')

     open(unit=12,file='OUTPUT_FILES/Uc_file_single'//suffix,status='unknown')
     close(12,status='delete')

     open(unit=12,file='OUTPUT_FILES/Uc_file_double'//suffix,status='unknown')
     close(12,status='delete')

  endif


  if (save_binary_seismograms .and. myrank == 0) then

     ! write the new files
     if(save_binary_seismograms_single) then
     if(seismotype == 4 .or. seismotype == 6) then
        open(unit=12,file='OUTPUT_FILES/Up_file_single'//suffix,status='unknown',access='direct',recl=4)
     else if(.not.p_sv) then
        open(unit=12,file='OUTPUT_FILES/Uy_file_single'//suffix,status='unknown',access='direct',recl=4)
     else
        open(unit=12,file='OUTPUT_FILES/Ux_file_single'//suffix,status='unknown',access='direct',recl=4)
     endif
     endif

     if(save_binary_seismograms_double) then
     if(seismotype == 4 .or. seismotype == 6) then
     else if(.not.p_sv) then
        open(unit=13,file='OUTPUT_FILES/Uz_file_double'//suffix,status='unknown',access='direct',recl=8)
     else
        open(unit=13,file='OUTPUT_FILES/Ux_file_double'//suffix,status='unknown',access='direct',recl=8)
     endif
     endif

     ! no Z component seismogram if pressure
     if(seismotype /= 4 .and. seismotype /= 6 .and. p_sv) then
       if(save_binary_seismograms_single) &
        open(unit=14,file='OUTPUT_FILES/Uz_file_single'//suffix,status='unknown',access='direct',recl=4)
       if(save_binary_seismograms_double) &
        open(unit=15,file='OUTPUT_FILES/Uz_file_double'//suffix,status='unknown',access='direct',recl=8)
     endif

     ! curl output
     if(seismotype == 5) then
       if(save_binary_seismograms_single) &
        open(unit=16,file='OUTPUT_FILES/Uc_file_single'//suffix,status='unknown',access='direct',recl=4)
       if(save_binary_seismograms_double) &
        open(unit=17,file='OUTPUT_FILES/Uc_file_double'//suffix,status='unknown',access='direct',recl=8)
     endif

  endif


  irecloc = 0
  do irec = 1,nrec

     if ( myrank == 0 ) then

        if ( which_proc_receiver(irec) == myrank ) then
           irecloc = irecloc + 1
           buffer_binary(:,irec,1) = sisux(:,irecloc)
           if ( number_of_components == 2 ) then
              buffer_binary(:,irec,2) = sisuz(:,irecloc)
           else if ( number_of_components == 3 ) then
              buffer_binary(:,irec,2) = sisuz(:,irecloc)
              buffer_binary(:,irec,3) = siscurl(:,irecloc)
           endif

#ifdef USE_MPI
        else
           call MPI_RECV(buffer_binary(1,irec,1),NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,MPI_DOUBLE_PRECISION,&
                which_proc_receiver(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)
           if ( number_of_components == 2 ) then
              call MPI_RECV(buffer_binary(1,irec,2),NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,MPI_DOUBLE_PRECISION,&
                   which_proc_receiver(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)
           endif
           if ( number_of_components == 3 ) then
              call MPI_RECV(buffer_binary(1,irec,2),NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,MPI_DOUBLE_PRECISION,&
                   which_proc_receiver(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)
              call MPI_RECV(buffer_binary(1,irec,3),NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,MPI_DOUBLE_PRECISION,&
                   which_proc_receiver(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)
           endif
#endif
        endif

    endif

  enddo

  irecloc = 0
  do irec = 1,nrec

     if ( myrank == 0 ) then

        if(.not. SU_FORMAT) then

          if(save_ASCII_seismograms) then

          ! write trace
          do iorientation = 1,number_of_components

             if(iorientation == 1) then
                channel = 'BXX'
             else if(iorientation == 2) then
                channel = 'BXZ'
             else if(iorientation == 3) then
                channel = 'cur'
             else
                call exit_MPI('incorrect channel value')
             endif

             ! in case of pressure, use different abbreviation
             if(seismotype == 4 .or. seismotype == 6) channel = 'PRE'
             ! in case of SH (membrane) waves, use different abbreviation
             if(.not.p_sv) channel = 'BXY'

             ! create the name of the seismogram file for each slice
             ! file name includes the name of the station, the network and the component
             length_station_name = len_trim(station_name(irec))
             length_network_name = len_trim(network_name(irec))

             ! check that length conforms to standard
             if(length_station_name < 1 .or. length_station_name > MAX_LENGTH_STATION_NAME) then
               call exit_MPI('wrong length of station name')
            endif
             if(length_network_name < 1 .or. length_network_name > MAX_LENGTH_NETWORK_NAME) then
               call exit_MPI('wrong length of network name')
            endif

             write(sisname,"('OUTPUT_FILES/',a,'.',a,'.',a3,'.sem',a1)") &
                  network_name(irec)(1:length_network_name),station_name(irec)(1:length_station_name),channel,component

             ! save seismograms in text format with no subsampling.
             ! Because we do not subsample the output, this can result in large files
             ! if the simulation uses many time steps. However, subsampling the output
             ! here would result in a loss of accuracy when one later convolves
             ! the results with the source time function
             if ( seismo_offset == 0 ) then
               open(unit=11,file=sisname(1:len_trim(sisname)),status='unknown')
               close(11,status='delete')
             endif
#ifndef PAUL_SAVE_ASCII_IN_BINARY
             open(unit=11,file=sisname(1:len_trim(sisname)),status='unknown',position='append')
#else
             open(unit=11,file=sisname(1:len_trim(sisname)),status='unknown',form='unformatted',position='append')
#endif

             ! make sure we never write more than the maximum number of time steps
             ! subtract offset of the source to make sure travel time is correct
#ifndef PAUL_SAVE_ASCII_IN_BINARY
             do isample = 1,seismo_current
                 write(11,*) sngl(dble(seismo_offset+isample-1)*deltat - t0),' ', &
                              sngl(buffer_binary(isample,irec,iorientation))
             enddo
#else
                 write(11) sngl(buffer_binary(:,irec,iorientation))
#endif

             close(11)
          enddo

          endif

          ! write binary seismogram
          if(save_binary_seismograms) then
          do isample = 1, seismo_current
            if(save_binary_seismograms_single) &
              write(12,rec=(irec-1)*NSTEP+seismo_offset+isample) sngl(buffer_binary(isample,irec,1))
              if(save_binary_seismograms_double) &
              write(13,rec=(irec-1)*NSTEP+seismo_offset+isample) buffer_binary(isample,irec,1)
            if ( seismotype /= 4 .and. seismotype /= 6 .and. p_sv) then
              if(save_binary_seismograms_single) &
                write(14,rec=(irec-1)*NSTEP+seismo_offset+isample) sngl(buffer_binary(isample,irec,2))
              if(save_binary_seismograms_double) &
                write(15,rec=(irec-1)*NSTEP+seismo_offset+isample) buffer_binary(isample,irec,2)
            endif
            if ( seismotype == 5 ) then
              if(save_binary_seismograms_single) &
                write(16,rec=(irec-1)*NSTEP+seismo_offset+isample) sngl(buffer_binary(isample,irec,3))
              if(save_binary_seismograms_double) &
                write(17,rec=(irec-1)*NSTEP+seismo_offset+isample) buffer_binary(isample,irec,3)
            endif
          enddo
          endif

        else ! if SU_FORMAT

            call write_output_SU(x_source,z_source,irec,buffer_binary,number_of_components)

        endif

#ifdef USE_MPI
     else
        if (which_proc_receiver(irec) == myrank) then
           irecloc = irecloc + 1
           call MPI_SEND(sisux(1,irecloc),NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,MPI_DOUBLE_PRECISION,0,irec, &
                           MPI_COMM_WORLD,ierror)
           if (number_of_components >= 2) then
              call MPI_SEND(sisuz(1,irecloc),NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,MPI_DOUBLE_PRECISION,0,irec, &
                           MPI_COMM_WORLD,ierror)
           endif
           if (number_of_components == 3) then
              call MPI_SEND(siscurl(1,irecloc),NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,MPI_DOUBLE_PRECISION,0,irec, &
                           MPI_COMM_WORLD,ierror)
           endif
        endif
#endif

     endif

  enddo

  if(save_binary_seismograms_single) close(12)
  if(save_binary_seismograms_double) close(13)
  if (seismotype /= 4 .and. seismotype /= 6 .and. p_sv) then
    if(save_binary_seismograms_single) close(14)
    if(save_binary_seismograms_double) close(15)
  endif
  if (seismotype == 5) then
    if(save_binary_seismograms_single) close(16)
    if(save_binary_seismograms_double) close(17)
  endif

!----

  deallocate(buffer_binary)

  end subroutine write_seismograms_to_file

