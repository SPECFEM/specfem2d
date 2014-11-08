
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

  subroutine write_seismograms(sisux,sisuz,siscurl,station_name,network_name, &
      NSTEP,nrecloc,which_proc_receiver,nrec,myrank,deltat,seismotype,st_xval,t0, &
      NSTEP_BETWEEN_OUTPUT_SEISMOS,seismo_offset,seismo_current,p_sv, &
      st_zval,x_source,z_source,SU_FORMAT,save_ASCII_seismograms, &
      save_binary_seismograms_single,save_binary_seismograms_double,subsamp_seismos)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

! uncomment this to save the ASCII *.sem* seismograms in binary instead, to save disk space and/or writing time
! we could/should move this flag to DATA/Par_file one day.
!
! #define PAUL_SAVE_ASCII_IN_BINARY

  include "constants.h"

  integer :: nrec,NSTEP,seismotype,subsamp_seismos
  integer :: NSTEP_BETWEEN_OUTPUT_SEISMOS,seismo_offset,seismo_current
  double precision :: t0,deltat

! output seismograms in Seismic Unix format (adjoint traces will be read in the same format)
  logical :: SU_FORMAT

  logical :: p_sv,save_ASCII_seismograms,save_binary_seismograms,save_binary_seismograms_single,save_binary_seismograms_double

  integer, intent(in) :: nrecloc,myrank
  integer, dimension(nrec),intent(in) :: which_proc_receiver

  double precision, dimension(NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,nrecloc), intent(in) :: sisux,sisuz,siscurl

  double precision :: st_xval(nrec)

  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name

  integer irec,length_station_name,length_network_name,iorientation,isample,number_of_components

  character(len=4) chn
  character(len=1) component
  character(len=150) sisname

! to write seismograms in single precision SEP and double precision binary format
  double precision, dimension(:,:), allocatable :: buffer_binary

! scaling factor for Seismic Unix xsu dislay
  double precision, parameter :: FACTORXSU = 1.d0


  integer  :: irecloc

!<SU_FORMAT
  double precision :: st_zval(nrec),x_source,z_source
  integer(kind=2) :: header2(2)
!>SU_FORMAT

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

  allocate(buffer_binary(NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,number_of_components))

  if (save_binary_seismograms .and. myrank == 0 .and. seismo_offset == 0) then

! delete the old files
     open(unit=12,file='OUTPUT_FILES/Ux_file_single.bin',status='unknown')
     close(12,status='delete')

     open(unit=12,file='OUTPUT_FILES/Ux_file_double.bin',status='unknown')
     close(12,status='delete')

     open(unit=12,file='OUTPUT_FILES/pressure_file_single.bin',status='unknown')
     close(12,status='delete')

     open(unit=12,file='OUTPUT_FILES/pressure_file_double.bin',status='unknown')
     close(12,status='delete')

     open(unit=12,file='OUTPUT_FILES/Uz_file_single.bin',status='unknown')
     close(12,status='delete')

     open(unit=12,file='OUTPUT_FILES/Uz_file_double.bin',status='unknown')
     close(12,status='delete')

     open(unit=12,file='OUTPUT_FILES/Curl_file_single.bin',status='unknown')
     close(12,status='delete')

     open(unit=12,file='OUTPUT_FILES/Curl_file_double.bin',status='unknown')
     close(12,status='delete')

  endif

  if (save_binary_seismograms .and. myrank == 0) then

! write the new files
     if(save_binary_seismograms_single) then
     if(seismotype == 4 .or. seismotype == 6) then
        open(unit=12,file='OUTPUT_FILES/pressure_file_single.bin',status='unknown',access='direct',recl=4)
     else if(.not.p_sv) then
        open(unit=12,file='OUTPUT_FILES/Uy_file_single.bin',status='unknown',access='direct',recl=4)
     else
        open(unit=12,file='OUTPUT_FILES/Ux_file_single.bin',status='unknown',access='direct',recl=4)
     endif
     endif

     if(save_binary_seismograms_double) then
     if(seismotype == 4 .or. seismotype == 6) then
        open(unit=13,file='OUTPUT_FILES/pressure_file_double.bin',status='unknown',access='direct',recl=8)
     else if(.not.p_sv) then
        open(unit=13,file='OUTPUT_FILES/Uz_file_double.bin',status='unknown',access='direct',recl=8)
     else
        open(unit=13,file='OUTPUT_FILES/Ux_file_double.bin',status='unknown',access='direct',recl=8)
     endif
     endif

! no Z component seismogram if pressure
     if(seismotype /= 4 .and. seismotype /= 6 .and. p_sv) then
       if(save_binary_seismograms_single) &
        open(unit=14,file='OUTPUT_FILES/Uz_file_single.bin',status='unknown',access='direct',recl=4)
       if(save_binary_seismograms_double) &
        open(unit=15,file='OUTPUT_FILES/Uz_file_double.bin',status='unknown',access='direct',recl=8)
     endif

! curl output
     if(seismotype == 5) then
       if(save_binary_seismograms_single) &
        open(unit=16,file='OUTPUT_FILES/Curl_file_single.bin',status='unknown',access='direct',recl=4)
       if(save_binary_seismograms_double) &
        open(unit=17,file='OUTPUT_FILES/Curl_file_double.bin',status='unknown',access='direct',recl=8)
     endif

  endif


  irecloc = 0
  do irec = 1,nrec

     if ( myrank == 0 ) then

        if ( which_proc_receiver(irec) == myrank ) then
           irecloc = irecloc + 1
           buffer_binary(:,1) = sisux(:,irecloc)
           if ( number_of_components == 2 ) then
              buffer_binary(:,2) = sisuz(:,irecloc)
           else if ( number_of_components == 3 ) then
              buffer_binary(:,2) = sisuz(:,irecloc)
              buffer_binary(:,3) = siscurl(:,irecloc)
           endif

#ifdef USE_MPI
        else
           call MPI_RECV(buffer_binary(1,1),NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,MPI_DOUBLE_PRECISION,&
                which_proc_receiver(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)
           if ( number_of_components == 2 ) then
              call MPI_RECV(buffer_binary(1,2),NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,MPI_DOUBLE_PRECISION,&
                   which_proc_receiver(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)
           endif
           if ( number_of_components == 3 ) then
              call MPI_RECV(buffer_binary(1,2),NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,MPI_DOUBLE_PRECISION,&
                   which_proc_receiver(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)
              call MPI_RECV(buffer_binary(1,3),NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,MPI_DOUBLE_PRECISION,&
                   which_proc_receiver(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)
           endif
#endif
        endif

        if(.not. SU_FORMAT) then

          if(save_ASCII_seismograms) then

          ! write trace
          do iorientation = 1,number_of_components

             if(iorientation == 1) then
                chn = 'BXX'
             else if(iorientation == 2) then
                chn = 'BXZ'
             else if(iorientation == 3) then
                chn = 'cur'
             else
                call exit_MPI('incorrect channel value')
             endif

             ! in case of pressure, use different abbreviation
             if(seismotype == 4 .or. seismotype == 6) chn = 'PRE'
             ! in case of SH (membrane) waves, use different abbreviation
             if(.not.p_sv) chn = 'BXY'

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
                  network_name(irec)(1:length_network_name),station_name(irec)(1:length_station_name),chn,component

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
                              sngl(buffer_binary(isample,iorientation))
             enddo
#else
                 write(11) sngl(buffer_binary(:,iorientation))
#endif

             close(11)
          enddo

          endif

          ! write binary seismogram
          if(save_binary_seismograms) then
          do isample = 1, seismo_current
            if(save_binary_seismograms_single) &
              write(12,rec=(irec-1)*NSTEP+seismo_offset+isample) sngl(buffer_binary(isample,1))
              if(save_binary_seismograms_double) &
              write(13,rec=(irec-1)*NSTEP+seismo_offset+isample) buffer_binary(isample,1)
            if ( seismotype /= 4 .and. seismotype /= 6 .and. p_sv) then
              if(save_binary_seismograms_single) &
                write(14,rec=(irec-1)*NSTEP+seismo_offset+isample) sngl(buffer_binary(isample,2))
              if(save_binary_seismograms_double) &
                write(15,rec=(irec-1)*NSTEP+seismo_offset+isample) buffer_binary(isample,2)
            endif
            if ( seismotype == 5 ) then
              if(save_binary_seismograms_single) &
                write(16,rec=(irec-1)*NSTEP+seismo_offset+isample) sngl(buffer_binary(isample,3))
              if(save_binary_seismograms_double) &
                write(17,rec=(irec-1)*NSTEP+seismo_offset+isample) buffer_binary(isample,3)
            endif
          enddo
          endif

        else ! if SU_FORMAT

          if (seismo_offset==0) then
             ! write SU headers (refer to Seismic Unix for details)
             write(12,rec=(irec-1)*60+(irec-1)*NSTEP+1)  irec                          ! receiver ID
             write(12,rec=(irec-1)*60+(irec-1)*NSTEP+10) NINT(st_xval(irec)-x_source)  ! offset
             write(12,rec=(irec-1)*60+(irec-1)*NSTEP+19) NINT(x_source)                ! source location xs
             write(12,rec=(irec-1)*60+(irec-1)*NSTEP+20) NINT(z_source)                ! source location zs
             write(12,rec=(irec-1)*60+(irec-1)*NSTEP+21) NINT(st_xval(irec))           ! receiver location xr
             write(12,rec=(irec-1)*60+(irec-1)*NSTEP+22) NINT(st_zval(irec))           ! receiver location zr
             if (nrec>1) write(12,rec=(irec-1)*60+(irec-1)*NSTEP+48) SNGL(st_xval(2)-st_xval(1)) ! receiver interval
             header2(1)=0  ! dummy
             header2(2)=int(NSTEP, kind=2)
             write(12,rec=(irec-1)*60+(irec-1)*NSTEP+29) header2
             header2(1)=NINT(deltat*1.0d6, kind=2)  ! deltat (unit: 10^{-6} second)
             header2(2)=0  ! dummy
             write(12,rec=(irec-1)*60+(irec-1)*NSTEP+30) header2
             if ( seismotype /= 4 .and. seismotype /= 6 .and. p_sv) then
                ! headers
                if (seismo_offset==0) then
                   write(14,rec=(irec-1)*60+(irec-1)*NSTEP+1)  irec
                   write(14,rec=(irec-1)*60+(irec-1)*NSTEP+10) NINT(st_xval(irec)-x_source)
                   write(14,rec=(irec-1)*60+(irec-1)*NSTEP+19) NINT(x_source)
                   write(14,rec=(irec-1)*60+(irec-1)*NSTEP+20) NINT(z_source)
                   write(14,rec=(irec-1)*60+(irec-1)*NSTEP+21) NINT(st_xval(irec))
                   write(14,rec=(irec-1)*60+(irec-1)*NSTEP+22) NINT(st_zval(irec))
                   if(nrec>1) write(14,rec=(irec-1)*60+(irec-1)*NSTEP+48) SNGL(st_xval(2)-st_xval(1))
                   header2(1)=0  ! dummy
                   header2(2)=int(NSTEP, kind=2)
                   write(14,rec=(irec-1)*60+(irec-1)*NSTEP+29) header2
                   header2(1)=NINT(deltat*1.0d6, kind=2)
                   header2(2)=0  ! dummy
                   write(14,rec=(irec-1)*60+(irec-1)*NSTEP+30) header2
                endif
             endif
          endif

          ! the "60" in the following corresponds to 240 bytes header (note the reclength is 4 bytes)
          do isample = 1, seismo_current
             write(12,rec=irec*60+(irec-1)*NSTEP+seismo_offset+isample) sngl(buffer_binary(isample,1))
             if ( seismotype /= 4 .and. seismotype /= 6 .and. p_sv) then
                write(14,rec=irec*60+(irec-1)*NSTEP+seismo_offset+isample) sngl(buffer_binary(isample,2))
             endif
          enddo

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

  end subroutine write_seismograms

