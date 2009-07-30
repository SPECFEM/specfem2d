
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
! Copyright Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
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
      NTSTEP_BETWEEN_OUTPUT_SEISMO,seismo_offset,seismo_current &
      )

  implicit none

  include "constants.h"
#ifdef USE_MPI
  include "mpif.h"
#endif

  integer :: nrec,NSTEP,seismotype
  integer :: NTSTEP_BETWEEN_OUTPUT_SEISMO,seismo_offset,seismo_current
  double precision :: t0,deltat

  integer, intent(in) :: nrecloc,myrank
  integer, dimension(nrec),intent(in) :: which_proc_receiver

  double precision, dimension(NTSTEP_BETWEEN_OUTPUT_SEISMO,nrecloc), intent(in) :: sisux,sisuz,siscurl

  double precision st_xval(nrec)

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

#ifdef USE_MPI
  integer  :: ierror
  integer, dimension(MPI_STATUS_SIZE)  :: status
#endif

!----

! write seismograms in ASCII format

! save displacement, velocity, acceleration or pressure
  if(seismotype == 1) then
    component = 'd'
  else if(seismotype == 2) then
    component = 'v'
  else if(seismotype == 3) then
    component = 'a'
  else if(seismotype == 4) then
    component = 'p'
  else if(seismotype == 5) then
    component = 'c'
  else
    call exit_MPI('wrong component to save for seismograms')
  endif


! only one seismogram if pressurs
  if(seismotype == 4) then
     number_of_components = 1
  else if(seismotype == 5) then
     number_of_components = NDIM+1
  else
     number_of_components = NDIM
  endif

  allocate(buffer_binary(NTSTEP_BETWEEN_OUTPUT_SEISMO,number_of_components))


  if ( myrank == 0 .and. seismo_offset == 0 ) then

! delete the old files
     open(unit=11,file='OUTPUT_FILES/Ux_file_single.bin',status='unknown')
     close(11,status='delete')

     open(unit=11,file='OUTPUT_FILES/Ux_file_double.bin',status='unknown')
     close(11,status='delete')

     open(unit=11,file='OUTPUT_FILES/pressure_file_single.bin',status='unknown')
     close(11,status='delete')

     open(unit=11,file='OUTPUT_FILES/pressure_file_double.bin',status='unknown')
     close(11,status='delete')

     open(unit=11,file='OUTPUT_FILES/Uz_file_single.bin',status='unknown')
     close(11,status='delete')

     open(unit=11,file='OUTPUT_FILES/Uz_file_double.bin',status='unknown')
     close(11,status='delete')

     open(unit=11,file='OUTPUT_FILES/Curl_file_single.bin',status='unknown')
     close(11,status='delete')

     open(unit=11,file='OUTPUT_FILES/Curl_file_double.bin',status='unknown')
     close(11,status='delete')

   endif

   if ( myrank == 0 ) then

! write the new files
     if(seismotype == 4) then
        open(unit=12,file='OUTPUT_FILES/pressure_file_single.bin',status='unknown',access='direct',recl=4)
     else
        open(unit=12,file='OUTPUT_FILES/Ux_file_single.bin',status='unknown',access='direct',recl=4)
     endif

     if(seismotype == 4) then
        open(unit=13,file='OUTPUT_FILES/pressure_file_double.bin',status='unknown',access='direct',recl=8)
     else
        open(unit=13,file='OUTPUT_FILES/Ux_file_double.bin',status='unknown',access='direct',recl=8)
     endif

! no Z component seismogram if pressure
     if(seismotype /= 4) then
        open(unit=14,file='OUTPUT_FILES/Uz_file_single.bin',status='unknown',access='direct',recl=4)
        open(unit=15,file='OUTPUT_FILES/Uz_file_double.bin',status='unknown',access='direct',recl=8)

     end if

! curl output
     if(seismotype == 5) then
        open(unit=16,file='OUTPUT_FILES/Curl_file_single.bin',status='unknown',access='direct',recl=4)
        open(unit=17,file='OUTPUT_FILES/Curl_file_double.bin',status='unknown',access='direct',recl=8)

     end if

  end if


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
           end if

#ifdef USE_MPI
        else
           call MPI_RECV(buffer_binary(1,1),NTSTEP_BETWEEN_OUTPUT_SEISMO,MPI_DOUBLE_PRECISION,&
                which_proc_receiver(irec),irec,MPI_COMM_WORLD,status,ierror)
           if ( number_of_components == 2 ) then
              call MPI_RECV(buffer_binary(1,2),NTSTEP_BETWEEN_OUTPUT_SEISMO,MPI_DOUBLE_PRECISION,&
                   which_proc_receiver(irec),irec,MPI_COMM_WORLD,status,ierror)
           end if
           if ( number_of_components == 3 ) then
              call MPI_RECV(buffer_binary(1,2),NTSTEP_BETWEEN_OUTPUT_SEISMO,MPI_DOUBLE_PRECISION,&
                   which_proc_receiver(irec),irec,MPI_COMM_WORLD,status,ierror)
              call MPI_RECV(buffer_binary(1,3),NTSTEP_BETWEEN_OUTPUT_SEISMO,MPI_DOUBLE_PRECISION,&
                   which_proc_receiver(irec),irec,MPI_COMM_WORLD,status,ierror)
           end if


#endif
        end if

! write trace
        do iorientation = 1,number_of_components

           if(iorientation == 1) then
              chn = 'BHX'
           else if(iorientation == 2) then
              chn = 'BHZ'
           else if(iorientation == 3) then
              chn = 'cur'
           else
              call exit_MPI('incorrect channel value')
           endif

           ! in case of pressure, use different abbreviation
           if(seismotype == 4) chn = 'PRE'

           ! create the name of the seismogram file for each slice
           ! file name includes the name of the station, the network and the component
           length_station_name = len_trim(station_name(irec))
           length_network_name = len_trim(network_name(irec))

           ! check that length conforms to standard
           if(length_station_name < 1 .or. length_station_name > MAX_LENGTH_STATION_NAME) then
             call exit_MPI('wrong length of station name')
          end if
           if(length_network_name < 1 .or. length_network_name > MAX_LENGTH_NETWORK_NAME) then
             call exit_MPI('wrong length of network name')
          end if

           write(sisname,"('OUTPUT_FILES/',a,'.',a,'.',a3,'.sem',a1)") station_name(irec)(1:length_station_name),&
                network_name(irec)(1:length_network_name),chn,component

           ! save seismograms in text format with no subsampling.
           ! Because we do not subsample the output, this can result in large files
           ! if the simulation uses many time steps. However, subsampling the output
           ! here would result in a loss of accuracy when one later convolves
           ! the results with the source time function
           if ( seismo_offset == 0 ) then
             open(unit=11,file=sisname(1:len_trim(sisname)),status='unknown')
             close(11,status='delete')
           endif
           open(unit=11,file=sisname(1:len_trim(sisname)),status='unknown',position='append')

           ! make sure we never write more than the maximum number of time steps
           ! subtract offset of the source to make sure travel time is correct
           do isample = 1,seismo_current
              if(iorientation == 1) then
                 write(11,*) sngl(dble(seismo_offset+isample-1)*deltat - t0),' ',sngl(buffer_binary(isample,iorientation))
              else
                 write(11,*) sngl(dble(seismo_offset+isample-1)*deltat - t0),' ',sngl(buffer_binary(isample,iorientation))
              endif
           enddo

           close(11)
        end do

! write binary seismogram
        do isample = 1, seismo_current
           write(12,rec=(irec-1)*NSTEP+seismo_offset+isample) sngl(buffer_binary(isample,1))
           write(13,rec=(irec-1)*NSTEP+seismo_offset+isample) buffer_binary(isample,1)
        if ( seismotype /= 4 ) then
           write(14,rec=(irec-1)*NSTEP+seismo_offset+isample) sngl(buffer_binary(isample,2))
           write(15,rec=(irec-1)*NSTEP+seismo_offset+isample) buffer_binary(isample,2)
        end if
        if ( seismotype == 5 ) then
           write(16,rec=(irec-1)*NSTEP+seismo_offset+isample) sngl(buffer_binary(isample,3))
           write(17,rec=(irec-1)*NSTEP+seismo_offset+isample) buffer_binary(isample,3)
        end if
        enddo
#ifdef USE_MPI

     else
        if ( which_proc_receiver(irec) == myrank ) then
           irecloc = irecloc + 1
           call MPI_SEND(sisux(1,irecloc),NTSTEP_BETWEEN_OUTPUT_SEISMO,MPI_DOUBLE_PRECISION,0,irec,MPI_COMM_WORLD,ierror)
           if ( number_of_components >= 2 ) then
              call MPI_SEND(sisuz(1,irecloc),NTSTEP_BETWEEN_OUTPUT_SEISMO,MPI_DOUBLE_PRECISION,0,irec,MPI_COMM_WORLD,ierror)
           end if
           if ( number_of_components == 3 ) then
              call MPI_SEND(siscurl(1,irecloc),NTSTEP_BETWEEN_OUTPUT_SEISMO,MPI_DOUBLE_PRECISION,0,irec,MPI_COMM_WORLD,ierror)
           end if
        end if

#endif

     end if

  enddo

  close(12)
  close(13)
  if ( seismotype /= 4 ) then
     close(14)
     close(15)
  end if
  if ( seismotype == 5 ) then
     close(16)
     close(17)
  end if

!----

  deallocate(buffer_binary)

!----
   if ( myrank == 0 ) then

! ligne de recepteurs pour Xsu
  open(unit=11,file='OUTPUT_FILES/receiver_line_Xsu_XWindow',status='unknown')

! subtract t0 from seismograms to get correct zero time
  write(11,110) FACTORXSU,NSTEP,deltat,-t0,nrec

  do irec=1,nrec
    write(11,"(f12.5)") st_xval(irec)
    if(irec < nrec) write(11,*) ','
  enddo

  if(seismotype == 1) then
    write(11,*) '@title="Ux@displacement@component"@<@Ux_file_single.bin'
  else if(seismotype == 2) then
    write(11,*) '@title="Ux@velocity@component"@<@Ux_file_single.bin'
  else
    write(11,*) '@title="Ux@acceleration@component"@<@Ux_file_single.bin'
  endif

  close(11)

! script de visualisation
  open(unit=11,file='OUTPUT_FILES/show_receiver_line_Xsu',status='unknown')
  write(11,"('#!/bin/csh')")
  write(11,*)
  write(11,*) '/bin/rm -f tempfile receiver_line_Xsu_postscript'
  write(11,*) '# concatener toutes les lignes'
  write(11,*) 'tr -d ''\012'' <receiver_line_Xsu_XWindow >tempfile'
  write(11,*) '# remettre fin de ligne'
  write(11,*) 'echo " " >> tempfile'
  write(11,*) '# supprimer espaces, changer arobas, dupliquer'
  write(11,120)
  write(11,*) '/bin/rm -f tempfile'
  write(11,*) '# copier fichier pour sortie postscript'
  write(11,130)
  write(11,*) '/bin/rm -f tempfile'
  write(11,*) 'echo ''rm -f uxpoly.ps uzpoly.ps'' > tempfile'
  write(11,*) 'cat tempfile receiver_line_Xsu_postscript > tempfile2'
  write(11,*) '/bin/mv -f tempfile2 receiver_line_Xsu_postscript'
  write(11,*) '/bin/rm -f tempfile'
  write(11,*) '# executer commande xsu'
  write(11,*) 'sh receiver_line_Xsu_XWindow'
  write(11,*) '/bin/rm -f tempfile tempfile2'
  close(11)

end if

! formats
  110 format('xwigb@xcur=',f8.2,'@n1=',i6,'@d1=',f15.8,'@f1=',f15.8,'@label1="Time@(s)"@label2="x@(m)"@n2=',i6,'@x2=')

  120 format('sed -e ''1,$s/ //g'' -e ''1,$s/@/ /g'' -e ''1,1p'' -e ''$,$s/Ux/Uz/g'' <tempfile > receiver_line_Xsu_XWindow')

  130 format('sed -e ''1,$s/xwigb/pswigp/g'' ', &
        '-e ''1,$s/Ux_file_single.bin/Ux_file_single.bin > uxpoly.ps/g'' ', &
        '-e ''1,$s/Uz_file_single.bin/Uz_file_single.bin > uzpoly.ps/g'' receiver_line_Xsu_XWindow > receiver_line_Xsu_postscript')

  end subroutine write_seismograms

