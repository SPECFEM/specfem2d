
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.1
!                   ------------------------------
!
!                         Dimitri Komatitsch
!          Universite de Pau et des Pays de l'Adour, France
!
!                          (c) January 2005
!
!========================================================================

! write seismograms to text files

  subroutine write_seismograms(sisux,sisuz,station_name,network_name,NSTEP,nrec,deltat,sismostype,st_xval,it,t0)

  implicit none

  include "constants.h"

  integer nrec,NSTEP,it,sismostype
  double precision t0,deltat

  double precision,dimension(NSTEP,nrec) :: sisux,sisuz

  double precision st_xval(nrec)

  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name

  integer irec,length_station_name,length_network_name,iorientation,isample

  character(len=4) chn
  character(len=1) component
  character(len=150) sisname

! scaling factor for Seismic Unix xsu dislay
  double precision, parameter :: FACTORXSU = 1.d0

! write seismograms in ASCII format

! save displacement, velocity or acceleration
  if(sismostype == 1) then
    component = 'd'
  else if(sismostype == 2) then
    component = 'v'
  else if(sismostype == 3) then
    component = 'a'
  else
    stop 'wrong component to save for seismograms'
  endif

  do irec = 1,nrec

    do iorientation = 1,NDIM

      if(iorientation == 1) then
        chn = 'BHX'
      else if(iorientation == 2) then
        chn = 'BHZ'
      else
        stop 'incorrect channel value'
      endif

! create the name of the seismogram file for each slice
! file name includes the name of the station, the network and the component
      length_station_name = len_trim(station_name(irec))
      length_network_name = len_trim(network_name(irec))

! check that length conforms to standard
      if(length_station_name < 1 .or. length_station_name > MAX_LENGTH_STATION_NAME) stop 'wrong length of station name'

      if(length_network_name < 1 .or. length_network_name > MAX_LENGTH_NETWORK_NAME) stop 'wrong length of network name'

      write(sisname,"('OUTPUT_FILES/',a,'.',a,'.',a3,'.sem',a1)") station_name(irec)(1:length_station_name),&
           network_name(irec)(1:length_network_name),chn,component

! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function
      open(unit=IOUT,file=sisname(1:len_trim(sisname)),status='unknown')

! make sure we never write more than the maximum number of time steps
! subtract offset of the source to make sure travel time is correct
      do isample = 1,min(it,NSTEP)
        if(iorientation == 1) then
          write(IOUT,*) sngl(dble(isample-1)*deltat - t0),' ',sngl(sisux(isample,irec))
        else
          write(IOUT,*) sngl(dble(isample-1)*deltat - t0),' ',sngl(sisuz(isample,irec))
        endif
      enddo

      close(IOUT)

      enddo

  enddo

!----

! write seismograms in single precision binary format

! X component
  open(unit=11,file='OUTPUT_FILES/Ux_file.bin',status='unknown',access='direct',recl=NSTEP*4)
  do irec=1,nrec
    write(11,rec=irec) (sngl(sisux(isample,irec)),isample=1,NSTEP)
  enddo
  close(11)

! Z component
  open(unit=11,file='OUTPUT_FILES/Uz_file.bin',status='unknown',access='direct',recl=NSTEP*4)
  do irec=1,nrec
    write(11,rec=irec) (sngl(sisuz(isample,irec)),isample=1,NSTEP)
  enddo
  close(11)

!----

! ligne de recepteurs pour Xsu
  open(unit=11,file='OUTPUT_FILES/receiver_line_Xsu_XWindow',status='unknown')

! subtract t0 from seismograms to get correct zero time
  write(11,110) FACTORXSU,NSTEP,deltat,-t0,nrec

  do irec=1,nrec
    write(11,140) st_xval(irec)
    if(irec < nrec) write(11,*) ','
  enddo

  if(sismostype == 1) then
    write(11,*) '@title="Ux@displacement@component"@<@Ux_file.bin'
  else if(sismostype == 2) then
    write(11,*) '@title="Ux@velocity@component"@<@Ux_file.bin'
  else
    write(11,*) '@title="Ux@acceleration@component"@<@Ux_file.bin'
  endif

  close(11)

! script de visualisation
  open(unit=11,file='OUTPUT_FILES/show_receiver_line_Xsu',status='unknown')
  write(11,100)
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

! formats
  100 format('#!/bin/csh -f')

  110 format('xwigb@xcur=',f8.2,'@n1=',i5,'@d1=',f15.8,'@f1=',f15.8,'@label1="Time@(s)"@label2="x@(m)"@n2=',i5,'@x2=')

  120 format('sed -e ''1,$s/ //g'' -e ''1,$s/@/ /g'' -e ''1,1p'' -e ''$,$s/Ux/Uz/g'' <tempfile > receiver_line_Xsu_XWindow')

  130 format('sed -e ''1,$s/xwigb/pswigp/g'' ', &
        '-e ''1,$s/Ux_file.bin/Ux_file.bin > uxpoly.ps/g'' ', &
        '-e ''1,$s/Uz_file.bin/Uz_file.bin > uzpoly.ps/g'' receiver_line_Xsu_XWindow > receiver_line_Xsu_postscript')

  140 format(f12.5)

  end subroutine write_seismograms

