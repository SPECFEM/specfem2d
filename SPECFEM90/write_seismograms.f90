
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.1
!                   ------------------------------
!
!                         Dimitri Komatitsch
!          Universite de Pau et des Pays de l'Adour, France
!
!                          (c) December 2004
!
!========================================================================

  subroutine write_seismograms(sisux,sisuz,nt,nrec,deltat,sismostype,iglob_rec,coord,npoin)

! save the seismograms in ASCII format

  implicit none

  include "constants.h"

  integer nt,nrec,sismostype,npoin
  double precision deltat

  integer iglob_rec(nrec)

  double precision sisux(nt,nrec)
  double precision sisuz(nt,nrec)
  double precision coord(NDIME,npoin)

  integer irec,it

  character(len=100) name

! scaling factor for Seismic Unix xsu dislay
  double precision, parameter :: FACTORXSU = 1.d0

! write seismograms in ASCII format

! X component
  do irec=1,nrec
    write(name,221) irec
    open(unit=11,file=name,status='unknown')
    do it=1,nt
      write(11,*) sngl(dble(it-1)*deltat),' ',sngl(sisux(it,irec))
    enddo
    close(11)
  enddo

! Z component
  do irec=1,nrec
    write(name,222) irec
    open(unit=11,file=name,status='unknown')
    do it=1,nt
      write(11,*) sngl(dble(it-1)*deltat),' ',sngl(sisuz(it,irec))
    enddo
    close(11)
  enddo

!----

! write seismograms in single precision binary format

! X component
  open(unit=11,file='Ux_file.bin',status='unknown',access='direct',recl=nt*nrec*4)
  write(11,rec=1) (sngl(sisux(irec,1)),irec=1,nt*nrec)
  close(11)

! Z component
  open(unit=11,file='Uz_file.bin',status='unknown',access='direct',recl=nt*nrec*4)
  write(11,rec=1) (sngl(sisuz(irec,1)),irec=1,nt*nrec)
  close(11)

!----

! ligne de recepteurs pour Xsu
  open(unit=11,file='receiver_line_Xsu_XWindow',status='unknown')

  write(11,110) FACTORXSU,nt,deltat,nrec

  do irec=1,nrec
    write(11,140) coord(1,iglob_rec(irec))
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
  open(unit=11,file='show_receiver_line_Xsu',status='unknown')
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

  110 format('xwigb@xcur=',f8.2,'@n1=',i5,'@d1=',f15.8,'@label1="Time@(s)"@label2="x@(m)"@n2=',i5,'@x2=')

  120 format('sed -e ''1,$s/ //g'' -e ''1,$s/@/ /g'' -e ''1,1p'' -e ''$,$s/Ux/Uz/g'' <tempfile > receiver_line_Xsu_XWindow')

  130 format('sed -e ''1,$s/xwigb/pswigp/g'' ', &
        '-e ''1,$s/Ux_file.bin/Ux_file.bin > uxpoly.ps/g'' ', &
        '-e ''1,$s/Uz_file.bin/Uz_file.bin > uzpoly.ps/g'' receiver_line_Xsu_XWindow > receiver_line_Xsu_postscript')

  140 format(f12.5)

  221 format('Ux_file_',i3.3,'.dat')
  222 format('Uz_file_',i3.3,'.dat')

  end subroutine write_seismograms

