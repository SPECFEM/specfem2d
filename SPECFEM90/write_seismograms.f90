
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

  subroutine write_seismograms(sisux,sisuz,nt,nrec,deltat)

! save the seismograms in ASCII format

  implicit none

  integer nt,nrec
  double precision deltat

  double precision sisux(nt,nrec)
  double precision sisuz(nt,nrec)

  integer irec,it

  character(len=100) name

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

 221 format('Ux_file_',i3.3,'.dat')
 222 format('Uz_file_',i3.3,'.dat')

  end subroutine write_seismograms

