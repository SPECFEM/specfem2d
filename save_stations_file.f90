
!========================================================================
!
!                   S P E C F E M 2 D  Version 6.1
!                   ------------------------------
!
! Copyright Universite de Pau, CNRS and INRIA, France,
! and Princeton University / California Institute of Technology, USA.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
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

  subroutine save_stations_file(nreceiverlines,nrec,xdeb,zdeb,xfin,zfin,enreg_surf, &
                            xinterface_top,zinterface_top,coefs_interface_top, &
                            npoints_interface_top,max_npoints_interface)

  implicit none

  integer :: nreceiverlines
  integer, dimension(nreceiverlines) :: nrec
  double precision, dimension(nreceiverlines) :: xdeb,zdeb,xfin,zfin
  logical, dimension(nreceiverlines) :: enreg_surf

  integer :: max_npoints_interface
  double precision, dimension(max_npoints_interface) :: xinterface_top, &
    zinterface_top,coefs_interface_top
  integer :: npoints_interface_top

  !local parameters
  integer :: ireceiverlines,irec,irec_global_number,ios
  integer :: nrec_total
  double precision :: xrec,zrec
  double precision, external :: value_spline

  print *
  print *,'writing the DATA/STATIONS_target file'
  print *

  ! total number of receivers in all the receiver lines
  nrec_total = sum(nrec)

  print *
  print *,'There are ',nrec_total,' receivers'

  print *
  print *,'Position (x,z) of the ',nrec_total,' receivers'
  print *

  open(unit=15,file='DATA/STATIONS_target',status='unknown',iostat=ios)
  if( ios /= 0 ) stop 'error saving STATIONS file'

  irec_global_number = 0

  ! loop on all the receiver lines
  do ireceiverlines = 1,nreceiverlines

    ! loop on all the receivers of this receiver line
    do irec = 1,nrec(ireceiverlines)

       ! compute global receiver number
       irec_global_number = irec_global_number + 1

       ! compute coordinates of the receiver
       if(nrec(ireceiverlines) > 1) then
          xrec = xdeb(ireceiverlines) + dble(irec-1)*(xfin(ireceiverlines) &
                                  -xdeb(ireceiverlines))/dble(nrec(ireceiverlines)-1)
          zrec = zdeb(ireceiverlines) + dble(irec-1)*(zfin(ireceiverlines) &
                                  -zdeb(ireceiverlines))/dble(nrec(ireceiverlines)-1)
       else
          xrec = xdeb(ireceiverlines)
          zrec = zdeb(ireceiverlines)
       endif

       ! modify position of receiver if we must record exactly at the surface
       if(enreg_surf(ireceiverlines)) &
            zrec = value_spline(xrec,xinterface_top,zinterface_top, &
                            coefs_interface_top,npoints_interface_top)

       ! display position of the receiver
       print *,'Receiver ',irec_global_number,' = ',xrec,zrec

       write(15,"('S',i4.4,'    AA ',f20.7,1x,f20.7,'       0.0         0.0')") irec_global_number,xrec,zrec

    enddo
  enddo

  close(15)

  end subroutine save_stations_file

