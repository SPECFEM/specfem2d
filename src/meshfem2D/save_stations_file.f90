!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

  subroutine save_stations_file(nreceiversets,nrec_line,xdeb,zdeb,xfin,zfin,record_at_surface_same_vertical, &
                                xinterface_top,zinterface_top,coefs_interface_top, &
                                npoints_interface_top,max_npoints_interface)

  use constants, only: IOUT,IMAIN,MAX_STRING_LEN,mygroup,IN_DATA_FILES
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS

  implicit none

  integer :: nreceiversets
  integer, dimension(nreceiversets) :: nrec_line
  double precision, dimension(nreceiversets) :: xdeb,zdeb,xfin,zfin
  logical, dimension(nreceiversets) :: record_at_surface_same_vertical

  integer :: max_npoints_interface
  double precision, dimension(max_npoints_interface) :: xinterface_top, &
    zinterface_top,coefs_interface_top
  integer :: npoints_interface_top

  !local parameters
  integer :: ireceiverlines,irec,irec_global_number,ios
  integer :: nrec_total
  double precision :: xrec,zrec
  double precision, external :: value_spline
  character(len=MAX_STRING_LEN) :: stations_filename,path_to_add

  stations_filename = trim(IN_DATA_FILES)//'STATIONS'

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    stations_filename = path_to_add(1:len_trim(path_to_add))//stations_filename(1:len_trim(stations_filename))
  endif

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) 'writing the '//trim(stations_filename)//' file'
  write(IMAIN,*)

  ! total number of receivers in all the receiver lines
  nrec_total = sum(nrec_line(:))

  write(IMAIN,*)
  write(IMAIN,*) 'There are ',nrec_total,' receivers'

  write(IMAIN,*)
  write(IMAIN,*) 'Target positions (x,z) of the ',nrec_total,' receivers'
  write(IMAIN,*)

  open(unit=IOUT,file=trim(stations_filename),status='unknown',iostat=ios)
  if (ios /= 0 ) call stop_the_code('error saving STATIONS file')

  irec_global_number = 0

  ! loop on all the receiver lines
  do ireceiverlines = 1,nreceiversets

    ! loop on all the receivers of this receiver line
    do irec = 1,nrec_line(ireceiverlines)

       ! compute global receiver number
       irec_global_number = irec_global_number + 1

       ! compute coordinates of the receiver
       if (nrec_line(ireceiverlines) > 1) then
          xrec = xdeb(ireceiverlines) + dble(irec-1)*(xfin(ireceiverlines) &
                                  -xdeb(ireceiverlines))/dble(nrec_line(ireceiverlines)-1)
          zrec = zdeb(ireceiverlines) + dble(irec-1)*(zfin(ireceiverlines) &
                                  -zdeb(ireceiverlines))/dble(nrec_line(ireceiverlines)-1)
       else
          xrec = xdeb(ireceiverlines)
          zrec = zdeb(ireceiverlines)
       endif

       ! modify position of receiver if we must record exactly at the surface
       if (record_at_surface_same_vertical(ireceiverlines)) &
            zrec = value_spline(xrec,xinterface_top,zinterface_top,coefs_interface_top,npoints_interface_top)

       ! display position of the receiver
       write(IMAIN,*) 'Receiver ',irec_global_number,' = ',xrec,zrec

       write(IOUT,"('S',i4.4,'    AA ',f20.7,1x,f20.7,'       0.0         0.0')") irec_global_number,xrec,zrec

    enddo
  enddo

  close(IOUT)

  end subroutine save_stations_file

