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
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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

program adj_seismogram

  ! This program cuts a certain portion of the seismograms and convert it
  ! into the adjoint source for generating banana-dougnut kernels
  !    rm -rf xadj_seismogram ; gfortran adj_seismogram.f90 -o xadj_seismogram

  implicit none

  integer, parameter :: NSTEP = 3000
  integer, parameter :: nrec = 1
  double precision, parameter :: t0 = 48.0   ! labeled as 'time delay'
  double precision, parameter :: deltat = 0.06
  double precision, parameter :: EPS = 1.d-40

  integer :: itime,icomp,istart,iend,nlenval,irec,NDIM,NDIMr,adj_comp
  double precision :: timeval,tstart(nrec),tend(nrec)
  character(len=150), dimension(nrec) :: station_name
  double precision, dimension(NSTEP) :: time_window
  double precision :: seism(NSTEP,3),Nnorm,seism_win(NSTEP)
  double precision :: seism_veloc(NSTEP),seism_accel(NSTEP),ft_bar(NSTEP)
  character(len=3) :: compr(2),comp(3)
  character(len=150) :: filename
  integer :: ier

  NDIM = 3
  comp = (/"BXX","BXY","BXZ"/)

  ! number of components
  !NDIMr=2  ! P-SV
  NDIMr=1  ! SH (membrane)
  !compr = (/"BXX","BXZ"/)    ! P-SV
  compr = (/"BXY","tmp"/)  ! SH (membrane)
  ! list of stations
  station_name(1) = 'S0001'

  ! KEY: 'absolute' time interval for window
  tstart(1) = 95.d0 + t0
  tend(1) = 130.d0 + t0

  ! chose the component for the adjoint source (adj_comp = 1:X, 2:Y, 3:Z)
  adj_comp = 2

  do irec = 1,nrec

     do icomp = 1, NDIMr

        filename = 'OUTPUT_FILES/'//'AA.'//trim(station_name(irec))//'.'// compr(icomp) // '.semd'
        open(unit = 10, file = trim(filename), status='old', iostat=ier)
        if (ier /= 0) then
          print *,'Error: could not open trace file ',trim(filename)
          print *,'Please check if file exists...'
          stop 'Error opening trace file'
        endif

        do itime = 1,NSTEP
           read(10,*) timeval , seism(itime,icomp)
        enddo

     enddo

     if (NDIMr==2) then
        seism(:,3) = seism(:,2)
        seism(:,2) = 0.d0
     else
        seism(:,2) = seism(:,1)
        seism(:,1) = 0.d0
        seism(:,3) = 0.d0
     endif

     close(10)


     istart = max(floor(tstart(irec)/deltat),1)
     iend = min(floor(tend(irec)/deltat),NSTEP)
     print *,'istart =',istart, 'iend =', iend
     print *,'tstart =',istart*deltat, 'tend =', iend*deltat
     if (istart >= iend) stop 'check istart,iend'
     nlenval = iend - istart +1

     do icomp = 1, NDIM

        print *,comp(icomp)

        filename = 'SEM/'//'AA.'//trim(station_name(irec))//'.'// comp(icomp) // '.adj'
        open(unit = 11, file = trim(filename), status='unknown', iostat=ier)
        if (ier /= 0) then
          print *,'Error: could not open adjoint trace file ',trim(filename)
          print *,'Please check if directory SEM/ exists...'
          stop 'Error opening trace file'
        endif


        time_window(:) = 0.d0
        seism_win(:) = seism(:,icomp)
        seism_veloc(:) = 0.d0
        seism_accel(:) = 0.d0

        do itime =istart,iend
           !time_window(itime) = 1.d0 - cos(pi*(itime-1)/NSTEP+1)**10   ! cosine window
           time_window(itime) = 1.d0 - (2* (dble(itime) - istart)/(iend-istart) -1.d0)**2  ! Welch window
        enddo

        do itime = 2,NSTEP-1
           seism_veloc(itime) = (seism_win(itime+1) - seism_win(itime-1))/(2*deltat)
        enddo
        seism_veloc(1) = (seism_win(2) - seism_win(1))/deltat
        seism_veloc(NSTEP) = (seism_win(NSTEP) - seism_win(NSTEP-1))/deltat

        do itime = 2,NSTEP-1
           seism_accel(itime) = (seism_veloc(itime+1) - seism_veloc(itime-1))/(2*deltat)
        enddo
        seism_accel(1) = (seism_veloc(2) - seism_veloc(1))/deltat
        seism_accel(NSTEP) = (seism_veloc(NSTEP) - seism_veloc(NSTEP-1))/deltat

        ! cross-correlation traveltime adjoint source
        Nnorm = deltat * sum(time_window(:) * seism_win(:) * seism_accel(:))
        !Nnorm = deltat * sum(time_window(:) * seism_veloc(:) * seism_veloc(:))
        if (abs(Nnorm) > EPS) then
           !ft_bar(:) = -seism_veloc(:) * time_window(:) / Nnorm
           ft_bar(:) = seism_veloc(:) * time_window(:) / Nnorm
           print *,'Norm =', Nnorm
        else
           print *, 'norm < EPS for file '
           print *,'Norm =', Nnorm
           ft_bar(:) = 0.d0
        endif

        do itime = 1,NSTEP
           if (icomp == adj_comp) then
              write(11,*) (itime-1)*deltat - t0, ft_bar(itime)
           else
              write(11,*) (itime-1)*deltat - t0, 0.d0
           endif
        enddo

     enddo
     close(11)

  enddo
  print *,'*************************'
  print *,'The input files (AA.S****.BXX/BXY/BXZ.adj) needed to run the adjoint simulation are in SEM'
  print *,'*************************'

end program adj_seismogram
