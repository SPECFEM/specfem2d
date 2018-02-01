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

program adj_seismogram

! This program cuts a certain portion of the seismograms and convert it
! into the adjoint source for generating banana-dougnut kernels
!
! The calculated adjoint source corresponds to a (cross-correlation) traveltime measurement
! within a given time window
!
! compilation:
!    rm -r -f xadj_seismogram ; gfortran adj_seismogram.f90 -o xadj_seismogram
!
! usage example:
!    ./xadj_seismogram 27.0 32.  AA.S0001  1

  implicit none

  ! number of time steps
  integer :: NSTEP

  ! time step size
  double precision :: deltat

  integer :: itime,icomp,istart,iend,nlen,irec,i
  integer :: NDIM,NDIMr,adj_comp
  integer :: ier

  double precision :: t0,time,time0,dummy,Nnorm

  double precision, dimension(:),allocatable :: time_window
  double precision, dimension(:,:),allocatable :: seism
  double precision, dimension(:),allocatable :: seism_veloc,seism_accel,ft_bar,seism_win

  character(len=3) :: compr(2),comp(3)
  character(len=150) :: filename

  ! single station setup
  integer, parameter :: nrec = 1
  double precision,dimension(nrec) :: tstart(nrec),tend(nrec)
  character(len=150), dimension(nrec) :: station_name

  ! input
  character(len=256) :: arg

  ! tolerance
  double precision, parameter :: EPS = 1.d-40
  double precision, parameter :: SMALLVAL = 1.d-9
  double precision, parameter :: PI = 3.141592653589793d0

  ! window taper: 1 == Welch / 2 == cosine
  integer,parameter :: window_taper_type = 1


!--------------------------------------------------------------------------------
! USER PARAMETERS
!--------------------------------------------------------------------------------

  ! reads in file arguments
  do i = 1,4
    call get_command_argument(i,arg)
    if (trim(arg) == '' .or. command_argument_count() /= 4) then
      print *,'Usage: '
      print *,'  xadj_seismogram t1 t2 station-name adj_comp[1-3]'
      print *,'with'
      print *,'  t1: window start time'
      print *,'  t2: window end time'
      print *,'  station-name : adjoint station (e.g. AA.S00001)'
      print *,'  adj_comp: 1 = adjoint source given by X component'
      print *,'  adj_comp: 2 = adjoint source given by Y component (SH waves)'
      print *,'  adj_comp: 3 = adjoint source given by Z component'
      stop 1
    endif
    select case (i)
    case (1)
      read(arg,*,iostat=ier) tstart(1)
      if (ier /= 0) stop 'Error reading tstart'
    case (2)
      read(arg,*,iostat=ier) tend(1)
      if (ier /= 0) stop 'Error reading tend'
    case (3)
      station_name(1) = trim(arg)
    case (4)
      read(arg,*,iostat=ier) adj_comp
      if (ier /= 0) stop 'Error reading adjoint component'
    end select
  enddo

  ! checks input arguments
  if (adj_comp /= 1 .and. adj_comp /= 2 .and. adj_comp /= 3) then
    print *,'Invalid adj_comp number, must be 1 (X-comp), 2 (Y-comp) or 3 (Z-comp)'
    stop 1
  endif

  !---------------------------------
  ! trace component setup
  !---------------------------------
  ! number of components
  !NDIMr = 2  ! P-SV
  !NDIMr = 1  ! SH (membrane)

  !---------------------------------
  ! single station setup
  !---------------------------------
  ! list of stations
  !station_name(1) = 'AA.'//'S0001'

  ! KEY: 'absolute' time interval for window
  !tstart(1) = 27.0 + 8.0
  !tend(1) = 32.0 + 8.0

  ! chose the component for the adjoint source (adj_comp = 1:X, 2:Y, 3:Z)
  !adj_comp = 1
  !---------------------------------

  ! all components
  NDIM = 3
  comp = (/"BXX","BXY","BXZ"/)

  ! select number of dimensions of input trace
  if (adj_comp == 2) then
    ! Y-component only for SH-waves
    NDIMr = 1
  else
    ! X- or Z-component for P-SV waves
    NDIMr = 2
  endif

  ! select wave type components
  if (NDIMr == 1) then
    ! SH (membrane)
    compr = (/"BXY","tmp"/)
  else
    ! P-SV
    compr = (/"BXX","BXZ"/)
  endif

  ! user output
  print *,'adjoint source - seismogram:'
  print *
  print *,'setup:'
  print *,'  number of adjoint stations (nrec)   = ',nrec
  do irec = 1,nrec
    print *,'  station name ',irec,':    ',station_name(irec)
  enddo
  print *
  print *,'  seismogram components   = ',NDIMr
  if (NDIMr == 1) then
    print *,'  seismogram label        = ',compr(1)
  else
    print *,'  seismogram labels       = ',compr(1),' / ',compr(2)
  endif
  print *
  print *,'  time window start/end                           = ',tstart(1),tend(1)
  print *,'  adjoint source trace component (1 == X / 2 == Y / 3 == Z) = ',adj_comp
  print *

  ! basic check
  do irec = 1,nrec
    if (tstart(irec) >= tend(irec)) then
      print *,"Invalid start/end time chosen for station ",irec,": tstart = ",tstart(irec), &
              " should be less then tend = ",tend(irec)
      stop 1
    endif
  enddo

  ! initializes
  NSTEP = 0
  t0 = 0.d0

  ! gets trace length and start time t0
  do irec = 1,nrec

    do icomp = 1,NDIMr

      filename = 'OUTPUT_FILES/'//trim(station_name(irec))//'.'// compr(icomp) // '.semd'
      open(unit = 10, file = trim(filename),status='old',iostat=ier)
      if (ier /= 0) then
        print *,'Error opening file: ',trim(filename)
        !stop 'Error opening trace file' ! note: stop '**comment**' finishes with an error code 0 on some systems
        stop 10
      endif

      ! counts lines
      nlen = 0
      do while (ier == 0)
        read(10,*,iostat=ier) time , dummy
        if (ier == 0) nlen = nlen + 1
      enddo
      rewind(10)

      ! sets length NSTEP from first file
      if (irec == 1 .and. icomp == 1) then
        NSTEP = nlen
      else
        ! checks with settings
        if (nlen /= NSTEP) then
          print *,'Error: number of lines ',nlen,' in file ',trim(filename),' do not match setting NSTEP = ',NSTEP
          stop 20
        endif
      endif

      ! check
      if (NSTEP <= 1) then
        print *,'Error: invalid number of lines (time steps) encountered: NSTEP = ',NSTEP
        stop 30
      endif

      ! sets time step size
      do itime = 1,NSTEP
        read(10,*) time , dummy

        ! stores start time and time step size
        if (irec == 1 .and. icomp == 1) then
          ! sets start time
          if (itime == 1) t0 = time
          ! sets time step
          if (itime == NSTEP) deltat = (time - t0)/dble(NSTEP-1)
        else
          ! checks start time
          if (itime == 1) then
            time0 = time
            if (abs(t0 - time0) > SMALLVAL) then
              print *,'Error: start time',time0,' in file ',trim(filename),' does not match t0 = ',t0
              stop 40
            endif
          endif

          ! checks time step size
          if (itime == NSTEP) then
            if (abs( (time - time0)/dble(NSTEP-1) - deltat) > SMALLVAL) then
              print *,'Error: time step size ',(time-time0)/dble(NSTEP-1),' in file ',trim(filename), &
                      ' does not match deltat = ',deltat
              stop 50
            endif
          endif
        endif
      enddo
      close(10)

    enddo
  enddo

  ! user output
  print *,'reading input traces:'
  print *,'  number of time steps (NSTEP)  = ',NSTEP
  print *
  print *,'  start time (t0)               = ',t0
  print *,'  end time                      = ',(NSTEP-1)*deltat + t0
  print *,'  time step size (deltat)       = ',deltat
  print *

  ! allocates trace arrays
  allocate( time_window(NSTEP), &
            seism(NSTEP,3), &
            seism_win(NSTEP), &
            seism_veloc(NSTEP), &
            seism_accel(NSTEP), &
            ft_bar(NSTEP),stat=ier)
  if (ier /= 0) stop 2

  ! creates adjoint sources
  do irec = 1,nrec

    ! reads in all trace components
    do icomp = 1,NDIMr

      filename = 'OUTPUT_FILES/'//trim(station_name(irec))//'.'// compr(icomp) // '.semd'
      open(unit = 10, file = trim(filename),status='old',iostat=ier)
      if (ier /= 0) then
        print *,'Error opening file: ',trim(filename)
        stop 'Error opening trace file'
      endif

      ! reads in trace
      do itime = 1,NSTEP
        read(10,*) time , seism(itime,icomp)
      enddo
      close(10)
    enddo

    if (NDIMr == 2) then
      seism(:,3) = seism(:,2)
      seism(:,2) = 0.d0
    else
      seism(:,2) = seism(:,1)
      seism(:,1) = 0.d0
      seism(:,3) = 0.d0
    endif

    ! start/end index
    ! (note that early start times have negative t0. it needs to be added to find correct index)
    istart = floor((tstart(irec) - t0)/deltat) + 1
    iend = ceiling((tend(irec) - t0)/deltat) + 1

    ! user output
    if (istart < 1 .or. istart > NSTEP -1) &
      print *,'*** warning: time window start is out of trace length! start time will be moved ***'
    if (iend < 1 .or. iend > NSTEP) &
      print *,'*** warning: time window end is out of trace length! end time will be moved ***'

    ! limits
    if (istart < 1) istart = 1
    if (istart > NSTEP) istart = NSTEP - 1
    if (iend < istart) iend = istart + 1
    if (iend > NSTEP) iend = NSTEP

    ! user output
    print *,'time window:'
    print *,'  index      : istart =',istart, 'iend =', iend
    print *,'  time window: tstart =',(istart-1)*deltat + t0, 'tend =',(iend-1)*deltat + t0
    print *

    if (istart >= iend) then
      print *,"Error start/end index: ",istart,iend
      stop 11
    endif

    ! window length
    nlen = iend - istart + 1

    ! outputs all 3-component X/Y/Z (needed for kernel simulations)
    do icomp = 1, NDIM

      print *,'component: ',comp(icomp)

      filename = 'SEM/'//trim(station_name(irec))//'.'// comp(icomp) // '.adj'
      open(unit = 11, file = trim(filename),status='unknown',iostat=ier)
      if (ier /= 0) then
        print *,'Error opening file: ',trim(filename)
        stop 60
      endif

      ! time window
      time_window(:) = 0.d0
      select case (window_taper_type)
      case (1)
        ! Welch window
        do itime = istart,iend
          time_window(itime) = 1.d0 - (2* (dble(itime) - istart)/(iend-istart) -1.d0)**2
        enddo

      case (2)
        ! cosine window
        do itime = istart,iend
          time_window(itime) = 1.d0 - cos(PI*(itime-1)/NSTEP+1)**10
        enddo

      case default
        print *,'Invalid window taper type ',window_taper_type
        print *,'Please choose 1 == Welch or 2 == cosine and re-compile before running'
        stop 60
      end select

      seism_win(:) = seism(:,icomp)

      ! first time derivative (by finite-differences)
      seism_veloc(:) = 0.d0
      do itime = 2,NSTEP-1
         seism_veloc(itime) = (seism_win(itime+1) - seism_win(itime-1))/(2*deltat)
      enddo
      seism_veloc(1) = (seism_win(2) - seism_win(1))/deltat
      seism_veloc(NSTEP) = (seism_win(NSTEP) - seism_win(NSTEP-1))/deltat

      ! second time derivative
      seism_accel(:) = 0.d0
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
         print *,'Norm < EPS for file, zeroeing out trace'
         print *,'Norm =', Nnorm
         ft_bar(:) = 0.d0
      endif
      print *

      if (icomp == adj_comp) then
        do itime = 1,NSTEP
          write(11,*) (itime-1)*deltat + t0, ft_bar(itime)
        enddo
      else
        do itime = 1,NSTEP
          write(11,*) (itime-1)*deltat + t0, 0.d0
        enddo
      endif
      close(11)

    enddo

  enddo

  ! user output
  print *,'*************************'
  print *,'The input files (AA.S****.BXX/BXY/BXZ.adj) needed to run the adjoint simulation are in folder SEM/'
  print *,'*************************'

end program adj_seismogram
