
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
!               Pieyre Le Loher, pieyre DOT le-loher aT inria.fr
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


  subroutine set_sources(myrank,NSOURCE,source_type,time_function_type, &
                      x_source,z_source,Mxx,Mzz,Mxz,f0,tshift_src,factor,angleforce,aval, &
                      t0,initialfield,ipass,deltat)

! gets source parameters

  implicit none
  include "constants.h"

  integer :: myrank
  integer :: NSOURCE
  integer, dimension(NSOURCE) :: source_type,time_function_type
  double precision, dimension(NSOURCE) :: x_source,z_source, &
    Mxx,Mzz,Mxz,f0,tshift_src,factor,angleforce
  double precision, dimension(NSOURCE) :: aval
  double precision :: t0
  double precision :: deltat
  integer :: ipass
  logical :: initialfield

  ! local parameters
  integer :: i_source

  ! checks the input
  do i_source=1,NSOURCE

    ! checks source type
    if(.not. initialfield) then
      if (source_type(i_source) == 1) then
        if ( myrank == 0 .and. ipass == 1 ) then

          write(IOUT,212) x_source(i_source),z_source(i_source),f0(i_source),tshift_src(i_source), &
                       factor(i_source),angleforce(i_source)

        endif
      else if(source_type(i_source) == 2) then
        if ( myrank == 0 .and. ipass == 1 ) then

          write(IOUT,222) x_source(i_source),z_source(i_source),f0(i_source),tshift_src(i_source), &
                       factor(i_source),Mxx(i_source),Mzz(i_source),Mxz(i_source)

        endif
      else
        call exit_MPI('Unknown source type number !')
      endif
    endif

    ! note: this is slightly different than in meshfem2D.f90,
    !          t0 will only be set within this if statement, i.e. only for type 4 or 5 sources
    !          (since we set f0 to a new values for these two types of sources)

    ! if Dirac source time function, use a very thin Gaussian instead
    ! if Heaviside source time function, use a very thin error function instead
    if(time_function_type(i_source) == 4 .or. time_function_type(i_source) == 5) then
      f0(i_source) = 1.d0 / (10.d0 * deltat)

      ! time delay of the source in seconds, use a 20 % security margin (use 2 / f0 if error function)
      if(time_function_type(i_source)== 5) then
        tshift_src(i_source) = 2.0d0 / f0(i_source) + tshift_src(i_source)
      else
        tshift_src(i_source) = 1.20d0 / f0(i_source) + tshift_src(i_source)
      endif
    endif

    ! for the source time function
    aval(i_source) = PI*PI*f0(i_source)*f0(i_source)

    ! convert angle from degrees to radians
    angleforce(i_source) = angleforce(i_source) * PI / 180.d0

  enddo ! do i_source=1,NSOURCE

  ! initializes simulation start time
  t0 = tshift_src(1)

  ! checks if user set USER_T0 to fix simulation start time
  ! note: USER_T0 has to be positive
  if( USER_T0 > 0.d0 ) then
    ! user cares about origin time and time shifts of the CMTSOLUTION
    ! and wants to fix simulation start time to a constant start time
    ! time 0 on time axis will correspond to given origin time

    ! notifies user
    if( myrank == 0 .and. ipass == 1) then
      write(IOUT,*)
      write(IOUT,*) '    USER_T0: ',USER_T0,'initial t0_start: ',t0
    endif

    ! checks if automatically set t0 is too small
    ! note: times in seismograms are shifted by t0(1)
    if( t0 <= USER_T0 ) then
      ! sets new simulation start time such that
      ! simulation starts at t = - t0 = - USER_T0
      t0 = USER_T0

      ! notifies user
      if( myrank == 0 .and. ipass == 1) then
        write(IOUT,*) '    fix new simulation start time. . . . . = ', - t0
        write(IOUT,*)
      endif

      ! loops over all sources
      do i_source=1,NSOURCE
        ! gets the given, initial time shifts
        if( time_function_type(i_source) == 5 ) then
          tshift_src(i_source) = tshift_src(i_source) - 2.0d0 / f0(i_source)
        else
          tshift_src(i_source) = tshift_src(i_source) - 1.20d0 / f0(i_source)
        endif

        ! sets new t0 according to simulation start time,
        ! using absolute time shifts for each source such that
        ! a zero time shift would have the maximum gaussian at time t = (it-1)*DT - t0_start = 0
        tshift_src(i_source) = USER_T0 + tshift_src(i_source)

        if( myrank == 0 .and. ipass == 1) then
          write(IOUT,*) '    source ',i_source,'uses tshift . . . . . = ',tshift_src(i_source)
        endif

      enddo

    else
      ! start time needs to be at least t0 for numerical stability
      ! notifies user
      if( myrank == 0 .and. ipass == 1) then
        write(IOUT,*) 'error: USER_T0 is too small'
        write(IOUT,*) '       must make one of three adjustements:'
        write(IOUT,*) '       - increase USER_T0 to be at least: ',t0
        write(IOUT,*) '       - decrease time shift tshift_src in SOURCE file'
        write(IOUT,*) '       - increase frequency f0 in SOURCE file'
      endif
      call exit_MPI('error USER_T0 is set but too small')
    endif
  else if( USER_T0 < 0.d0 ) then
    if( myrank == 0 .and. ipass == 1 ) then
      write(IOUT,*) 'error: USER_T0 is negative, must be set zero or positive!'
    endif
    call exit_MPI('error negative USER_T0 parameter in constants.h')
  endif

  ! output formats
212 format(//,5x,'Source Type. . . . . . . . . . . . . . = Collocated Force',/5x, &
                  'X-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
                  'Y-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
                  'Fundamental frequency (Hz) . . . . . . =',1pe20.10,/5x, &
                  'Time delay (s) . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Multiplying factor . . . . . . . . . . =',1pe20.10,/5x, &
                  'Angle from vertical direction (deg). . =',1pe20.10,/5x)

222 format(//,5x,'Source Type. . . . . . . . . . . . . . = Moment-tensor',/5x, &
                  'X-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
                  'Y-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
                  'Fundamental frequency (Hz) . . . . . . =',1pe20.10,/5x, &
                  'Time delay (s) . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Multiplying factor . . . . . . . . . . =',1pe20.10,/5x, &
                  'Mxx. . . . . . . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Mzz. . . . . . . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Mxz. . . . . . . . . . . . . . . . . . =',1pe20.10)

  end subroutine set_sources
