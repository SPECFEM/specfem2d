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
!               Pieyre Le Loher, pieyre DOT le-loher aT inria.fr
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


  subroutine set_sources()

! gets source parameters

  use specfem_par, only : myrank,NSOURCES,source_type,time_function_type, &
                          x_source,z_source,Mxx,Mzz,Mxz,f0,tshift_src,factor,anglesource,aval, &
                          t0,initialfield,deltat,USER_T0

  implicit none
  include "constants.h"

  ! local parameters
  integer :: i_source
  double precision, dimension(NSOURCES) :: t0_source,hdur
  double precision :: min_tshift_src_original

  ! checks the input
  do i_source= 1,NSOURCES

    ! checks source type
    if (.not. initialfield) then
      if (source_type(i_source) == 1) then
        if (myrank == 0) then
          ! user output
          write(IMAIN,212) x_source(i_source),z_source(i_source),f0(i_source),tshift_src(i_source), &
                       factor(i_source),anglesource(i_source)
        endif
      else if (source_type(i_source) == 2) then
        if (myrank == 0) then
          ! user output
          write(IMAIN,222) x_source(i_source),z_source(i_source),f0(i_source),tshift_src(i_source), &
                       factor(i_source),Mxx(i_source),Mzz(i_source),Mxz(i_source)
        endif
      else
        call exit_MPI(myrank,'Unknown source type number !')
      endif
    endif

    ! if Dirac source time function, use a very thin Gaussian instead
    ! if Heaviside source time function, use a very thin error function instead
    if (time_function_type(i_source) == 4 .or. time_function_type(i_source) == 5) &
      f0(i_source) = 1.d0 / (10.d0 * deltat)

    ! checks source frequency
    if (abs(f0(i_source)) < TINYVAL) then
      call exit_MPI(myrank,'Error source frequency is zero')
    endif

    ! half-duration of source
    hdur(i_source) = 1.d0 / f0(i_source)

    ! sets source start times, shifted by the given (non-zero) time-shift
    if (time_function_type(i_source)== 5) then
      t0_source(i_source) = 2.0d0 * hdur(i_source) + tshift_src(i_source)
    else
      t0_source(i_source) = 1.20d0 * hdur(i_source) + tshift_src(i_source)
    endif

    ! for the source time function
    aval(i_source) = PI*PI*f0(i_source)*f0(i_source)

    ! convert angle from degrees to radians
    anglesource(i_source) = anglesource(i_source) * PI / 180.d0

  enddo ! do i_source= 1,NSOURCES

  ! initializes simulation start time
  if (NSOURCES == 1) then
    ! simulation start time
    t0 = t0_source(1)
    ! sets source time shift relative to simulation start time
    min_tshift_src_original = tshift_src(1)
    tshift_src(1) = 0.d0
  else
    ! starts with earliest start time
    t0 = minval( t0_source(:) )
    ! sets source time shifts relative to simulation start time
    min_tshift_src_original = minval( tshift_src(:) )
    tshift_src(:) = t0_source(:) - t0
  endif

  ! checks if user set USER_T0 to fix simulation start time
  ! note: USER_T0 has to be positive
  if (USER_T0 > 0.d0) then
    ! user cares about origin time and time shifts of the CMTSOLUTION
    ! and wants to fix simulation start time to a constant start time
    ! time 0 on time axis will correspond to given origin time

    ! notifies user
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '    using USER_T0 . . . . . . . . . = ',USER_T0
      write(IMAIN,*) '      original t0 . . . . . . . . . = ',t0
      write(IMAIN,*) '      min_tshift_src_original . . . = ',min_tshift_src_original
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! checks if automatically set t0 is too small
    ! note: times in seismograms are shifted by t0(1)
    if (t0 <= USER_T0 + min_tshift_src_original) then

      ! sets new simulation start time such that
      ! simulation starts at t = - t0 = - USER_T0
      t0 = USER_T0

      ! notifies user
      if (myrank == 0) then
        write(IMAIN,*) '    fix new simulation start time . = ', - t0
      endif

      ! loops over all sources
      do i_source= 1,NSOURCES
        ! sets the given, initial time shifts
        if (time_function_type(i_source) == 5) then
          tshift_src(i_source) = t0_source(i_source) - 2.0d0 * hdur(i_source)
        else
          tshift_src(i_source) = t0_source(i_source) - 1.20d0 * hdur(i_source)
        endif
        ! user output
        if (myrank == 0) then
          write(IMAIN,*) '    source ',i_source,'uses tshift = ',tshift_src(i_source)
        endif
      enddo
      ! user output
      if (myrank == 0) then
        write(IMAIN,*)
      endif

    else
      ! start time needs to be at least t0 for numerical stability
      ! notifies user
      if (myrank == 0) then
        write(IMAIN,*) 'Error: USER_T0 is too small'
        write(IMAIN,*) '       must make one of three adjustements:'
        write(IMAIN,*) '       - increase USER_T0 to be at least: ',t0
        write(IMAIN,*) '       - decrease time shift tshift_src in SOURCE file'
        write(IMAIN,*) '       - increase frequency f0 in SOURCE file'
      endif
      call exit_MPI(myrank,'Error USER_T0 is set but too small')
    endif
  else if (USER_T0 < 0.d0) then
    if (myrank == 0) then
      write(IMAIN,*) 'error: USER_T0 is negative, must be set zero or positive!'
    endif
    call exit_MPI(myrank,'Error negative USER_T0 parameter in constants.h')
  endif

  ! checks onset times
  if (.not. initialfield) then

    ! loops over sources
    do i_source = 1,NSOURCES

      ! excludes Dirac and Heaviside sources
      if (time_function_type(i_source) /= 4 .and. time_function_type(i_source) /= 5) then

        ! user output
        if (myrank == 0) then
          write(IMAIN,*) '    Onset time. . . . . . = ',- (t0+tshift_src(i_source))
          write(IMAIN,*) '    Fundamental period. . = ',1.d0/f0(i_source)
          write(IMAIN,*) '    Fundamental frequency = ',f0(i_source)
        endif

        ! checks source onset time
        if (t0+tshift_src(i_source) < 1.d0/f0(i_source)) then
          call exit_MPI(myrank,'Onset time too small')
        else
          if (myrank == 0) then
            write(IMAIN,*) '    --> onset time ok'
          endif
        endif
      endif
    enddo

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
