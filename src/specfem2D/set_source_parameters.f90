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

  subroutine set_source_parameters()

! gets source parameters and determines simulation start/onset time

  use constants, only: NGLLX,NGLLZ,NDIM,IMAIN,TINYVAL,PI,F_ORMSBY

  use specfem_par, only: myrank,NSOURCES,source_type,time_function_type, &
                         x_source,z_source,Mxx,Mzz,Mxz,f0_source,tshift_src,factor,anglesource, &
                         t0,initialfield,USER_T0,vx_source,vz_source,SOURCE_IS_MOVING

  implicit none

  ! local parameters
  integer :: i_source
  double precision, dimension(NSOURCES) :: t0_source
  double precision :: min_tshift_src_original,hdur
  double precision :: M0,Mw

  ! note: see also routine in read_source_file.f90 which reads in source parameters

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  Total number of sources: ', NSOURCES
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! checks the input
  do i_source = 1,NSOURCES

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  Setting parameters for source',i_source
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! half-duration of the marmousi_ormsby_wavelet whose freq signature is 5-10-60-80
    ! so here we choose 35Hz as dominant frequency
    if (time_function_type(i_source) == 11) then
      ! warning
      if (f0_source(i_source) /= F_ORMSBY) then
        print *,'warning: Marmousi Ormsby wavelet dominant frequency will be set to 35 Hz'
        f0_source(i_source) = F_ORMSBY
      endif
    endif

    ! user output for source type
    if (.not. initialfield) then
      if (source_type(i_source) == 1) then
        ! force
        if (myrank == 0) then
          ! user output
          if (SOURCE_IS_MOVING) then
            write(IMAIN,213) x_source(i_source),z_source(i_source),vx_source(i_source),vz_source(i_source), &
                             f0_source(i_source), tshift_src(i_source), factor(i_source),anglesource(i_source)
          else
            write(IMAIN,212) x_source(i_source),z_source(i_source),f0_source(i_source),tshift_src(i_source), &
                             factor(i_source),anglesource(i_source)
          endif
          write(IMAIN,*)
        endif
      else if (source_type(i_source) == 2) then
        ! moment tensor
        if (myrank == 0) then
          ! user output
          if (SOURCE_IS_MOVING) then
            write(IMAIN,223) x_source(i_source),z_source(i_source),vx_source(i_source),vz_source(i_source), &
                             f0_source(i_source),tshift_src(i_source), factor(i_source),Mxx(i_source), &
                             Mzz(i_source),Mxz(i_source)
          else
            write(IMAIN,222) x_source(i_source),z_source(i_source),f0_source(i_source),tshift_src(i_source), &
                             factor(i_source),Mxx(i_source),Mzz(i_source),Mxz(i_source)
          endif
          write(IMAIN,*)
        endif
      else
        call exit_MPI(myrank,'Unknown source type number !')
      endif
    endif

    ! half-duration of source
    hdur = 1.d0 / f0_source(i_source)

    ! sets source start times, shifted by the given (non-zero) time-shift
    if (time_function_type(i_source) == 5 .or. time_function_type(i_source) == 11) then
      ! Heaviside or Ormsby
      t0_source(i_source) = 2.0d0 * hdur + tshift_src(i_source)
    else if (time_function_type(i_source) == 10) then
      t0_source(i_source) = 0.0
    else
      t0_source(i_source) = 1.20d0 * hdur + tshift_src(i_source)
    endif

    ! moment-tensor source
    ! note: (mentioned also in SEM2DPACK by J.-P. Ampuero)
    !       Moment tensor source in 2D is equivalent to a line source in 3D.
    !       The scaling from 3D to 2D should thus be corrected to have similar amplitudes between 2D and 3D simulations.
    !       This scaling factor is related to the characteristic size of the rupture area in the off-plane dimension.
    !       According to Pablo Ampuero, an approximate scaling for a fault area with aspect ratio close to 1 is
    !         M0_2D = dtau * (M0_3D/dtau)**2/3
    !       with dtau being the stress drop (typically a few MPa) of the rupture event.
    if (.not. initialfield) then
      if (source_type(i_source) == 2) then
        ! scalar moment
        ! simplified for 2D: Myy == 0, Mxy == 0, Myz == 0
        ! euclidean (or Frobenius) norm of a matrix: M0**2 = sum( Mij**2 )
        !  M0 = Mxx**2 + Myy**2 + Mzz**2 + 2.d0 * ( Mxy**2 + Mxz**2 + Myz**2 )
        ! for 2D:
        M0 = Mxx(i_source)**2 + Mzz(i_source)**2 + 2.d0 * (Mxz(i_source)**2)
        ! adds 1/2 to be coherent with double couple or point sources
        M0 = dsqrt(M0/2.0d0)
        ! note: CMTSOLUTION file values are in dyne.cm
        !       1 dyne is 1 gram * 1 cm / (1 second)^2
        !       1 Newton is 1 kg * 1 m / (1 second)^2
        ! for SPECFEM2D we assume M given in N-m: thus 1 Newton = 100,000 dynes
        !                                         and  1 N-m    = 1.e7 dyne-cm    (or  1 dyne.cm = 1e-7 Newton.m)
        ! scalar moment in dyne-cm
        M0 = M0 * 1.d7

        ! moment magnitude
        ! moment magnitude by Hanks & Kanamori, 1979
        ! Mw = 2/3 log( M0 ) - 10.7       (dyne-cm)
        if (M0 > 0.d0) then
          Mw = 2.d0/3.d0 * log10( M0 ) - 10.7d0
        else
          Mw = 0.d0
        endif
        if (myrank == 0) then
          write(IMAIN,*) '    2D seismic moment M0 = ',sngl(M0),'(dyne-cm)'
          write(IMAIN,*) '        moment magnitude = ',sngl(Mw)
          write(IMAIN,*)
          call flush_IMAIN()
        endif
      endif
    endif

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
      do i_source = 1,NSOURCES
        ! half-duration of source
        hdur = 1.d0 / f0_source(i_source)

        ! sets the given, initial time shifts
        if (time_function_type(i_source) == 5) then
          tshift_src(i_source) = t0_source(i_source) - 2.0d0 * hdur
        else
          tshift_src(i_source) = t0_source(i_source) - 1.20d0 * hdur
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
      write(IMAIN,*) 'Error: USER_T0 is negative, must be set zero or positive!'
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
          write(IMAIN,*) '  Source: ',i_source
          write(IMAIN,*) '    Onset time. . . . . . = ',- (t0+tshift_src(i_source))
          write(IMAIN,*) '    Fundamental period. . = ',1.d0/f0_source(i_source)
          write(IMAIN,*) '    Fundamental frequency = ',f0_source(i_source)
        endif

        ! checks source onset time
        if (t0+tshift_src(i_source) < 1.d0/f0_source(i_source) .and. time_function_type(i_source) /= 10) then
          call exit_MPI(myrank,'Onset time too small')
        else
          if (myrank == 0) then
            write(IMAIN,*) '    The onset time is ok'
            write(IMAIN,*)
            call flush_IMAIN()
          endif
        endif
      endif
    enddo

  endif


  ! output formats
212 format(5x,'Source Type. . . . . . . . . . . . . . = Collocated Force',/5x, &
                  'X-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
                  'Y-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
                  'Fundamental frequency (Hz) . . . . . . =',1pe20.10,/5x, &
                  'Time delay (s) . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Multiplying factor . . . . . . . . . . =',1pe20.10,/5x, &
                  'Angle from vertical direction (deg). . =',1pe20.10,/5x)

213 format(5x,'Source Type. . . . . . . . . . . . . . = Collocated Force',/5x, &
                  'Initial X-position (meters). . . . . . =',1pe20.10,/5x, &
                  'Initial Y-position (meters). . . . . . =',1pe20.10,/5x, &
                  'X-velocity (meters/second). . .  . . . =',1pe20.10,/5x, &
                  'Z-velocity (meters/second). . .  . . . =',1pe20.10,/5x, &
                  'Fundamental frequency (Hz) . . . . . . =',1pe20.10,/5x, &
                  'Time delay (s) . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Multiplying factor . . . . . . . . . . =',1pe20.10,/5x, &
                  'Angle from vertical direction (deg). . =',1pe20.10,/5x)

222 format(5x,'Source Type. . . . . . . . . . . . . . = Moment-tensor',/5x, &
                  'X-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
                  'Y-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
                  'Fundamental frequency (Hz) . . . . . . =',1pe20.10,/5x, &
                  'Time delay (s) . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Multiplying factor . . . . . . . . . . =',1pe20.10,/5x, &
                  'Mxx. . . . . . . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Mzz. . . . . . . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Mxz. . . . . . . . . . . . . . . . . . =',1pe20.10)

223 format(5x,'Source Type. . . . . . . . . . . . . . = Moment-tensor',/5x, &
                  'Initial X-position (meters). . . . . . =',1pe20.10,/5x, &
                  'Initial Y-position (meters). . . . . . =',1pe20.10,/5x, &
                  'X-velocity (meters/second). . .  . . . =',1pe20.10,/5x, &
                  'Z-velocity (meters/second). . .  . . . =',1pe20.10,/5x, &
                  'Fundamental frequency (Hz) . . . . . . =',1pe20.10,/5x, &
                  'Time delay (s) . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Multiplying factor . . . . . . . . . . =',1pe20.10,/5x, &
                  'Mxx. . . . . . . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Mzz. . . . . . . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Mxz. . . . . . . . . . . . . . . . . . =',1pe20.10)


  end subroutine set_source_parameters
