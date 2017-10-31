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

  subroutine check_stability()

! checks simulation stability and outputs timerun infos

#ifdef USE_MPI
  use mpi
#endif

  use constants, only: IMAIN,STABILITY_THRESHOLD,CUSTOM_REAL

  use specfem_par, only: myrank,timeval,it,NSTEP,GPU_MODE, &
                         ELASTIC_SIMULATION,any_elastic,displ_elastic, &
                         POROELASTIC_SIMULATION,any_poroelastic, &
                         displs_poroelastic,displw_poroelastic, &
                         ACOUSTIC_SIMULATION,any_acoustic,potential_acoustic, &
                         timestamp_seconds_start

  use specfem_par_noise, only: NOISE_TOMOGRAPHY

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: displnorm_all,displnorm_all_glob

  ! timer to count elapsed time
  double precision :: tCPU,t_remain,t_total,timestamp_seconds_current
  integer :: ihours,iminutes,iseconds,int_tCPU, &
             ihours_remain,iminutes_remain,iseconds_remain,int_t_remain, &
             ihours_total,iminutes_total,iseconds_total,int_t_total

  ! to determine date and time at which the run will finish
  character(len=8) :: datein
  character(len=10) :: timein
  character(len=5) :: zone
  integer, dimension(8) :: time_values
  character(len=3), dimension(12) :: month_name
  character(len=3), dimension(0:6) :: weekday_name
  data month_name /'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/
  data weekday_name /'Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat'/
  integer :: year,mon,day,hr,minutes,timestamp,julian_day_number,day_of_week
  integer, external :: idaywk

  ! checks if anything to do
  if (GPU_MODE) then
    ! todo: wavefields are on GPU and not transferred onto CPU yet for the following norm checks
    return
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '******************************************************************'
    if (timeval >= 1.d-3 .and. timeval < 1000.d0) then
      write(IMAIN,"('Time step number ',i7,'   t = ',f9.4,' s out of ',i7)") it,timeval,NSTEP
    else
      write(IMAIN,"('Time step number ',i7,'   t = ',1pe13.6,' s out of ',i7)") it,timeval,NSTEP
    endif
    write(IMAIN,*) '******************************************************************'
    write(IMAIN,*) 'We have done ',sngl(100.d0*dble(it-1)/dble(NSTEP-1)),'% of the total'
  endif


  ! elastic wavefield
  if (ELASTIC_SIMULATION) then
    if (any_elastic) then
      displnorm_all = maxval(sqrt(displ_elastic(1,:)**2 + displ_elastic(2,:)**2))
    else
      displnorm_all = 0._CUSTOM_REAL
    endif

    if (NOISE_TOMOGRAPHY /= 0) then
      if (myrank == 0) write(*,*) 'Noise simulation ', NOISE_TOMOGRAPHY, ' of 3'
    endif

    ! master collects norm from all processes
    call max_all_cr(displnorm_all, displnorm_all_glob)
    if (myrank == 0) &
      write(IMAIN,*) 'Max norm of vector field in solid (elastic) = ', displnorm_all_glob

    ! check stability of the code in solid, exit if unstable
    ! negative values can occur with some compilers when the unstable value is greater
    ! than the greatest possible floating-point number of the machine
    ! (is-not-a-number check is more robust when done on actual array values rather than return values from sqrt and amxval)
    if (displnorm_all > STABILITY_THRESHOLD .or. displnorm_all < 0._CUSTOM_REAL .or. &
! this trick checks for NaN (Not a Number), which is not even equal to itself
        displ_elastic(1,1) /= displ_elastic(1,1)) &
      call exit_MPI(myrank,'code became unstable and blew up in solid (elastic)')
  endif

  ! poroelastic wavefield
  if (POROELASTIC_SIMULATION) then
    if (any_poroelastic) then
      displnorm_all = maxval(sqrt(displs_poroelastic(1,:)**2 + displs_poroelastic(2,:)**2))
    else
      displnorm_all = 0._CUSTOM_REAL
    endif

    ! master collects norm from all processes
    call max_all_cr(displnorm_all, displnorm_all_glob)
    if (myrank == 0) &
      write(IMAIN,*) 'Max norm of vector field in solid (poroelastic) = ',displnorm_all_glob

    ! check stability of the code in solid, exit if unstable
    ! negative values can occur with some compilers when the unstable value is greater
    ! than the greatest possible floating-point number of the machine
    if (displnorm_all > STABILITY_THRESHOLD .or. displnorm_all < 0._CUSTOM_REAL .or. &
! this trick checks for NaN (Not a Number), which is not even equal to itself
        displs_poroelastic(1,1) /= displs_poroelastic(1,1)) &
      call exit_MPI(myrank,'code became unstable and blew up in solid (poroelastic)')

    if (any_poroelastic) then
      displnorm_all = maxval(sqrt(displw_poroelastic(1,:)**2 + displw_poroelastic(2,:)**2))
    else
      displnorm_all = 0._CUSTOM_REAL
    endif

    ! master collects norm from all processes
    call max_all_cr(displnorm_all, displnorm_all_glob)
    if (myrank == 0) &
      write(IMAIN,*) 'Max norm of vector field in fluid (poroelastic) = ',displnorm_all_glob

    ! check stability of the code in solid, exit if unstable
    ! negative values can occur with some compilers when the unstable value is greater
    ! than the greatest possible floating-point number of the machine
    if (displnorm_all > STABILITY_THRESHOLD .or. displnorm_all < 0._CUSTOM_REAL .or. &
! this trick checks for NaN (Not a Number), which is not even equal to itself
        displw_poroelastic(1,1) /= displw_poroelastic(1,1)) &
      call exit_MPI(myrank,'code became unstable and blew up in fluid (poroelastic)')
  endif


  ! acoustic wavefield
  if (ACOUSTIC_SIMULATION) then
    if (any_acoustic) then
      displnorm_all = maxval(abs(potential_acoustic(:)))
    else
      displnorm_all = 0._CUSTOM_REAL
    endif

    ! master collects norm from all processes
    call max_all_cr(displnorm_all, displnorm_all_glob)
    if (myrank == 0) &
      write(IMAIN,*) 'Max absolute value of scalar field in fluid (acoustic) = ',displnorm_all_glob

    ! check stability of the code in fluid, exit if unstable
    ! negative values can occur with some compilers when the unstable value is greater
    ! than the greatest possible floating-point number of the machine
    if (displnorm_all > STABILITY_THRESHOLD .or. displnorm_all < 0._CUSTOM_REAL .or. &
! this trick checks for NaN (Not a Number), which is not even equal to itself
        potential_acoustic(1) /= potential_acoustic(1)) &
      call exit_MPI(myrank,'code became unstable and blew up in fluid (acoustic)')
  endif

  ! count elapsed wall-clock time
  call date_and_time(datein,timein,zone,time_values)

  ! time_values(1): year
  ! time_values(2): month of the year
  ! time_values(3): day of the month
  ! time_values(5): hour of the day
  ! time_values(6): minutes of the hour
  ! time_values(7): seconds of the minute
  ! time_values(8): milliseconds of the second

  ! get timestamp in minutes of current date and time
  year = time_values(1)
  mon = time_values(2)
  day = time_values(3)
  hr = time_values(5)
  minutes = time_values(6)
  call convtime(timestamp,year,mon,day,hr,minutes)

  ! convert to seconds instead of minutes, to be more precise for 2D runs, which can be fast
  timestamp_seconds_current = timestamp*60.d0 + time_values(7) + time_values(8)/1000.d0

  ! elapsed time since beginning of the simulation
  if (myrank == 0) then
    tCPU = timestamp_seconds_current - timestamp_seconds_start
    int_tCPU = int(tCPU)
    ihours = int_tCPU / 3600
    iminutes = (int_tCPU - 3600*ihours) / 60
    iseconds = int_tCPU - 3600*ihours - 60*iminutes
    write(IMAIN,*) 'Elapsed time in seconds = ',tCPU
    write(IMAIN,"(' Elapsed time in hh:mm:ss = ',i6,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    write(IMAIN,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it)

    ! compute estimated remaining simulation time
    t_remain = (NSTEP - it) * (tCPU/dble(it))
    int_t_remain = int(t_remain)
    ihours_remain = int_t_remain / 3600
    iminutes_remain = (int_t_remain - 3600*ihours_remain) / 60
    iseconds_remain = int_t_remain - 3600*ihours_remain - 60*iminutes_remain
    write(IMAIN,*) 'Time steps remaining = ',NSTEP - it
    write(IMAIN,*) 'Estimated remaining time in seconds = ',t_remain
    write(IMAIN,"(' Estimated remaining time in hh:mm:ss = ',i6,' h ',i2.2,' m ',i2.2,' s')") &
           ihours_remain,iminutes_remain,iseconds_remain

    ! compute estimated total simulation time
    t_total = t_remain + tCPU
    int_t_total = int(t_total)
    ihours_total = int_t_total / 3600
    iminutes_total = (int_t_total - 3600*ihours_total) / 60
    iseconds_total = int_t_total - 3600*ihours_total - 60*iminutes_total
    write(IMAIN,*) 'Estimated total run time in seconds = ',t_total
    write(IMAIN,"(' Estimated total run time in hh:mm:ss = ',i6,' h ',i2.2,' m ',i2.2,' s')") &
           ihours_total,iminutes_total,iseconds_total

    if (it < NSTEP) then
      ! compute date and time at which the run should finish (useful for long runs)
      ! add remaining minutes and get date and time of that future timestamp in minutes
      timestamp = nint((timestamp_seconds_current + t_remain) / 60.d0)
      call invtime(timestamp,year,mon,day,hr,minutes)

      ! convert to Julian day to get day of the week
      call calndr(day,mon,year,julian_day_number)
      day_of_week = idaywk(julian_day_number)

      write(IMAIN,"(' The run will finish approximately on: ',a3,' ',a3,' ',i2.2,', ',i4.4,' ',i2.2,':',i2.2)") &
          weekday_name(day_of_week),month_name(mon),day,year,hr,minutes

    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine check_stability

