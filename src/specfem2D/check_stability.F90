
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

  subroutine check_stability()

! checks simulation stability and outputs timerun infos

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par, only: myrank,timeval,it,NSTEP,NOISE_TOMOGRAPHY, &
                         any_elastic_glob,any_elastic,displ_elastic, &
                         any_poroelastic_glob,any_poroelastic, &
                         displs_poroelastic,displw_poroelastic, &
                         any_acoustic_glob,any_acoustic,potential_acoustic, &
                         timestamp_seconds_start

  implicit none
  include "constants.h"

  ! local parameters
  double precision displnorm_all,displnorm_all_glob
  ! timer to count elapsed time
  double precision :: tCPU,t_remain,t_total,timestamp_seconds_current
  integer :: ihours,iminutes,iseconds,int_tCPU, &
             ihours_remain,iminutes_remain,iseconds_remain,int_t_remain, &
             ihours_total,iminutes_total,iseconds_total,int_t_total
  ! to determine date and time at which the run will finish
  character(len=8) datein
  character(len=10) timein
  character(len=5)  :: zone
  integer, dimension(8) :: time_values
  character(len=3), dimension(12) :: month_name
  character(len=3), dimension(0:6) :: weekday_name
  data month_name /'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/
  data weekday_name /'Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat'/
  integer :: year,mon,day,hr,minutes,timestamp,julian_day_number,day_of_week
  integer, external :: idaywk
#ifdef USE_MPI
  integer :: ier
#endif


  ! user output
  if (myrank == 0) then
    write(IOUT,*)
    write(IOUT,*) '******************************************************************'
    if(timeval >= 1.d-3 .and. timeval < 1000.d0) then
      write(IOUT,"('Time step number ',i7,'   t = ',f9.4,' s out of ',i7)") it,timeval,NSTEP
    else
      write(IOUT,"('Time step number ',i7,'   t = ',1pe13.6,' s out of ',i7)") it,timeval,NSTEP
    endif
    write(IOUT,*) '******************************************************************'
    write(IOUT,*) 'We have done ',sngl(100.d0*dble(it-1)/dble(NSTEP-1)),'% of the total'
  endif


  ! elastic wavefield
  if(any_elastic_glob) then
    if(any_elastic) then
      displnorm_all = maxval(sqrt(displ_elastic(1,:)**2 &
                                + displ_elastic(2,:)**2 &
                                + displ_elastic(3,:)**2))
    else
      displnorm_all = 0.d0
    endif

    displnorm_all_glob = displnorm_all
#ifdef USE_MPI
    call MPI_ALLREDUCE (displnorm_all, displnorm_all_glob, 1, MPI_DOUBLE_PRECISION, &
                      MPI_MAX, MPI_COMM_WORLD, ier)
#endif

      if (NOISE_TOMOGRAPHY /= 0) then
        if (myrank == 0) write(*,*) 'Noise simulation ', NOISE_TOMOGRAPHY, ' of 3'
      endif


    if (myrank == 0) &
      write(IOUT,*) 'Max norm of vector field in solid (elastic) = ', displnorm_all_glob

    ! check stability of the code in solid, exit if unstable
    ! negative values can occur with some compilers when the unstable value is greater
    ! than the greatest possible floating-point number of the machine
    if(displnorm_all_glob > STABILITY_THRESHOLD .or. displnorm_all_glob < 0) &
      call exit_MPI('code became unstable and blew up in solid (elastic)')

  endif

  ! poroelastic wavefield
  if(any_poroelastic_glob) then
    if(any_poroelastic) then
      displnorm_all = maxval(sqrt(displs_poroelastic(1,:)**2 &
                                + displs_poroelastic(2,:)**2))
    else
      displnorm_all = 0.d0
    endif

    displnorm_all_glob = displnorm_all
#ifdef USE_MPI
    call MPI_ALLREDUCE (displnorm_all, displnorm_all_glob, 1, MPI_DOUBLE_PRECISION, &
                      MPI_MAX, MPI_COMM_WORLD, ier)
#endif

    if (myrank == 0) &
      write(IOUT,*) 'Max norm of vector field in solid (poroelastic) = ',displnorm_all_glob

    ! check stability of the code in solid, exit if unstable
    ! negative values can occur with some compilers when the unstable value is greater
    ! than the greatest possible floating-point number of the machine
    if(displnorm_all_glob > STABILITY_THRESHOLD .or. displnorm_all_glob < 0) &
      call exit_MPI('code became unstable and blew up in solid (poroelastic)')

    if(any_poroelastic) then
      displnorm_all = maxval(sqrt(displw_poroelastic(1,:)**2 &
                                + displw_poroelastic(2,:)**2))
    else
      displnorm_all = 0.d0
    endif

    displnorm_all_glob = displnorm_all
#ifdef USE_MPI
    call MPI_ALLREDUCE (displnorm_all, displnorm_all_glob, 1, MPI_DOUBLE_PRECISION, &
                      MPI_MAX, MPI_COMM_WORLD, ier)
#endif

    if (myrank == 0) &
      write(IOUT,*) 'Max norm of vector field in fluid (poroelastic) = ',displnorm_all_glob

    ! check stability of the code in solid, exit if unstable
    ! negative values can occur with some compilers when the unstable value is greater
    ! than the greatest possible floating-point number of the machine
    if(displnorm_all_glob > STABILITY_THRESHOLD .or. displnorm_all_glob < 0) &
      call exit_MPI('code became unstable and blew up in fluid (poroelastic)')

  endif


  ! acoustic wavefield
  if(any_acoustic_glob) then
    if(any_acoustic) then
      displnorm_all = maxval(abs(potential_acoustic(:)))
    else
      displnorm_all = 0.d0
    endif

    displnorm_all_glob = displnorm_all
#ifdef USE_MPI
    call MPI_ALLREDUCE (displnorm_all, displnorm_all_glob, 1, MPI_DOUBLE_PRECISION, &
                      MPI_MAX, MPI_COMM_WORLD, ier)
#endif

    if (myrank == 0) &
      write(IOUT,*) 'Max absolute value of scalar field in fluid (acoustic) = ',displnorm_all_glob

    ! check stability of the code in fluid, exit if unstable
    ! negative values can occur with some compilers when the unstable value is greater
    ! than the greatest possible floating-point number of the machine
    if(displnorm_all_glob > STABILITY_THRESHOLD .or. displnorm_all_glob < 0) &
      call exit_MPI('code became unstable and blew up in fluid (acoustic)')

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
      write(IOUT,*) 'Elapsed time in seconds = ',tCPU
      write(IOUT,"(' Elapsed time in hh:mm:ss = ',i6,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
      write(IOUT,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it)

      ! compute estimated remaining simulation time
      t_remain = (NSTEP - it) * (tCPU/dble(it))
      int_t_remain = int(t_remain)
      ihours_remain = int_t_remain / 3600
      iminutes_remain = (int_t_remain - 3600*ihours_remain) / 60
      iseconds_remain = int_t_remain - 3600*ihours_remain - 60*iminutes_remain
      write(IOUT,*) 'Time steps remaining = ',NSTEP - it
      write(IOUT,*) 'Estimated remaining time in seconds = ',t_remain
      write(IOUT,"(' Estimated remaining time in hh:mm:ss = ',i6,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_remain,iminutes_remain,iseconds_remain

      ! compute estimated total simulation time
      t_total = t_remain + tCPU
      int_t_total = int(t_total)
      ihours_total = int_t_total / 3600
      iminutes_total = (int_t_total - 3600*ihours_total) / 60
      iseconds_total = int_t_total - 3600*ihours_total - 60*iminutes_total
      write(IOUT,*) 'Estimated total run time in seconds = ',t_total
      write(IOUT,"(' Estimated total run time in hh:mm:ss = ',i6,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_total,iminutes_total,iseconds_total

      if(it < NSTEP) then
        ! compute date and time at which the run should finish (useful for long runs)
        ! add remaining minutes and get date and time of that future timestamp in minutes
        timestamp = nint((timestamp_seconds_current + t_remain) / 60.d0)
        call invtime(timestamp,year,mon,day,hr,minutes)

        ! convert to Julian day to get day of the week
        call calndr(day,mon,year,julian_day_number)
        day_of_week = idaywk(julian_day_number)

        write(IOUT,"(' The run will finish approximately on: ',a3,' ',a3,' ',i2.2,', ',i4.4,' ',i2.2,':',i2.2)") &
            weekday_name(day_of_week),month_name(mon),day,year,hr,minutes

      endif
      write(IOUT,*)
  endif

  end subroutine check_stability

