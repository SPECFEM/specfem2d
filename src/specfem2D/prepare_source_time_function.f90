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


  subroutine prepare_source_time_function()

  ! prepares source_time_function array

  use constants,only: IMAIN,ZERO,ONE,TWO,HALF,PI,QUARTER,SOURCE_DECAY_MIMIC_TRIANGLE

  use specfem_par, only: AXISYM,NSTEP,NSOURCES,source_time_function, &
                         time_function_type,name_of_source_file,burst_band_width,f0_source,tshift_src,factor, &
                         aval,t0,nb_proc_source,deltat,stage_time_scheme,C_LDDRK,is_proc_source, &
                         USE_TRICK_FOR_BETTER_PRESSURE,myrank

  implicit none

  ! local parameters
  double precision :: stf_used, timeval, DecT, Tc, omegat, omega_coa,time,coeff, t_used, Nc
  double precision, dimension(NSOURCES) :: hdur,hdur_gauss
  double precision, external :: netlib_specfun_erf
  integer :: it,i_source,ier,num_file
  integer :: i_stage
  double precision, dimension(4) :: c_RK
  character(len=27) :: error_msg1 = 'Error opening file source: '
  character(len=177) :: error_msg

  if (stage_time_scheme == 4) then
    c_RK(1) = 0.0d0 * deltat
    c_RK(2) = 0.5d0 * deltat
    c_RK(3) = 0.5d0 * deltat
    c_RK(4) = 1.0d0 * deltat
  endif

  ! user output
  if (is_proc_source(1) == 1) then
    write(IMAIN,*)
    write(IMAIN,*) 'Saving the source time function in a text file...'
    write(IMAIN,*)
    call flush_IMAIN()
    ! opens source time file for output
    open(unit=55,file='OUTPUT_FILES/plot_source_time_function.txt',status='unknown',iostat=ier)
    if (ier /= 0) stop 'Error opening source time function text-file'
  endif

  ! loop on all the sources
  do i_source = 1,NSOURCES

    if (AXISYM) then
      factor(i_source) = - factor(i_source)
    endif
    ! The following lines could be needed to set absolute amplitudes.
    ! In this case variables vpext,rhoext,density,poroelastcoef,assign_external_model,ispec_selected_source,kmato
    ! and double precision :: rho, cp logical :: already_done = .false. need to be introduced
    !    if(is_proc_source(i_source) == 1) then
    !      if (AXISYM) then
    !        if (.not. already_done) then
    !          if(  assign_external_model ) then
    !            cp = vpext(0,0,ispec_selected_source(i_source))
    !            TODO (above): We must interpolate to find the exact cp value at source location
    !          else
    !            rho = density(1,kmato(ispec_selected_source(i_source)))
    !            cp = sqrt(poroelastcoef(3,1,kmato(ispec_selected_source(i_source)))/rho)
    !          endif
    !          factor(i_source) = - factor(i_source)*2.0d0*cp**2*0.45d-5 !0.225d-5
    !          if (time_function_type (i_source) == 7)  factor(i_source) = factor(i_source) * 222066.1d0 !444132.2d0
    !          already_done = .true.
    !        endif
    !      endif
    !    endif

    num_file = 800 + i_source

    ! note: t0 is the simulation start time, tshift_src is the time shift of the source
    !          relative to this start time
    if (time_function_type(i_source) >= 5 .and. USE_TRICK_FOR_BETTER_PRESSURE) then
      call exit_MPI(myrank,'USE_TRICK_FOR_BETTER_PRESSURE is not compatible yet with the type of source you want to use')
    endif

    do i_stage = 1,stage_time_scheme
      ! loop on all the time steps
      do it = 1,NSTEP

        ! compute current time
        if (stage_time_scheme == 1) timeval = (it-1)*deltat
        if (stage_time_scheme == 4) timeval = (it-1)*deltat+c_RK(i_stage)*deltat
        if (stage_time_scheme == 6) timeval = (it-1)*deltat+C_LDDRK(i_stage)*deltat

        t_used = timeval - t0 - tshift_src(i_source)
        stf_used = 0.d0

        if (is_proc_source(i_source) == 1) then
          if (time_function_type(i_source) == 1) then
            ! Ricker type: second derivative
            if (USE_TRICK_FOR_BETTER_PRESSURE) then
              ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
              ! use the second derivative of the source for the source time function instead of the source itself,
              ! and then record -minus_int_int_pressure_acoustic() as pressure seismograms instead of -minus_pressure_acoustic();
              ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
              ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
              ! thus in fluid elements minus_pressure_acoustic() is accurate at zeroth order while minus_int_int_pressure_acoustic()
              ! is accurate at second order and thus contains significantly less numerical noise.
              ! Second derivative of Ricker source time function :
              source_time_function(i_source,it,i_stage) = factor(i_source) * &
                       2.0d0*aval(i_source) * (3.0d0 - 12.0d0*aval(i_source)*t_used**2 + 4.0d0*aval(i_source)**2*t_used**4) * &
                       exp(-aval(i_source)*t_used**2)
            else
              ! Ricker (second derivative of a Gaussian) source time function
              source_time_function(i_source,it,i_stage) = - factor(i_source) * &
                         (ONE-TWO*aval(i_source)*t_used**2) * &
                         exp(-aval(i_source)*t_used**2)

              ! source_time_function(i_source,it) = - factor(i_source) *  &
              !               TWO*aval(i_source)*sqrt(aval(i_source))*&
              !               t_used/pi * exp(-aval(i_source)*t_used**2)
            endif

          else if (time_function_type(i_source) == 2) then
            ! Ricker type: first derivative
            if (USE_TRICK_FOR_BETTER_PRESSURE) then
              ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
              ! use the second derivative of the source for the source time function instead of the source itself,
              ! and then record -minus_int_int_pressure_acoustic() as pressure seismograms instead of -minus_pressure_acoustic();
              ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
              ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
              ! thus in fluid elements minus_pressure_acoustic() is accurate at zeroth order while minus_int_int_pressure_acoustic()
              ! is accurate at second order and thus contains significantly less numerical noise.
              ! Third derivative of Gaussian source time function :
              source_time_function(i_source,it,i_stage) = factor(i_source) * &
                        4.0d0*aval(i_source)**2*t_used * (3.0d0 - 2.0d0*aval(i_source)*t_used**2) * &
                        exp(-aval(i_source)*t_used**2)
            else
              ! First derivative of a Gaussian source time function
              source_time_function(i_source,it,i_stage) = - factor(i_source) * &
                        TWO*aval(i_source)*t_used * &
                        exp(-aval(i_source)*t_used**2)
            endif

          else if (time_function_type(i_source) == 3 .or. time_function_type(i_source) == 4) then
            ! Gaussian/Dirac type
            if (USE_TRICK_FOR_BETTER_PRESSURE) then
              ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
              ! use the second derivative of the source for the source time function instead of the source itself,
              ! and then record -minus_int_int_pressure_acoustic() as pressure seismograms instead of -minus_pressure_acoustic();
              ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
              ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
              ! thus in fluid elements minus_pressure_acoustic() is accurate at zeroth order while minus_int_int_pressure_acoustic()
              ! is accurate at second order and thus contains significantly less numerical noise.
              ! Second derivative of Gaussian :
              source_time_function(i_source,it,i_stage) = factor(i_source) * &
                         2.0d0 * aval(i_source) * (2.0d0 * aval(i_source) * t_used**2 - 1.0d0) * &
                         exp(-aval(i_source)*t_used**2)
            else
              ! Gaussian or Dirac (we use a very thin Gaussian instead) source time function
              source_time_function(i_source,it,i_stage) = factor(i_source) * &
                        exp(-aval(i_source)*t_used**2)
            endif

          else if (time_function_type(i_source) == 5) then
            ! Heaviside source time function (we use a very thin error function instead)
            hdur(i_source) = 1.d0 / f0_source(i_source)
            hdur_gauss(i_source) = hdur(i_source) * 5.d0 / 3.d0
            source_time_function(i_source,it,i_stage) = factor(i_source) * 0.5d0*(1.0d0 + &
                netlib_specfun_erf(SOURCE_DECAY_MIMIC_TRIANGLE*t_used/hdur_gauss(i_source)))

          else if (time_function_type(i_source) == 6) then
            ! ocean acoustics type I
            DecT = t0 + tshift_src(i_source)
            Tc = 4.d0 / f0_source(i_source) + DecT
            omega_coa = TWO * PI * f0_source(i_source)

            if (timeval > DecT .and. timeval < Tc) then
              ! source time function from Computational Ocean Acoustics
              omegat =  omega_coa * ( timeval - DecT )
              source_time_function(i_source,it,i_stage) = factor(i_source) * HALF * &
                    sin( omegat ) * ( ONE - cos( QUARTER * omegat ) )
              !source_time_function(i_source,it,i_stage) = - factor(i_source) * HALF / omega_coa / omega_coa * &
              !      ( sin(omegat) - 8.d0 / 9.d0 * sin(3.d0/ 4.d0 * omegat) - 8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) )
            else
              source_time_function(i_source,it,i_stage) = ZERO
            endif

          else if (time_function_type(i_source) == 7) then
            ! ocean acoustics type II
            DecT = t0 + tshift_src(i_source)
            Tc = 4.d0 / f0_source(i_source) + DecT
            omega_coa = TWO * PI * f0_source(i_source)

            if (timeval > DecT .and. timeval < Tc) then
              ! source time function from Computational Ocean Acoustics
              omegat =  omega_coa * ( timeval - DecT )
              !source_time_function(i_source,it,i_stage) = factor(i_source) * HALF / omega_coa / omega_coa * &
              !      ( sin(omegat) - 8.d0 / 9.d0 * sin(3.d0/ 4.d0 * omegat) - &
              !     8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) -1./15.*( timeval - DecT ) + 1./15.*4./f0_source(i_source))
              source_time_function(i_source,it,i_stage) = factor(i_source) * HALF / omega_coa / omega_coa * &
                     ( - sin(omegat) + 8.d0 / 9.d0 * sin(3.d0 / 4.d0 * omegat) + &
                      8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) - 1.d0 / 15.d0 * omegat )
            else if (timeval > DecT) then
              source_time_function(i_source,it,i_stage) = &
              - factor(i_source) * HALF / omega_coa / 15.d0 * (4.d0 / f0_source(i_source))
            else
              source_time_function(i_source,it,i_stage) = ZERO
            endif

          else if (time_function_type(i_source) == 8) then
            ! external type
            ! opens external file to read in source time function
            if (it == 1) then
              coeff = factor(i_source)
              error_msg = error_msg1//name_of_source_file(i_source)
              open(unit=num_file,file=name_of_source_file(i_source),iostat=ier)
              if (ier /= 0 ) call exit_MPI(myrank,error_msg)
            endif

            ! format: #time #stf-value
            read(num_file,*) time, source_time_function(i_source,it,i_stage)

            ! closes external file
            if (it == NSTEP ) close(num_file)

          else if (time_function_type(i_source) == 9) then
            ! burst type
            DecT = t0 + tshift_src(i_source)
            t_used = (timeval-t0-tshift_src(i_source))
            Nc = TWO * f0_source(i_source) / burst_band_width(i_source)
            Tc = Nc / f0_source(i_source) + DecT
            if (timeval > DecT .and. timeval < Tc) then ! t_used > 0 t_used < Nc/f0_source(i_source)) then
              source_time_function(i_source,it,i_stage) = - factor(i_source) * &
                        0.5d0*(ONE-cos(TWO*PI*f0_source(i_source)*t_used/Nc))*sin(TWO*PI*f0_source(i_source)*t_used)
            !else if (timeval > DecT) then
            !  source_time_function(i_source,it,i_stage) = ZERO
            else
              source_time_function(i_source,it,i_stage) = ZERO
            endif
          else
            call exit_MPI(myrank,'unknown source time function')
          endif

          stf_used = stf_used + source_time_function(i_source,it,i_stage)

          ! output relative time in third column, in case user wants to check it as well
          !if (myrank == 0 .and. i_source == 1) write(55,*) &
          !  sngl(timeval-t0-tshift_src(1)),real(source_time_function(1,it),4),sngl(timeval)
          if (i_source == 1 .and. i_stage == 1) then
            ! note: earliest start time of the simulation is: (it-1)*deltat - t0
            write(55,*) sngl(timeval-t0),sngl(stf_used),sngl(timeval)
          endif
        endif
      enddo
    enddo
  enddo

  if (is_proc_source(1) == 1) close(55)

  ! nb_proc_source is the number of processes that own the source (the nearest point). It can be greater
  ! than one if the nearest point is on the interface between several partitions with an explosive source.
  ! since source contribution is linear, the source_time_function is cut down by that number (it would have been similar
  ! if we just had elected one of those processes).
  do i_source = 1,NSOURCES
    source_time_function(i_source,:,:) = source_time_function(i_source,:,:) / nb_proc_source(i_source)
  enddo

  end subroutine prepare_source_time_function

