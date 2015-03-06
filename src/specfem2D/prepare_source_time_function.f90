
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


  subroutine prepare_source_time_function()

! prepares source_time_function array

  use specfem_par, only: NSTEP,NSOURCES,source_time_function, &
                         time_function_type,f0,tshift_src,factor,aval, &
                         t0,nb_proc_source,deltat,stage_time_scheme,c_LDDRK,is_proc_source, &
                         USE_TRICK_FOR_BETTER_PRESSURE

  implicit none
  include "constants.h"

  ! local parameters
  double precision :: stf_used, timeval, DecT, Tc, omegat, omega_coa,time,coeff, t_used
  double precision, dimension(NSOURCES) :: hdur,hdur_gauss
  double precision, external :: netlib_specfun_erf
  integer :: it,i_source,ier,num_file
  integer :: i_stage
  double precision, dimension(4) :: c_RK
  character(len=256) :: name_of_file

  if(stage_time_scheme == 4)then
   c_RK(1)=0.0d0*deltat
   c_RK(2)=0.5d0*deltat
   c_RK(3)=0.5d0*deltat
   c_RK(4)=1.0d0*deltat
  endif

  ! user output
  if (is_proc_source(1) == 1) then
    write(IOUT,*)
    write(IOUT,*) 'Saving the source time function in a text file...'
    write(IOUT,*)
    open(unit=55,file='OUTPUT_FILES/source.txt',status='unknown')
  endif

    ! loop on all the sources
    do i_source=1,NSOURCES

    num_file = 800 + i_source

    ! note: t0 is the simulation start time, tshift_src is the time shift of the source
    !          relative to this start time
    
    if (time_function_type(i_source) >= 5 .and. USE_TRICK_FOR_BETTER_PRESSURE) then 
      call exit_MPI('USE_TRICK_FOR_BETTER_PRESSURE is not compatible yet with the type of source you want to use')
    endif

    do i_stage = 1,stage_time_scheme

 ! loop on all the time steps
    do it=1,NSTEP

! compute current time
    if(stage_time_scheme == 1) timeval = (it-1)*deltat

    if(stage_time_scheme == 4) timeval = (it-1)*deltat+c_RK(i_stage)*deltat

    if(stage_time_scheme == 6) timeval = (it-1)*deltat+c_LDDRK(i_stage)*deltat
    
    t_used =(timeval-t0-tshift_src(i_source))

    stf_used = 0.d0

    if(is_proc_source(i_source) == 1) then

      if( time_function_type(i_source) == 1 ) then
  
        if (USE_TRICK_FOR_BETTER_PRESSURE) then 
          ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
          ! use the second derivative of the source for the source time function instead of the source itself,
          ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
          ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
          ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
          ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
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

      else if( time_function_type(i_source) == 2 ) then
        if (USE_TRICK_FOR_BETTER_PRESSURE) then 
          ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
          ! use the second derivative of the source for the source time function instead of the source itself,
          ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
          ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
          ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
          ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
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

      else if(time_function_type(i_source) == 3 .or. time_function_type(i_source) == 4) then
        if (USE_TRICK_FOR_BETTER_PRESSURE) then 
          ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
          ! use the second derivative of the source for the source time function instead of the source itself,
          ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
          ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
          ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
          ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
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

      else if(time_function_type(i_source) == 5) then

        ! Heaviside source time function (we use a very thin error function instead)
        hdur(i_source) = 1.d0 / f0(i_source)
        hdur_gauss(i_source) = hdur(i_source) * 5.d0 / 3.d0
        source_time_function(i_source,it,i_stage) = factor(i_source) * 0.5d0*(1.0d0 + &
            netlib_specfun_erf(SOURCE_DECAY_MIMIC_TRIANGLE*t_used/hdur_gauss(i_source)))

      else if(time_function_type(i_source) == 6) then

        DecT = t0 + tshift_src(i_source)

        Tc = 4.d0 / f0(i_source) + DecT

        if ( timeval > DecT .and. timeval < Tc ) then

           ! source time function from Computational Ocean Acoustics
           omega_coa = TWO * PI * f0(i_source)
           omegat =  omega_coa * ( timeval - DecT )
           source_time_function(i_source,it,i_stage) = factor(i_source) * HALF * &
                 sin( omegat ) * ( ONE - cos( QUARTER * omegat ) )
           !source_time_function(i_source,it,i_stage) = - factor(i_source) * HALF / omega_coa / omega_coa * &
           !      ( sin(omegat) - 8.d0 / 9.d0 * sin(3.d0/ 4.d0 * omegat) - 8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) )

        else

           source_time_function(i_source,it,i_stage) = ZERO

        endif

       else if(time_function_type(i_source) == 7) then

        DecT = t0 + tshift_src(i_source)
        Tc = 4.d0 / f0(i_source) + DecT
        omega_coa = TWO * PI * f0(i_source)

        if ( timeval > DecT .and. timeval < Tc ) then
          ! source time function from Computational Ocean Acoustics
           omegat =  omega_coa * ( timeval - DecT )
           !source_time_function(i_source,it,i_stage) = factor(i_source) * HALF / omega_coa / omega_coa * &
           !      ( sin(omegat) - 8.d0 / 9.d0 * sin(3.d0/ 4.d0 * omegat) - &
           !     8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) -1./15.*( timeval - DecT ) + 1./15.*4./f0(i_source))

           source_time_function(i_source,it,i_stage) = factor(i_source) * HALF / omega_coa / omega_coa * &
                 ( - sin(omegat) + 8.d0 / 9.d0 * sin(3.d0 / 4.d0 * omegat) + &
                  8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) - 1.d0 / 15.d0 * omegat )

        else if ( timeval > DecT ) then

           source_time_function(i_source,it,i_stage) = - factor(i_source) * HALF / omega_coa / 15.d0 * (4.d0 / f0(i_source))

        else

           source_time_function(i_source,it,i_stage) = ZERO

        endif

       else if(time_function_type(i_source) == 8) then

        if (it == 1 ) then
          coeff = factor(i_source)
          write(name_of_file,"(a,i3.3,a)") 'DATA/source_custom',int(coeff),'.txt'
          open(unit=num_file,file=name_of_file,iostat=ier)
          if( ier /= 0 ) call exit_MPI('error opening file source_custom*****')
        endif

        read(num_file,*) time, source_time_function(i_source,it,i_stage)

        if (it == NSTEP ) close(num_file)

      else
        call exit_MPI('unknown source time function')
      endif

      stf_used = stf_used + source_time_function(i_source,it,i_stage)

      ! output relative time in third column, in case user wants to check it as well
! if (myrank == 0 .and. i_source == 1) write(55,*) sngl(timeval-t0-tshift_src(1)),real(source_time_function(1,it),4),sngl(timeval)
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
  do i_source=1,NSOURCES
    source_time_function(i_source,:,:) = source_time_function(i_source,:,:) / nb_proc_source(i_source)
  enddo

  end subroutine prepare_source_time_function
