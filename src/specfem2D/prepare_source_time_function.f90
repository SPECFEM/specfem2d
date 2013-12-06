
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
! Copyright CNRS, INRIA and University of Pau, France,
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


  subroutine prepare_source_time_function(myrank,NSTEP,NSOURCES,source_time_function, &
                          time_function_type,f0,tshift_src,factor,aval, &
                          t0,nb_proc_source,deltat,stage_time_scheme,c_LDDRK)

! prepares source_time_function array

  implicit none
  include "constants.h"

  integer :: myrank,NSTEP

  integer :: NSOURCES
  integer, dimension(NSOURCES) :: time_function_type
  double precision, dimension(NSOURCES) :: f0,tshift_src,factor
  double precision, dimension(NSOURCES) :: aval
  double precision :: t0
  integer,dimension(NSOURCES) :: nb_proc_source
  double precision :: deltat
  integer :: stage_time_scheme
  double precision, dimension(Nstages) :: c_LDDRK

  real(kind=CUSTOM_REAL),dimension(NSOURCES,NSTEP,stage_time_scheme) :: source_time_function

  ! local parameters
  double precision :: stf_used, time, DecT, Tc, omegat, omega_coa
  double precision, dimension(NSOURCES) :: hdur,hdur_gauss
  double precision, external :: netlib_specfun_erf
  integer :: it,i_source
  integer :: i_stage
  double precision, dimension(4) :: c_RK

  if(stage_time_scheme == 4)then
   c_RK(1)=0.0d0*deltat
   c_RK(2)=0.5d0*deltat
   c_RK(3)=0.5d0*deltat
   c_RK(4)=1.0d0*deltat
  endif

  ! user output
  if (myrank == 0) then
    write(IOUT,*)
    write(IOUT,*) 'Saving the source time function in a text file...'
    write(IOUT,*)
    open(unit=55,file='OUTPUT_FILES/source.txt',status='unknown')
  endif

  !    ! loop on all the sources
  !    do i_source=1,NSOURCES

  ! loop on all the time steps


  do it = 1,NSTEP

    ! note: t0 is the simulation start time, tshift_src is the time shift of the source
    !          relative to this start time
    do i_stage = 1,stage_time_scheme

    ! compute current time
    if(stage_time_scheme == 1)then
    time = (it-1)*deltat
    endif

    if(stage_time_scheme == 4)then
     time = (it-1)*deltat+c_RK(i_stage)*deltat
    endif

    if(stage_time_scheme == 6)then
     time = (it-1)*deltat+c_LDDRK(i_stage)*deltat
    endif


    stf_used = 0.d0

    ! loop on all the sources
    do i_source=1,NSOURCES

      if( time_function_type(i_source) == 1 ) then

        ! Ricker (second derivative of a Gaussian) source time function
        source_time_function(i_source,it,i_stage) = - factor(i_source) * &
                  (ONE-TWO*aval(i_source)*(time-t0-tshift_src(i_source))**2) * &
                  exp(-aval(i_source)*(time-t0-tshift_src(i_source))**2)

        ! source_time_function(i_source,it) = - factor(i_source) *  &
        !               TWO*aval(i_source)*sqrt(aval(i_source))*&
        !               (time-t0-tshift_src(i_source))/pi * exp(-aval(i_source)*(time-t0-tshift_src(i_source))**2)

      else if( time_function_type(i_source) == 2 ) then

        ! first derivative of a Gaussian source time function
        source_time_function(i_source,it,i_stage) = - factor(i_source) * &
                  TWO*aval(i_source)*(time-t0-tshift_src(i_source)) * &
                  exp(-aval(i_source)*(time-t0-tshift_src(i_source))**2)

      else if(time_function_type(i_source) == 3 .or. time_function_type(i_source) == 4) then

        ! Gaussian or Dirac (we use a very thin Gaussian instead) source time function
        source_time_function(i_source,it,i_stage) = factor(i_source) * &
                  exp(-aval(i_source)*(time-t0-tshift_src(i_source))**2)

      else if(time_function_type(i_source) == 5) then

        ! Heaviside source time function (we use a very thin error function instead)
        hdur(i_source) = 1.d0 / f0(i_source)
        hdur_gauss(i_source) = hdur(i_source) * 5.d0 / 3.d0
        source_time_function(i_source,it,i_stage) = factor(i_source) * 0.5d0*(1.0d0 + &
            netlib_specfun_erf(SOURCE_DECAY_MIMIC_TRIANGLE*(time-t0-tshift_src(i_source))/hdur_gauss(i_source)))

      else if(time_function_type(i_source) == 6) then

        DecT = t0 + tshift_src(i_source)

        Tc = 4.d0 / f0(i_source) + DecT

        if ( time > DecT .and. time < Tc ) then

           ! source time function from Computational Ocean Acoustics
           omega_coa = TWO * PI * f0(i_source)
           omegat =  omega_coa * ( time - DecT )
           source_time_function(i_source,it,i_stage) = factor(i_source) * HALF * &
                 sin( omegat ) * ( ONE - cos( QUART * omegat ) )
           !source_time_function(i_source,it,i_stage) = - factor(i_source) * HALF / omega_coa / omega_coa * &
           !      ( sin(omegat) - 8.d0 / 9.d0 * sin(3.d0/ 4.d0 * omegat) - 8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) )

        else

           source_time_function(i_source,it,i_stage) = ZERO

        endif

       else if(time_function_type(i_source) == 7) then

        DecT = t0 + tshift_src(i_source)
        Tc = 4.d0 / f0(i_source) + DecT
        omega_coa = TWO * PI * f0(i_source)

        if ( time > DecT .and. time < Tc ) then
          ! source time function from Computational Ocean Acoustics
           omegat =  omega_coa * ( time - DecT )
           !source_time_function(i_source,it,i_stage) = factor(i_source) * HALF / omega_coa / omega_coa * &
           !      ( sin(omegat) - 8.d0 / 9.d0 * sin(3.d0/ 4.d0 * omegat) - &
           !	   8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) -1./15.*( time - DecT ) + 1./15.*4./f0(i_source))

           source_time_function(i_source,it,i_stage) = factor(i_source) * HALF / omega_coa / omega_coa * &
                 ( - sin(omegat) + 8.d0 / 9.d0 * sin(3.d0 / 4.d0 * omegat) + &
                  8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) - 1.d0 / 15.d0 * omegat )

        else if ( time > DecT ) then

           source_time_function(i_source,it,i_stage) = - factor(i_source) * HALF / omega_coa / 15.d0 * (4.d0 / f0(i_source))

        else

           source_time_function(i_source,it,i_stage) = ZERO

        endif


      else
        call exit_MPI('unknown source time function')
      endif

      stf_used = stf_used + source_time_function(i_source,it,i_stage)

      ! output relative time in third column, in case user wants to check it as well
      ! if (myrank == 0 .and. i_source == 1) write(55,*) sngl(time-t0-tshift_src(1)),real(source_time_function(1,it),4),sngl(time)
      if (myrank == 0 .and. i_source == 1 .and. i_stage == 1) then
          ! note: earliest start time of the simulation is: (it-1)*deltat - t0
          write(55,*) sngl(time-t0),sngl(stf_used),sngl(time)
      endif

    enddo


    enddo

  enddo

  if (myrank == 0) close(55)

  ! nb_proc_source is the number of processes that own the source (the nearest point). It can be greater
  ! than one if the nearest point is on the interface between several partitions with an explosive source.
  ! since source contribution is linear, the source_time_function is cut down by that number (it would have been similar
  ! if we just had elected one of those processes).
  do i_source=1,NSOURCES
    source_time_function(i_source,:,:) = source_time_function(i_source,:,:) / nb_proc_source(i_source)
  enddo

  end subroutine prepare_source_time_function
