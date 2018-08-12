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

  subroutine prepare_source_time_function()

  ! prepares source_time_function array

  use constants, only: IMAIN,ZERO,ONE,TWO,HALF,PI,QUARTER,SOURCE_DECAY_MIMIC_TRIANGLE,SOURCE_IS_MOVING,C_LDDRK,OUTPUT_FILES

  use specfem_par, only: NSTEP,NSOURCES,source_time_function, &
                         time_function_type,name_of_source_file,burst_band_width,f0_source,tshift_src,factor, &
                         t0,deltat, &
                         time_stepping_scheme,stage_time_scheme,islice_selected_source, &
                         USE_TRICK_FOR_BETTER_PRESSURE,myrank,initialfield

  implicit none

  ! local parameters
  double precision :: stf_used, timeval, DecT, Tc, omegat, omega_coa,time,coeff, t_used, Nc
  double precision :: hdur,hdur_gauss

  integer :: it,i_source,ier,num_file
  integer :: i_stage

  ! RK
  double precision, dimension(4) :: c_RK

  character(len=150) :: error_msg1 = 'Error opening the file that contains the external source: '
  character(len=250) :: error_msg
  logical :: trick_ok

  ! external functions
  double precision, external :: comp_source_time_function_heaviside_hdur
  double precision, external :: comp_source_time_function_Gaussian,comp_source_time_function_dGaussian, &
    comp_source_time_function_d2Gaussian,comp_source_time_function_d3Gaussian,marmousi_ormsby_wavelet,cos_taper
  double precision, external :: comp_source_time_function_Ricker,comp_source_time_function_d2Ricker

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Preparing source time function'
    call flush_IMAIN()
  endif

  ! checks if anything to do
  if (initialfield) then
    ! uses an initialfield
    ! dummy allocation
    allocate(source_time_function(1,1,1))
    ! we're all done
    return
  else
    allocate(source_time_function(NSOURCES,NSTEP,stage_time_scheme),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating array source_time_function')
  endif
  ! initializes stf array
  source_time_function(:,:,:) = 0.d0

  ! checks time scheme
  ! Newmark: time_stepping_scheme == 1
  ! LDDRK  : time_stepping_scheme == 2
  ! RK     : time_stepping_scheme == 3
  if (time_stepping_scheme < 1 .or. time_stepping_scheme > 3) call stop_the_code('Error invalid time stepping scheme for STF')

  ! RK time scheme constants
  if (time_stepping_scheme == 3) then
    c_RK(1) = 0.0d0 * deltat
    c_RK(2) = 0.5d0 * deltat
    c_RK(3) = 0.5d0 * deltat
    c_RK(4) = 1.0d0 * deltat
  endif

  ! user output
  if (myrank == islice_selected_source(1)) then
    if (myrank == 0) then
      write(IMAIN,*) '  saving the source time function in a text file...'
      call flush_IMAIN()
    endif
    ! opens source time file for output
    open(unit=55,file=trim(OUTPUT_FILES)//'plot_source_time_function.txt',status='unknown',iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening source time function text-file')
  endif

  ! loop on all the sources
  do i_source = 1,NSOURCES

    ! The following lines could be needed to set absolute amplitudes.
    ! In this case variables vpext,rhoext,density,poroelastcoef,assign_external_model,ispec_selected_source,kmato
    ! and double precision :: rho, cp logical :: already_done = .false. need to be introduced
    !    if (myrank == islice_selected_source(i_source)) then
    !      if (AXISYM) then
    !        if (.not. already_done) then
    !          if (  assign_external_model ) then
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

    ! checks if trick for better pressure can be applied
    if (USE_TRICK_FOR_BETTER_PRESSURE) then
      trick_ok = (time_function_type(i_source) < 4) .or. (time_function_type(i_source) == 7) .or. &
                 (time_function_type(i_source) == 9) .or. (time_function_type(i_source) == 10)

      if (.not. trick_ok) then
        call exit_MPI(myrank,'USE_TRICK_FOR_BETTER_PRESSURE is not compatible yet with the type of source you want to use!')
      endif
    endif

    do i_stage = 1,stage_time_scheme

      ! loop on all the time steps
      do it = 1,NSTEP
        ! compute current time
        select case(time_stepping_scheme)
        case (1)
          ! Newmark
          timeval = dble(it-1)*deltat
        case (2)
          ! LDDRK: Low-Dissipation and low-dispersion Runge-Kutta
          timeval = dble(it-1)*deltat + dble(C_LDDRK(i_stage))*deltat
        case (3)
          ! RK: Runge-Kutta
          timeval = dble(it-1)*deltat + c_RK(i_stage)*deltat
        case default
          call exit_MPI(myrank,'Error invalid time stepping scheme chosen, please check...')
        end select

        t_used = timeval - t0 - tshift_src(i_source)

        ! only process/partition containing source must set STF
        if (myrank == islice_selected_source(i_source) .or. SOURCE_IS_MOVING) then

          ! determines source_time_function value for different source types
          select case (time_function_type(i_source))
          case (1)
            ! Ricker: second derivative of a Gaussian
            if (USE_TRICK_FOR_BETTER_PRESSURE) then
              ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
              ! use the second derivative of the source for the source time function instead of the source itself,
              ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
              ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
              ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
              ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
              ! is accurate at second order and thus contains significantly less numerical noise.
              ! Second derivative of Ricker source time function :
              source_time_function(i_source,it,i_stage) =  - factor(i_source) * &
                       comp_source_time_function_d2Ricker(t_used,f0_source(i_source))
            else
              ! Ricker (second derivative of a Gaussian) source time function
              source_time_function(i_source,it,i_stage) = - factor(i_source) * &
                      comp_source_time_function_Ricker(t_used,f0_source(i_source))
            endif

          case (2)
            ! first derivative of a Gaussian

            if (USE_TRICK_FOR_BETTER_PRESSURE) then
              ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
              ! use the second derivative of the source for the source time function instead of the source itself,
              ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
              ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
              ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
              ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
              ! is accurate at second order and thus contains significantly less numerical noise.
              ! Third derivative of Gaussian source time function :
              source_time_function(i_source,it,i_stage) = - factor(i_source) * &
                        comp_source_time_function_d3Gaussian(t_used,f0_source(i_source))
            else
              ! First derivative of a Gaussian source time function
              source_time_function(i_source,it,i_stage) = - factor(i_source) * &
                        comp_source_time_function_dGaussian(t_used,f0_source(i_source))
            endif

          case (3,4)
            ! Gaussian/Dirac type

            if (USE_TRICK_FOR_BETTER_PRESSURE) then
              ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
              ! use the second derivative of the source for the source time function instead of the source itself,
              ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
              ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
              ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
              ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
              ! is accurate at second order and thus contains significantly less numerical noise.
              ! Second derivative of Gaussian :
              source_time_function(i_source,it,i_stage) = - factor(i_source) * &
                         comp_source_time_function_d2Gaussian(t_used,f0_source(i_source))
            else
              ! Gaussian or Dirac (we use a very thin Gaussian instead) source time function
              source_time_function(i_source,it,i_stage) = - factor(i_source) * &
                          comp_source_time_function_Gaussian(t_used,f0_source(i_source))
            endif

          case (5)
            ! Heaviside source time function (we use a very thin error function instead)
            hdur = 1.d0 / f0_source(i_source)
            hdur_gauss = hdur * 5.d0 / 3.d0

            ! convert the half duration for triangle STF to the one for Gaussian STF
            hdur_gauss = hdur_gauss / SOURCE_DECAY_MIMIC_TRIANGLE

            ! quasi-Heaviside
            source_time_function(i_source,it,i_stage) = - factor(i_source) * &
                                    comp_source_time_function_heaviside_hdur(t_used,hdur_gauss)

          case (6)
            ! ocean acoustics type I
            DecT = t0 + tshift_src(i_source)
            Tc = 4.d0 / f0_source(i_source) + DecT
            omega_coa = TWO * PI * f0_source(i_source)

            if (timeval > DecT .and. timeval < Tc) then
              ! source time function from Computational Ocean Acoustics
              omegat =  omega_coa * ( timeval - DecT )
              source_time_function(i_source,it,i_stage) = factor(i_source) * HALF * &
                    sin( omegat ) * ( ONE - cos( QUARTER * omegat ) )
              !source_time_function(i_source,it,i_stage) = factor(i_source) * HALF / omega_coa / omega_coa * &
              !      ( sin(omegat) - 8.d0 / 9.d0 * sin(3.d0/ 4.d0 * omegat) - 8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) )
            else
              source_time_function(i_source,it,i_stage) = ZERO
            endif

          case (7)
            ! ocean acoustics type II
            DecT = t0 + tshift_src(i_source)
            Tc = 4.d0 / f0_source(i_source) + DecT
            omega_coa = TWO*PI*f0_source(i_source)
            if (USE_TRICK_FOR_BETTER_PRESSURE) then !TODO
              ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
              ! use the second derivative of the source for the source time function instead of the source itself,
              ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
              ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
              ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
              ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
              ! is accurate at second order and thus contains significantly less numerical noise.
              ! Second derivative of source 7 :
              if (timeval > DecT .and. timeval < Tc) then ! t_used > 0 t_used < Nc/f0_source(i_source)) then
                source_time_function(i_source,it,i_stage) = factor(i_source) * &
                          0.5d0*(ONE-cos(omega_coa*t_used/4.0d0))*sin(omega_coa*t_used)
              else
                source_time_function(i_source,it,i_stage) = ZERO
              endif
            else
              !Tc = 1.d0 / f0_source(i_source) + DecT ! For source 1 OASES
              !if (timeval > DecT .and. timeval < Tc) then ! t_used > 0 t_used < Nc/f0_source(i_source)) then
              !  source_time_function(i_source,it,i_stage) = factor(i_source) * ( &  ! Source 1 OASES
              !            0.75d0 - cos(omega_coa*t_used) + 0.25d0*cos(TWO*omega_coa*t_used))
              !else
              !  source_time_function(i_source,it,i_stage) = ZERO
              !endif
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
!              if (timeval > DecT .and. timeval < Tc) then
!                ! source time function from Computational Ocean Acoustics
!                omegat =  omega_coa * ( timeval - DecT )
!                !source_time_function(i_source,it,i_stage) = factor(i_source) * HALF / omega_coa / omega_coa * &
!                !      ( sin(omegat) - 8.d0 / 9.d0 * sin(3.d0/ 4.d0 * omegat) - &
!                !     8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) -1./15.*( timeval - DecT ) + 1./15.*4./f0_source(i_source))
!                source_time_function(i_source,it,i_stage) = factor(i_source) * HALF / omega_coa / omega_coa * &
!                       ( - sin(omegat) + 8.d0 / 9.d0 * sin(3.d0 / 4.d0 * omegat) + &
!                        8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) - 1.d0 / 15.d0 * omegat )
!              else if (timeval > DecT) then
!                source_time_function(i_source,it,i_stage) = &
!                       - factor(i_source) * HALF / omega_coa / 15.d0 * (4.d0 / f0_source(i_source))
!              else
!                source_time_function(i_source,it,i_stage) = ZERO
!              endif
            endif

          case (8)
            ! external type
            ! opens external file to read in source time function
            if (it == 1) then
              coeff = factor(i_source)
              open(unit=num_file,file=trim(name_of_source_file(i_source)),status='old',action='read',iostat=ier)
              if (ier /= 0) then
                error_msg = trim(error_msg1)//trim(name_of_source_file(i_source))
                call exit_MPI(myrank,error_msg)
              endif
            endif

            ! format: #time #stf-value
            read(num_file,*) time, source_time_function(i_source,it,i_stage)

            ! closes external file
            if (it == NSTEP ) close(num_file)

          case (9)
            ! burst type
            DecT = t0 + tshift_src(i_source)
            t_used = (timeval-t0-tshift_src(i_source))
            Nc = TWO * f0_source(i_source) / burst_band_width(i_source)
            Tc = Nc / f0_source(i_source) + DecT
            omega_coa = TWO*PI*f0_source(i_source)

            if (USE_TRICK_FOR_BETTER_PRESSURE) then
              ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
              ! use the second derivative of the source for the source time function instead of the source itself,
              ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
              ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
              ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
              ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
              ! is accurate at second order and thus contains significantly less numerical noise.
              ! Second derivative of Burst :
              if (timeval > DecT .and. timeval < Tc) then ! t_used > 0 t_used < Nc/f0_source(i_source)) then
                source_time_function(i_source,it,i_stage) = - factor(i_source) * (0.5d0 * (omega_coa)**2 * &
                          sin(omega_coa*t_used) * cos(omega_coa*t_used/Nc) / Nc**2 - &
                          0.5d0 * (omega_coa)**2 * sin(omega_coa*t_used) * &
                          (ONE-cos(omega_coa*t_used/Nc)) + &
                          (omega_coa)**2 * cos(omega_coa*t_used) * &
                          sin(omega_coa*t_used/Nc) / Nc)
              !else if (timeval > DecT) then
              !  source_time_function(i_source,it,i_stage) = ZERO
              else
                source_time_function(i_source,it,i_stage) = ZERO
              endif
              ! Integral of burst
              !if (timeval > DecT .and. timeval < Tc) then ! t_used > 0 t_used < Nc/f0_source(i_source)) then
              !  source_time_function(i_source,it,i_stage) = - factor(i_source) * ( &
              !  Nc*( (Nc+1.0d0)*cos((omega_coa*(Nc-1.0d0)*t_used)/Nc) + &
              !       (Nc-1.0d0)*cos((omega_coa*(Nc+1.0d0)*t_used)/Nc)) - &
              !  TWO*(Nc**2-1.0d0)*cos(omega_coa*t_used) &
              !  ) / (8.0d0*PI*f0_source(i_source)*(Nc-1)*(Nc+1))
              !else
              !  source_time_function(i_source,it,i_stage) = ZERO
              !endif
              ! Double integral of burst
              !if (timeval > DecT .and. timeval < Tc) then ! t_used > 0 t_used < Nc/f0_source(i_source)) then
              !  source_time_function(i_source,it,i_stage) = - factor(i_source) * ( &
              !      -sin(TWO*f0_source(i_source)*Pi*t_used)/(8.0d0*f0_source(i_source)**TWO*Pi**2) + &
              !      (Nc**2*sin((TWO*f0_source(i_source)*(Nc-1)*PI*t_used)/Nc))/(16.0d0*f0_source(i_source)**2*(Nc-1)**2*Pi**2) + &
              !      (Nc**2*sin((TWO*f0_source(i_source)*(Nc+1)*PI*t_used)/Nc))/(16.0d0*f0_source(i_source)**2*(Nc+1)**2*Pi**2) )
              !else
              !  source_time_function(i_source,it,i_stage) = ZERO
              !endif
            else
              if (timeval > DecT .and. timeval < Tc) then ! t_used > 0 t_used < Nc/f0_source(i_source)) then
                source_time_function(i_source,it,i_stage) = - factor(i_source) * &
                          0.5d0*(ONE-cos(omega_coa*t_used/Nc))*sin(omega_coa*t_used)
              !else if (timeval > DecT) then
              !  source_time_function(i_source,it,i_stage) = ZERO
              else
                source_time_function(i_source,it,i_stage) = ZERO
              endif
            endif

          case (10)
            ! Sinus source time function
            omega_coa = TWO*PI*f0_source(i_source)
            if (USE_TRICK_FOR_BETTER_PRESSURE) then
              ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
              ! use the second derivative of the source for the source time function instead of the source itself,
              ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
              ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
              ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
              ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
              ! is accurate at second order and thus contains significantly less numerical noise.
              ! Third derivative of Gaussian source time function :
              source_time_function(i_source,it,i_stage) = -TWO*PI*omega_coa*f0_source(i_source)* &
                                                          factor(i_source) * sin(omega_coa*t_used)
            else
              ! First derivative of a Gaussian source time function
              source_time_function(i_source,it,i_stage) = factor(i_source) * sin(omega_coa*t_used)
            endif

          case (11)
              ! Marmousi_ormsby_wavelet
              hdur = 1.0 / 35.0
              source_time_function(i_source,it,i_stage) = factor(i_source) * &
                        cos_taper(t_used,hdur) * marmousi_ormsby_wavelet(PI*t_used)

          case default
            call exit_MPI(myrank,'unknown source time function')

          end select

          ! outputs (first source) source time function
          if (i_source == 1 .and. i_stage == 1) then
            stf_used = source_time_function(i_source,it,i_stage)

            ! note: earliest start time of the simulation is: (it-1)*deltat - t0 - tshift_src(i_source)
            if (myrank == islice_selected_source(1)) write(55,*) timeval-t0,' ',stf_used

          endif
        endif
      enddo
    enddo
  enddo

  ! closes STF file
  if (myrank == islice_selected_source(1)) close(55)

  !print *,"myrank:",myrank,"Ok"

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_source_time_function

