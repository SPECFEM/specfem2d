!========================================================================
!
!                            S P E C F E M 2 D
!                            -----------------
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

  use constants, only: CUSTOM_REAL,IMAIN,ZERO,ONE,TWO,HALF,PI,QUARTER, &
                       SOURCE_DECAY_MIMIC_TRIANGLE, &
                       C_LDDRK,C_RK4,ALPHA_SYMPLECTIC

  use specfem_par, only: NSTEP, NSOURCES, source_time_function, &
                         time_function_type, burst_band_width, f0_source,tshift_src, &
                         factor, t0, DT, SOURCE_IS_MOVING, &
                         time_stepping_scheme, NSTAGE_TIME_SCHEME, islice_selected_source, &
                         USE_TRICK_FOR_BETTER_PRESSURE, myrank, initialfield, PRINT_SOURCE_TIME_FUNCTION

  implicit none

  ! local parameters
  double precision :: timeval, t_used
  double precision :: hdur, hdur_gauss, f0, f0_sampling
  double precision :: stf

  integer :: it,isource,ier
  integer :: i_stage
  logical :: trick_ok

  ! external functions
  double precision, external :: comp_source_time_function_heaviside_hdur
  double precision, external :: comp_source_time_function_Gaussian_norm,comp_source_time_function_d2Gaussian_norm
  double precision, external :: comp_source_time_function_Gaussian,comp_source_time_function_dGaussian, &
                                comp_source_time_function_d2Gaussian,comp_source_time_function_d3Gaussian
  double precision, external :: comp_source_time_function_marmousi
  double precision, external :: comp_source_time_function_Ricker,comp_source_time_function_d2Ricker
  double precision, external :: comp_source_time_function_mono,comp_source_time_function_d2mono
  double precision, external :: comp_source_time_function_ocean_I
  double precision, external :: comp_source_time_function_ocean_II,comp_source_time_function_d2ocean_II
  double precision, external :: comp_source_time_function_ext
  double precision, external :: comp_source_time_function_burst,comp_source_time_function_d2burst
  double precision, external :: comp_source_time_function_Brune,comp_source_time_function_Smooth_Brune
  double precision, external :: comp_source_time_function_Yoffe,comp_source_time_function_Yoffe_integrated
  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Preparing source time function'
    call flush_IMAIN()
  endif

  ! Newmark: time_stepping_scheme == 1
  ! LDDRK  : time_stepping_scheme == 2
  ! RK     : time_stepping_scheme == 3
  ! PEFRL  : time_stepping_scheme == 4
  ! user output
  select case(time_stepping_scheme)
  case (1)
    ! Newmark
    if (myrank == 0) write(IMAIN,*) '  time stepping scheme:   Newmark'
  case (2)
    ! LDDRK
    if (myrank == 0) write(IMAIN,*) '  time stepping scheme:   LDDRK'
  case (3)
    ! RK4
    if (myrank == 0) write(IMAIN,*) '  time stepping scheme:   RK4'
  case (4)
    ! symplectic PEFRL
    if (myrank == 0) write(IMAIN,*) '  time stepping scheme:   symplectic PEFRL'
  case default
    call stop_the_code('Error invalid time stepping scheme for STF')
  end select

  if (myrank == 0) then
    write(IMAIN,*) '  time stepping stages: ',NSTAGE_TIME_SCHEME
    write(IMAIN,*) '  time step size      : ',sngl(DT)
    write(IMAIN,*)
    write(IMAIN,*) '  number of time steps: ',NSTEP
    if (initialfield) then
      write(IMAIN,*) '  initital field      : ',initialfield
    else
      write(IMAIN,*) '  number of sources   : ',NSOURCES
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! checks if anything to do
  if (initialfield) then
    ! uses an initialfield
    ! dummy allocation
    allocate(source_time_function(1,1,1))
    ! we're all done
    return
  endif

  ! checks if trick for better pressure can be applied
  do isource = 1,NSOURCES
    if (USE_TRICK_FOR_BETTER_PRESSURE) then
      trick_ok = (time_function_type(isource) < 4) .or. (time_function_type(isource) == 7) .or. &
                 (time_function_type(isource) == 9) .or. (time_function_type(isource) == 10)
      if (.not. trick_ok) then
        print *,'Error: source ',isource,' has invalid source time function type ',time_function_type(isource), &
                ' for USE_TRICK_FOR_BETTER_PRESSURE'
        call exit_MPI(myrank,'USE_TRICK_FOR_BETTER_PRESSURE is not compatible yet with the type of source you want to use!')
      endif
    endif
  enddo

  ! initializes stf array
  allocate(source_time_function(NSOURCES,NSTEP,NSTAGE_TIME_SCHEME),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating array source_time_function')
  source_time_function(:,:,:) = 0.0_CUSTOM_REAL

  ! loop over all the sources
  do isource = 1,NSOURCES

    ! note: t0 is the simulation start time, tshift_src is the time shift of the source
    !          relative to this start time

    do i_stage = 1,NSTAGE_TIME_SCHEME

      ! loop on all the time steps
      do it = 1,NSTEP

        ! compute current time
        select case(time_stepping_scheme)
        case (1)
          ! Newmark
          timeval = dble(it-1)*DT
        case (2)
          ! LDDRK: Low-Dissipation and low-dispersion Runge-Kutta
          ! note: the LDDRK scheme updates displacement after the stiffness computations and
          !       after adding boundary/coupling/source terms.
          !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme.
          !       we therefore at an additional -DT to have the corresponding timing for the source.
          timeval = dble(it-1-1)*DT + dble(C_LDDRK(i_stage))*DT
        case (3)
          ! RK: Runge-Kutta
          ! note: similar like LDDRK above, displ(n+1) will be determined after stiffness/source/.. computations.
          !       thus, adding an additional -DT to have the same timing in seismogram as Newmark
          timeval = dble(it-1-1)*DT + dble(C_RK4(i_stage))*DT
        case (4)
          ! symplectic PEFRL
          ! note: similar like LDDRK above, displ(n+1) will be determined after final stage of stiffness/source/.. computations.
          !       thus, adding an additional -DT to have the same timing in seismogram as Newmark
          !
          !       for symplectic schemes, the current stage time step size is the sum of all previous and current coefficients
          !          sum( ALPHA_SYMPLECTIC(1:i_stage) ) * DT
          timeval = dble(it-1-1)*DT + dble(sum(ALPHA_SYMPLECTIC(1:i_stage))) * DT
        case default
          call exit_MPI(myrank,'Error invalid time stepping scheme chosen, please check...')
        end select

        ! adjust current source time source time shift
        timeval = timeval - tshift_src(isource)

        ! adjust current source time by simulation start time t0
        t_used = timeval - t0

        ! only process/partition containing source must set STF
        if (myrank == islice_selected_source(isource) .or. SOURCE_IS_MOVING) then

          ! source frequency
          f0 = f0_source(isource)

          ! note regarding USE_TRICK_FOR_BETTER_PRESSURE:
          !   this uses a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
          !   use the second derivative of the source for the source time function instead of the source itself,
          !   and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
          !   this is mathematically equivalent, but numerically significantly more accurate because in the explicit
          !   Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
          !   thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
          !   is accurate at second order and thus contains significantly less numerical noise.

          ! determines source_time_function value for different source types
          select case (time_function_type(isource))
          case (0)
            ! normalized Gaussian
            ! converts frequency to half-duration
            hdur = 1.d0 / f0
            ! convert the half duration for triangle STF to the one for Gaussian STF
            hdur_gauss = hdur / SOURCE_DECAY_MIMIC_TRIANGLE
            if (USE_TRICK_FOR_BETTER_PRESSURE) then
              stf = comp_source_time_function_d2Gaussian_norm(t_used,hdur_gauss)
            else
              stf = comp_source_time_function_Gaussian_norm(t_used,hdur_gauss)
            endif

          case (1)
            ! Ricker
            if (USE_TRICK_FOR_BETTER_PRESSURE) then
              ! Second derivative of Ricker source time function
              stf = comp_source_time_function_d2Ricker(t_used,f0)
            else
              ! Ricker (second derivative of a Gaussian) source time function
              stf = comp_source_time_function_Ricker(t_used,f0)
            endif

          case (2)
            ! first derivative of a Gaussian
            if (USE_TRICK_FOR_BETTER_PRESSURE) then
              ! Third derivative of Gaussian source time function
              stf = comp_source_time_function_d3Gaussian(t_used,f0)
            else
              ! First derivative of a Gaussian source time function
              stf = comp_source_time_function_dGaussian(t_used,f0)
            endif

          case (3,4)
            ! Gaussian/Dirac type
            if (USE_TRICK_FOR_BETTER_PRESSURE) then
              ! Second derivative of Gaussian
              stf = comp_source_time_function_d2Gaussian(t_used,f0)
            else
              ! Gaussian or Dirac (we use a very thin Gaussian instead) source time function
              stf = comp_source_time_function_Gaussian(t_used,f0)
            endif

          case (5)
            ! Heaviside source time function (we use a very thin error function instead)
            f0_sampling = 1.d0 / (10.d0 * DT)    ! empirical sampling frequency limit
            ! if Heaviside source time function, use a very thin error function instead
            if (f0 < f0_sampling) then
              ! can have a smooth quasi-Heaviside if 0 < f0 < sampling-frequency
              ! uses smooth Heaviside with half-duration hdur = 1/f0
              hdur = 1.d0 / f0   ! half-duration
            else
              ! limit frequency to sampling size frequency
              ! converts to half-duration
              hdur = 1.d0 / f0
              ! adds a factor 5/3 to half-duration
              hdur = hdur * 5.d0 / 3.d0
            endif
            ! convert the half duration for triangle STF to the one for Gaussian STF
            hdur_gauss = hdur / SOURCE_DECAY_MIMIC_TRIANGLE
            ! quasi-Heaviside
            stf = comp_source_time_function_heaviside_hdur(t_used,hdur_gauss)

          case (6)
            ! ocean acoustics type I
            stf = comp_source_time_function_ocean_I(t_used,f0)

          case (7)
            ! ocean acoustics type II
            if (USE_TRICK_FOR_BETTER_PRESSURE) then
              ! Second derivative of source 7
              stf = comp_source_time_function_d2ocean_II(t_used,f0)
            else
              ! ocean acoustics type II
              stf = comp_source_time_function_ocean_II(t_used,f0)
            endif

          case (8)
            ! external type
            stf = comp_source_time_function_ext(it,isource)

          case (9)
            ! burst type
            if (USE_TRICK_FOR_BETTER_PRESSURE) then
              ! Second derivative of Burst
              stf = comp_source_time_function_d2burst(t_used,f0,burst_band_width(isource))
            else
              ! Burst source time function
              stf = comp_source_time_function_burst(t_used,f0,burst_band_width(isource))
            endif

          case (10)
            ! Monochromatic/Sinus source time function
            if (USE_TRICK_FOR_BETTER_PRESSURE) then
              ! second derivative of Monochromatic/Sinus
              stf = comp_source_time_function_d2mono(t_used,f0)
            else
              ! Monochromatic/Sinus
              stf = comp_source_time_function_mono(t_used,f0)
            endif

          case (11)
            ! Marmousi Ormsby wavelet
            stf = comp_source_time_function_marmousi(t_used,f0)

          case (12)
            ! Brune source time function
            ! Frequency parameter: hdur == 1/f0 is the source duration or the rise time
            stf = comp_source_time_function_Brune(t_used,f0)

          case (13)
            ! Smoothed Brune source time function
            ! Frequency parameter: hdur == 1/f0 is the source duration or the rise time
            stf = comp_source_time_function_Smooth_Brune(t_used,f0)

          case (14)
            ! Regularized Yoffe
            ! Frequency parameter        -> T_acc == 1/f0  (acceleration duration)
            ! burst band width parameter -> T_eff == 1/bbw (effective duration)
            stf = comp_source_time_function_Yoffe(t_used,f0,burst_band_width(isource))

          case (15)
            ! Integrated Regularized Yoffe
            ! Frequency parameter        -> T_acc == 1/f0  (acceleration duration)
            ! burst band width parameter -> T_eff == 1/bbw (effective duration)
            stf = comp_source_time_function_Yoffe_integrated(t_used,f0,burst_band_width(isource))

          case default
            call exit_MPI(myrank,'unknown source time function')

          end select

          ! stores source time function values
          source_time_function(isource,it,i_stage) = real(stf, kind=CUSTOM_REAL)  ! converts to custom real (working precision)

        endif
      enddo
    enddo
  enddo

  ! prints source time function to file
  if (PRINT_SOURCE_TIME_FUNCTION) call print_stf_file()

  ! source amplification
  ! amplifies source time function (STF) by a factor
  ! note: the amplification factor will amplify the source time function values.
  !       in case this is not desired, one just needs to set the amplification factor to 1.0 in DATA/SOURCE:
  !         ..
  !         ## Amplification (factor to amplify source time function)
  !         factor                          = 1.d0          # amplification factor
  !         ..
  ! loop over all the sources
  do isource = 1,NSOURCES
    ! AXISYM - The following lines could be needed to set absolute amplitudes
    !   use specfem_par, only: rho_vpstore,rhostore,ispec_selected_source
    !   double precision :: rho, cp
    !   logical :: already_done = .false. need to be introduced
    !   ..
    !    if (myrank == islice_selected_source(isource)) then
    !      if (AXISYM) then
    !        if (.not. already_done) then
    !          cp = rho_vpstore(0,0,ispec_selected_source(isource)) / rhostore(0,0,ispec_selected_source(isource))
    !          TODO (above): We must interpolate to find the exact cp value at source location
    !
    !          factor(isource) = - factor(isource)*2.0d0*cp**2*0.45d-5 !0.225d-5
    !          if (time_function_type (isource) == 7)  factor(isource) = factor(isource) * 222066.1d0 !444132.2d0
    !          already_done = .true.
    !        endif
    !      endif
    !    endif

    ! applies source amplification factor
    source_time_function(isource,:,:) = factor(isource) * source_time_function(isource,:,:)
  enddo

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_source_time_function

!
!------------------------------------------------------------
!

  subroutine print_stf_file()

  use constants, only: CUSTOM_REAL,IMAIN,OUTPUT_FILES,MAX_STRING_LEN, &
    C_LDDRK,C_RK4,ALPHA_SYMPLECTIC

  use specfem_par, only: NSTEP,NSOURCES,source_time_function,time_function_type,factor, &
    tshift_src,t0,DT,time_stepping_scheme, &
    myrank,islice_selected_source,NOISE_TOMOGRAPHY

  implicit none

  ! local parameters
  double precision :: timeval
  real(kind=CUSTOM_REAL) :: t_used, stf_used, stf_used_amplified
  integer :: it,ier
  integer :: isource,i_stage
  character(len=MAX_STRING_LEN) :: plot_file

  ! only plot for non-noise simulations
  if (NOISE_TOMOGRAPHY > 0) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  saving the source time function in a text file...'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  do isource = 1,NSOURCES
    ! outputs first stage only for source time functions
    i_stage = 1

    if (myrank == islice_selected_source(isource)) then

      ! opens source time function file
      if (NSOURCES == 1) then
        plot_file = 'plot_source_time_function.txt'
      else
        write(plot_file,"('plot_source_time_function',i0,'.txt')") isource
      endif

      ! opens source time file for output
      open(unit=55,file=trim(OUTPUT_FILES)//trim(plot_file),status='unknown',iostat=ier)
      if (ier /= 0) call stop_the_code('Error opening source time function text-file')

      ! header
      write(55,'("# source time function")')
      write(55,'("# source              : ",i0)') isource
      write(55,'("# time function type  : ",i0)') time_function_type(isource)
      write(55,'("# time shift          : ",es12.5)') tshift_src(isource)
      write(55,'("# amplification factor: ",es12.5)') factor(isource)
      write(55,'("# DT                  : ",es12.5)') DT
      write(55,'("# t0                  : ",es12.5)') t0
      write(55,'("# format")')
      write(55,'("#time  #STF  #STF_amplified")')

      ! source function time values
      do it = 1,NSTEP

        ! compute current time
        select case(time_stepping_scheme)
        case (1)
          ! Newmark
          timeval = dble(it-1)*DT
        case (2)
          ! LDDRK: Low-Dissipation and low-dispersion Runge-Kutta
          ! note: the LDDRK scheme updates displacement after the stiffness computations and
          !       after adding boundary/coupling/source terms.
          !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme.
          !       we therefore at an additional -DT to have the corresponding timing for the source.
          timeval = dble(it-1-1)*DT + dble(C_LDDRK(i_stage))*DT
        case (3)
          ! RK: Runge-Kutta
          ! note: similar like LDDRK above, displ(n+1) will be determined after stiffness/source/.. computations.
          !       thus, adding an additional -DT to have the same timing in seismogram as Newmark
          timeval = dble(it-1-1)*DT + dble(C_RK4(i_stage))*DT
        case (4)
          ! symplectic PEFRL
          ! note: similar like LDDRK above, displ(n+1) will be determined after final stage of stiffness/source/.. computations.
          !       thus, adding an additional -DT to have the same timing in seismogram as Newmark
          !
          !       for symplectic schemes, the current stage time step size is the sum of all previous and current coefficients
          !          sum( ALPHA_SYMPLECTIC(1:i_stage) ) * DT
          timeval = dble(it-1-1)*DT + dble(sum(ALPHA_SYMPLECTIC(1:i_stage))) * DT
        case default
          call exit_MPI(myrank,'Error invalid time stepping scheme chosen, please check...')
        end select

        ! note: earliest start time of the simulation is: (it-1)*DT - t0 - tshift_src(isource)
        t_used = real(timeval - t0 - tshift_src(isource),kind=CUSTOM_REAL)

        ! initial STF
        stf_used = source_time_function(isource,it,i_stage)

        ! amplified STF
        stf_used_amplified = real(stf_used * factor(isource),kind=CUSTOM_REAL)

        ! we'll output both, the initial value given by the source time function (comp_source_time_function**)
        ! and the scaled value that uses the amplification factor as given in the SOURCE file.
        !
        ! note: we call this print-function before the array source_time_function(..) has been amplified.

        ! format: #time #initial_source_function_value #used_STF_amplified
        write(55,*) t_used, stf_used, stf_used_amplified

      enddo

      ! closes STF file
      close(55)

    endif  ! myrank == islice_selected_source(isource)
  enddo

  end subroutine print_stf_file
