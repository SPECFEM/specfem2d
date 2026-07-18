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


  subroutine prepare_attenuation()

  use constants, only: IMAIN,TWO,PI,FOUR_THIRDS,TWO_THIRDS,USE_A_STRONG_FORMULATION_FOR_E1, &
                       USE_OLD_C_ATTENUATION_ROUTINE_INSTEAD,USE_FIXED_ATTENUATION_ABSORPTION_BAND, &
                       TINYVAL
  use specfem_par

  implicit none

  ! local parameters
  integer :: i,j,ispec,ier

  ! for shifting of velocities if needed in the case of viscoelasticity
  double precision :: vp,vs,rhol,mul,kappal
  double precision :: qkappal,qmul

  ! temporary attenuation array factors
  real(kind=CUSTOM_REAL) :: Mu_nu1_sent,Mu_nu2_sent
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: tau_epsilon_nu1_sent,tau_epsilon_nu2_sent
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: inv_tau_sigma_nu1_sent,inv_tau_sigma_nu2_sent, &
                                                       phi_nu1_sent,phi_nu2_sent
  real(kind=CUSTOM_REAL), dimension(N_SLS) ::  phinu,tauinvnu,temp,coef

  ! EM attenuation
  double precision :: tauinv(2)

  ! info
  logical :: show_warning
  integer :: nelem_viscoacoustic,iproc

  ! absorption band
  double precision :: f_min_attenuation, f_max_attenuation, f_center
  double precision :: T_min_period, T_max_period, min_resolved_period

  ! attenuation
  ! user output
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) 'Attenuation:'
    write(IMAIN,*) '  viscoelastic  attenuation: ',ATTENUATION_VISCOELASTIC, &
                   '(shear & bulk attenuation in elastic domains)'
    write(IMAIN,*) '  viscoacoustic attenuation: ',ATTENUATION_VISCOACOUSTIC, &
                   '(bulk attenuation in acoustic domains)'
    write(IMAIN,*) '  poro-fluid    attenuation: ',ATTENUATION_PORO_FLUID_PART, &
                   '(viscous attenuation in poroelastic domains)'
    write(IMAIN,*) '  permittivity  attenuation: ',ATTENUATION_PERMITTIVITY, &
                   '(for electromagnetic materials)'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! attenuation array allocations
  call prepare_attenuation_global_arrays()

  ! checks if anything further to do
  if (.not. ATTENUATION_VISCOELASTIC .and. &
      .not. ATTENUATION_VISCOACOUSTIC .and. &
      .not. ATTENUATION_PORO_FLUID_PART .and. &
      .not. ATTENUATION_PERMITTIVITY) return

  ! physical dispersion: scales moduli from reference frequency to simulation (source) center frequency
  !
  ! if attenuation is on, shift the velocity model to right frequency;
  ! rescale mu to average frequency for attenuation
  !
  ! the formulas to implement the scaling can be found for instance in
  ! Liu, H. P., Anderson, D. L. and Kanamori, H., Velocity dispersion due to
  ! anelasticity: implications for seismology and mantle composition,
  ! Geophys. J. R. Astron. Soc., vol. 47, pp. 41-58 (1976)
  !
  ! and in Aki, K. and Richards, P. G., Quantitative seismology, theory and methods,
  ! W. H. Freeman, (1980), second edition, sections 5.5 and 5.5.2, eq. (5.81) p. 170.
  !
  ! Beware that in the book of Aki and Richards eq. (5.81) is given for velocities
  ! while we need an equation for "mu" and thus we have an additional factor of 2
  ! in the scaling factor below and in equation (49) of Komatitsch and Tromp, Geophys. J. Int. (2002) 149, 390-412,
  ! because "mu" is related to the square of velocity.
  !
  ! mu(omega_c) = mu(omega_0)[ 1 + 2/(pi Q_mu) ln(omega_c / omega_0) ]
  !
  ! determines reference frequency
  if (READ_VELOCITIES_AT_f0) then
    ! velocity model given at user-specified reference frequency ATTENUATION_f0_REFERENCE
    ! ATTENUATION_f0_REFERENCE given already
    continue
  else
    ! shift reference frequency to dominant frequency of (first) source;
    ! if source is not a Dirac or Heavyside then ATTENUATION_f0_REFERENCE is f0 of the first source
    if (.not. (time_function_type(1) == 4 .or. time_function_type(1) == 5)) then
      ATTENUATION_f0_REFERENCE = f0_source(1)
      if (ATTENUATION_PERMITTIVITY) f0_electromagnetic = f0_source(1)
    endif
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "Preparing attenuation"
    write(IMAIN,*) "  The code uses a constant Q quality factor, but approximated"
    write(IMAIN,*) "  based on a series of Zener standard linear solids (SLS)."
    write(IMAIN,*)
    write(IMAIN,*) "  number of SLS bodies: ",N_SLS
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! setup attenuation
  if (ATTENUATION_VISCOELASTIC .or. ATTENUATION_VISCOACOUSTIC) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  Attenuation in viscoelastic or viscoacoustic parts of the model:'
      if (READ_VELOCITIES_AT_f0) then
        write(IMAIN,*) '  reading velocity at f0 : ',READ_VELOCITIES_AT_f0
        write(IMAIN,*) '  assuming velocity model given at reference frequency'
        write(IMAIN,*)
        write(IMAIN,*) "  Reference frequency set by the user (Hz):",sngl(ATTENUATION_f0_REFERENCE), &
                                                      " period (s):",sngl(1.0/ATTENUATION_f0_REFERENCE)
      else
        write(IMAIN,*) '  assuming velocity model given at unrelaxed state (infinite frequency)'
        write(IMAIN,*)
        if (.not. (time_function_type(1) == 4 .or. time_function_type(1) == 5)) then
          write(IMAIN,*) "  Reference frequency set by (first) source frequency (Hz):",sngl(ATTENUATION_f0_REFERENCE), &
                                                                        " period (s):",sngl(1.0/ATTENUATION_f0_REFERENCE)
        else
          write(IMAIN,*) "  Reference frequency set by the user (Hz):",sngl(ATTENUATION_f0_REFERENCE), &
                                                        " period (s):",sngl(1.0/ATTENUATION_f0_REFERENCE)
        endif
      endif
      write(IMAIN,*)
      write(IMAIN,*) "  Approximation is performed in the following frequency band:"
      if (USE_FIXED_ATTENUATION_ABSORPTION_BAND) then
        write(IMAIN,*) "  Using fixed frequency band selection"
      else
        if (COMPUTE_FREQ_BAND_AUTOMATIC) then
          write(IMAIN,*) "  The following values are computed automatically by the code"
          write(IMAIN,*) "  based on the estimated maximum frequency resolution of your mesh"
          write(IMAIN,*) "  and can thus vary from what you have requested."
        else
          write(IMAIN,*) "  Using user specified min/max attenuation periods"
        endif
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif
    call synchronize_all()

    ! safety check
    if (N_SLS < 1) call stop_the_code('Invalid N_SLS value, must be at least 1')

    ! attenuation constants for standard linear solids
    ! nu1 is the dilatation/incompressibility mode (QKappa)
    ! nu2 is the shear mode (Qmu)
    ! array index (1) is the first standard linear solid, (2) is the second etc.

    ! initializes frequency band
    f_min_attenuation = 1.d0
    f_max_attenuation = 1.d0

    ! determine minimum/maximum frequency of attenuation band
    if (USE_FIXED_ATTENUATION_ABSORPTION_BAND) then
      ! old way - uses fixed decades
      ! note: in both cases, fmin/max are chosen such that f0 is log-centered.
      !       the spread between fmin and fmax in the old C-attenuation routine was ~1.08 decades
      !       (given fmax = 12*fmin -> band width was log10(12)=1.08 decades),
      !       in the newer case, it was 2 decades
      !       (given fmax=f0*10,fmin=f0/10 and fmax = 10**2 fmin -> band width was log(10**2)=2 decades)
      if (USE_OLD_C_ATTENUATION_ROUTINE_INSTEAD .and. .not. USE_SOLVOPT) then
        ! f_min and f_max are computed as : f_max/f_min=12 and (log(f_min)+log(f_max))/2 = log(f0)
        f_min_attenuation = exp(log(ATTENUATION_f0_REFERENCE) - log(12.d0)/2.d0)
        f_max_attenuation = 12.d0 * f_min_attenuation
        f_center = ATTENUATION_f0_REFERENCE
      else
        ! use a wide bandwidth (always OK when using three or more Standard Linear Solids, can be a bit inaccurate if using only two)
        f_min_attenuation = ATTENUATION_f0_REFERENCE / 10.d0
        f_max_attenuation = ATTENUATION_f0_REFERENCE * 10.d0
        f_center = ATTENUATION_f0_REFERENCE
      endif
    else
      ! new way - compatible w/ SPECFEM3D Cartesian
      ! note: absorption band gets determined by number of SLS and minimum period the mesh resolves.
      !       the minimum period is estimated in the check_grid() routine. a standard table of the band width (in decades)
      !       is then used to determine min/max frequencies of the absorption band.
      if (COMPUTE_FREQ_BAND_AUTOMATIC) then ! otherwise they were entered as input values by the user in the Par_file
        ! mesh resolved minimum period
        min_resolved_period = mesh_T_min
        ! uses a small margin
        min_resolved_period = 0.99 * min_resolved_period

        ! determines absorption band (fmin/fmax) by number of SLS and minimum period
        call get_attenuation_periods(N_SLS,min_resolved_period,T_min_period,T_max_period)

        ! frequency band
        f_min_attenuation = 1.d0 / T_max_period  ! minimum frequency
        f_max_attenuation = 1.d0 / T_min_period  ! maximum frequency
      else
        ! user specified min/max attenuation periods in Par_file
        ! check user values
        if (MIN_ATTENUATION_PERIOD < TINYVAL) &
          call stop_the_code("Invalid MIN_ATTENUATION_PERIOD in Par_file, must be > 0.0")
        if (MAX_ATTENUATION_PERIOD < TINYVAL) &
          call stop_the_code("Invalid MAX_ATTENUATION_PERIOD in Par_file, must be > 0.0")

        ! frequency band
        f_min_attenuation = 1.d0 / MAX_ATTENUATION_PERIOD
        f_max_attenuation = 1.d0 / MIN_ATTENUATION_PERIOD
      endif

      ! determines center frequency of absorption band
      call get_attenuation_center_freq(f_center,f_min_attenuation,f_max_attenuation)
    endif

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  Attenuation frequency band min/max (Hz):",sngl(f_min_attenuation), &
                                                              '/',sngl(f_max_attenuation)
      write(IMAIN,*) "              period band    min/max (s) :",sngl(1.0/f_max_attenuation), &
                                                              '/',sngl(1.0/f_min_attenuation)
      write(IMAIN,*) "  Logarithmic center frequency (Hz):",sngl(f_center)
      write(IMAIN,*) "                     period     (s):",sngl(1.0/f_center)
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! temporary arrays for function argument
    allocate(tau_epsilon_nu1_sent(N_SLS), &
             tau_epsilon_nu2_sent(N_SLS), &
             inv_tau_sigma_nu1_sent(N_SLS), &
             inv_tau_sigma_nu2_sent(N_SLS), &
             phi_nu1_sent(N_SLS), &
             phi_nu2_sent(N_SLS),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating attenuation coefficient arrays')
    tau_epsilon_nu1_sent(:) = 0.0_CUSTOM_REAL
    tau_epsilon_nu2_sent(:) = 0.0_CUSTOM_REAL
    inv_tau_sigma_nu1_sent(:) = 0.0_CUSTOM_REAL
    inv_tau_sigma_nu2_sent(:) = 0.0_CUSTOM_REAL
    phi_nu1_sent(:) = 0.0_CUSTOM_REAL
    phi_nu2_sent(:) = 0.0_CUSTOM_REAL

    show_warning = .false.
    nelem_viscoacoustic = 0

    ! define the attenuation quality factors.
    do ispec = 1,nspec

      ! checks element type for velocity shifts
      if (READ_VELOCITIES_AT_f0) then
        ! safety check
        if (ispec_is_anisotropic(ispec) .or. ispec_is_poroelastic(ispec)) &
          call stop_the_code('READ_VELOCITIES_AT_f0 only implemented for non anisotropic, non-poroelastic materials for now')

        if (ispec_is_acoustic(ispec)) then
          show_warning = .true.
          nelem_viscoacoustic = nelem_viscoacoustic + 1
        endif
      endif

      do j = 1,NGLLZ
        do i = 1,NGLLX
          ! determines relaxation factors
          ! bulk attenuation
          qkappal = qkappa_attenuation_store(i,j,ispec)
          ! shear attenuation
          qmul = qmu_attenuation_store(i,j,ispec)

          ! if no attenuation in that elastic element
          if (qkappal > 9998.999d0 .and. qmul > 9998.999d0) cycle

          ! determines attenuation factors
          call attenuation_model(qkappal,qmul,ATTENUATION_f0_REFERENCE,f_min_attenuation,f_max_attenuation, &
                                 N_SLS, &
                                 tau_epsilon_nu1_sent,inv_tau_sigma_nu1_sent,phi_nu1_sent,Mu_nu1_sent, &
                                 tau_epsilon_nu2_sent,inv_tau_sigma_nu2_sent,phi_nu2_sent,Mu_nu2_sent)

          ! stores attenuation values
          ! bulk attenuation (Qkappa)
          inv_tau_sigma_nu1(i,j,ispec,:) = inv_tau_sigma_nu1_sent(:)
          tau_epsilon_nu1(i,j,ispec,:) = tau_epsilon_nu1_sent(:)
          phi_nu1(i,j,ispec,:) = phi_nu1_sent(:)

          ! shear attenuation (Qmu)
          inv_tau_sigma_nu2(i,j,ispec,:) = inv_tau_sigma_nu2_sent(:)
          tau_epsilon_nu2(i,j,ispec,:) = tau_epsilon_nu2_sent(:)
          phi_nu2(i,j,ispec,:) = phi_nu2_sent(:)

          Mu_nu1(i,j,ispec) = Mu_nu1_sent
          Mu_nu2(i,j,ispec) = Mu_nu2_sent

          ! acoustic domains
          if (ATTENUATION_VISCOACOUSTIC .and. USE_A_STRONG_FORMULATION_FOR_E1 .and. time_stepping_scheme == 1 ) then
            ! bulk attenuation (Qkappa)
            phinu(:)    = phi_nu1(i,j,ispec,:)
            tauinvnu(:) = inv_tau_sigma_nu1(i,j,ispec,:)
            temp(:)      = exp(- 0.5d0 * tauinvnu(:) * DT)
            coef(:)     = (1.d0 - temp(:)) / tauinvnu(:)
            A_newmark_e1_sf(:,i,j,ispec) = temp(:)
            B_newmark_e1_sf(:,i,j,ispec) = phinu(:) * coef(:)
          endif

          ! elastic domains
          if (ATTENUATION_VISCOELASTIC .and. time_stepping_scheme == 1 ) then
            ! bulk attenuation (Qkappa)
            phinu(:)    = phi_nu1(i,j,ispec,:)
            tauinvnu(:) = inv_tau_sigma_nu1(i,j,ispec,:)
            temp(:)      = exp(- 0.5d0 * tauinvnu(:) * DT)
            coef(:)     = (1.d0 - temp(:)) / tauinvnu(:)
            A_newmark_nu1(:,i,j,ispec) = temp(:)
            B_newmark_nu1(:,i,j,ispec) = phinu(:) * coef(:)

            ! shear attenuation (Qmu)
            phinu(:)    = phi_nu2(i,j,ispec,:)
            tauinvnu(:) = inv_tau_sigma_nu2(i,j,ispec,:)
            temp(:)      = exp(- 0.5d0 * tauinvnu(:) * DT)
            coef(:)     = (1.d0 - temp(:)) / tauinvnu(:)
            A_newmark_nu2(:,i,j,ispec) = temp(:)
            B_newmark_nu2(:,i,j,ispec) = phinu(:) * coef(:)
          endif

          ! shifts velocities
          if (READ_VELOCITIES_AT_f0) then
            ! shifts velocity model
            rhol = dble(rhostore(i,j,ispec))
            vp = dble(rho_vpstore(i,j,ispec)/rhol)
            vs = dble(rho_vsstore(i,j,ispec)/rhol)

            ! shifts vp and vs (according to f0 and attenuation band)
            call shift_velocities_from_f0(vp,vs,rhol,qkappal,qmul, &
                                          ATTENUATION_f0_REFERENCE,f_center, &
                                          N_SLS,tau_epsilon_nu1_sent,tau_epsilon_nu2_sent, &
                                          inv_tau_sigma_nu1_sent,inv_tau_sigma_nu2_sent)

            ! stores shifted values
            ! determines mu and kappa
            mul = rhol * vs * vs
            if (AXISYM) then ! CHECK kappa
              kappal = rhol * vp * vp - FOUR_THIRDS * mul
            else
              kappal = rhol * vp * vp - mul
            endif
            ! to compare:
            !lambdal = rhol * vp*vp - TWO * mul
            !if (AXISYM) then ! CHECK kappa
            !  kappal = lambdal + TWO_THIRDS * mul
            !  vp = sqrt((kappal + FOUR_THIRDS * mul)/rhol)
            !else
            !  kappal = lambdal + mul
            !  vp = sqrt((kappal + mul)/rhol)
            !endif

            ! stores unrelaxed moduli
            mustore(i,j,ispec) = mul
            kappastore(i,j,ispec) = kappal

            ! stores density times vp and vs
            rho_vpstore(i,j,ispec) = rhol * vp
            rho_vsstore(i,j,ispec) = rhol * vs
          endif
        enddo
      enddo
    enddo

    ! warning
    do iproc = 0,NPROC
      if (iproc == myrank) then
        if (show_warning) then
          if (myrank == 0) then
            ! warning text
            print *
            print *,'******************************'
            print *,'WARNING: READ_VELOCITIES_AT_f0 in viscoacoustic elements may imply having to rebuild the mass matrix &
                 &with the shifted velocities, since the fluid mass matrix contains Kappa; not implemented yet, BEWARE!!'
            print *,'******************************'
            print *
          endif
          print *,'  slice ',myrank,' has ',nelem_viscoacoustic,' viscoacoustic elements'
        endif
      endif
      call synchronize_all()
    enddo

    ! for PMLs
    if (PML_BOUNDARY_CONDITIONS) call prepare_attenuation_with_PML()

  endif ! of if (ATTENUATION_VISCOELASTIC .or. ATTENUATION_VISCOACOUSTIC)

  ! allocate memory variables for viscous attenuation (poroelastic media)
  if (ATTENUATION_PORO_FLUID_PART) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  Attenuation in poroelastic parts of the model:'
      write(IMAIN,*) '  f0 poroelastic                 : ',freq0_poroelastic
      write(IMAIN,*) '  Q0 poroelastic                 : ',Q0_poroelastic
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! precompute Runge Kutta coefficients if viscous attenuation
    ! viscous attenuation is implemented following the memory variable formulation of
    ! J. M. Carcione Wave fields in real media: wave propagation in anisotropic,
    ! anelastic and porous media, Elsevier, p. 304-305, 2007
    theta_e = (sqrt(Q0_poroelastic**2+1.d0) +1.d0)/(2.d0*pi*freq0_poroelastic*Q0_poroelastic)
    theta_s = (sqrt(Q0_poroelastic**2+1.d0) -1.d0)/(2.d0*pi*freq0_poroelastic*Q0_poroelastic)

    thetainv = - 1.d0 / theta_s
    alphaval = 1.d0 + DT * thetainv + DT**2 * thetainv**2 / 2.d0 &
                    + DT**3 * thetainv**3 / 6.d0 + DT**4 * thetainv**4 / 24.d0
    betaval = DT / 2.d0 + DT**2 * thetainv / 3.d0 + DT**3 * thetainv**2 / 8.d0 + DT**4 * thetainv**3 / 24.d0
    gammaval = DT / 2.d0 + DT**2 * thetainv / 6.d0 + DT**3 * thetainv**2 / 24.d0
  endif

  ! allocate memory variables for permittivity attenuation (electromagnetic media)
  if (ATTENUATION_PERMITTIVITY) then
    ! if source is not a Dirac or Heavyside then f0_attenuation is f0 of the first
    ! source
    if (.not. (time_function_type(1) == 4 .or. time_function_type(1) == 5)) then
      f0_electromagnetic = f0_source(1)
    endif

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  Attenuation permittivity in electromagnetic parts of the model:'
      write(IMAIN,*) '  f0 electromagnetic                 : ',f0_electromagnetic
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! precompute Runge Kutta coefficients if viscous attenuation
    ! viscous attenuation is implemented following the memory variable
    ! formulation of
    ! J. M. Carcione Wave fields in real media: wave propagation in anisotropic,
    ! anelastic and porous media, Elsevier, p. 304-305, 2007
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX

          tau_e(i,j,ispec,1) = (sqrt(Qe11_electromagnetic(kmato(ispec))**2+1.d0) +1.d0)/ &
                                (2.d0*pi*f0_electromagnetic*Qe11_electromagnetic(kmato(ispec)))
          tau_e(i,j,ispec,2) = (sqrt(Qe33_electromagnetic(kmato(ispec))**2+1.d0) +1.d0)/ &
                                (2.d0*pi*f0_electromagnetic*Qe33_electromagnetic(kmato(ispec)))
          tau_d(i,j,ispec,1) = (sqrt(Qe11_electromagnetic(kmato(ispec))**2+1.d0) -1.d0)/ &
                                (2.d0*pi*f0_electromagnetic*Qe11_electromagnetic(kmato(ispec)))
          tau_d(i,j,ispec,2) = (sqrt(Qe33_electromagnetic(kmato(ispec))**2+1.d0) -1.d0)/ &
                               (2.d0*pi*f0_electromagnetic*Qe33_electromagnetic(kmato(ispec)))

          tauinv(:) = - 1.d0 / tau_d(i,j,ispec,:)
          alphaval_em(i,j,ispec,:) = 1.d0 + deltat*tauinv(:) + deltat**2*tauinv(:)**2 / 2.d0 &
                          + deltat**3*tauinv(:)**3 / 6.d0 + deltat**4*tauinv(:)**4 / 24.d0
          betaval_em(i,j,ispec,:) = deltat / 2.d0 + deltat**2*tauinv(:) / 3.d0 + deltat**3*tauinv(:)**2 / 8.d0 &
                          + deltat**4*tauinv(:)**3 / 24.d0
          gammaval_em(i,j,ispec,:) = deltat / 2.d0 + deltat**2*tauinv(:) / 6.d0 + deltat**3*tauinv(:)**2 / 24.d0

        enddo
      enddo
    enddo
  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_attenuation

!
!-------------------------------------------------------------------------------------
!

  subroutine prepare_attenuation_global_arrays()

! allocates all attenuation arrays (required with and without attenuation)

  use constants, only: IMAIN,TWO,PI,FOUR_THIRDS,TWO_THIRDS,USE_A_STRONG_FORMULATION_FOR_E1
  use specfem_par

  implicit none

  ! local parameters
  integer :: ier

  ! elastic domains
  if (ATTENUATION_VISCOELASTIC) then
    nspec_ATT_el = nspec
  else
    nspec_ATT_el = 1
  endif

  ! allocate memory variables for attenuation
  allocate(e1(N_SLS,NGLLX,NGLLZ,nspec_ATT_el), &
           e11(N_SLS,NGLLX,NGLLZ,nspec_ATT_el), &
           e13(N_SLS,NGLLX,NGLLZ,nspec_ATT_el), &
           dux_dxl_old(NGLLX,NGLLZ,nspec_ATT_el), &
           duz_dzl_old(NGLLX,NGLLZ,nspec_ATT_el), &
           dux_dzl_plus_duz_dxl_old(NGLLX,NGLLZ,nspec_ATT_el), &
           A_newmark_nu1(N_SLS,NGLLX,NGLLZ,nspec_ATT_el), &
           B_newmark_nu1(N_SLS,NGLLX,NGLLZ,nspec_ATT_el), &
           A_newmark_nu2(N_SLS,NGLLX,NGLLZ,nspec_ATT_el), &
           B_newmark_nu2(N_SLS,NGLLX,NGLLZ,nspec_ATT_el), stat=ier)

  e1(:,:,:,:) = 0._CUSTOM_REAL
  e11(:,:,:,:) = 0._CUSTOM_REAL
  e13(:,:,:,:) = 0._CUSTOM_REAL
  dux_dxl_old(:,:,:) = 0._CUSTOM_REAL
  duz_dzl_old(:,:,:) = 0._CUSTOM_REAL
  dux_dzl_plus_duz_dxl_old(:,:,:) = 0._CUSTOM_REAL

  ! acoustic domains
  if (ATTENUATION_VISCOACOUSTIC) then
    nglob_ATT = nglob
    nspec_ATT_ac = nspec
  else
    nglob_ATT = 1
    nspec_ATT_ac = 1
  endif

  allocate(e1_acous_sf(N_SLS,NGLLX,NGLLZ,nspec_ATT_ac), &
           sum_forces_old(NGLLX,NGLLZ,nspec_ATT_ac), stat=ier)

  if (ATTENUATION_VISCOACOUSTIC .and. .not. USE_A_STRONG_FORMULATION_FOR_E1) then
    allocate(e1_acous(nglob_acoustic,N_SLS), &
             dot_e1(nglob_acoustic,N_SLS), &
             A_newmark_e1_sf(1,1,1,1), &
             B_newmark_e1_sf(1,1,1,1),stat=ier)
    if (time_stepping_scheme == 1) then
      ! Newmark scheme
      allocate(dot_e1_old(nglob_acoustic,N_SLS), &
               A_newmark_e1(nglob_acoustic,N_SLS), &
               B_newmark_e1(nglob_acoustic,N_SLS),stat=ier)
    else
      ! dummy arrays
      allocate(dot_e1_old(1,N_SLS), &
               A_newmark_e1(1,N_SLS), &
               B_newmark_e1(1,N_SLS),stat=ier)
    endif

    if (time_stepping_scheme == 2) then
      ! LDDRK scheme
      allocate(e1_acous_temp(nglob_acoustic,N_SLS),stat=ier)
    else
      ! dummy array
      allocate(e1_acous_temp(1,N_SLS),stat=ier)
    endif

  else if (ATTENUATION_VISCOACOUSTIC .and. USE_A_STRONG_FORMULATION_FOR_E1) then
    allocate(e1_acous(1,N_SLS), &
             e1_acous_temp(1,N_SLS), &
             dot_e1(1,N_SLS), &
             dot_e1_old(1,N_SLS), &
             A_newmark_e1(1,N_SLS), &
             B_newmark_e1(1,N_SLS),stat=ier)
    allocate(A_newmark_e1_sf(N_SLS,NGLLX,NGLLZ,nspec), &
             B_newmark_e1_sf(N_SLS,NGLLX,NGLLZ,nspec),stat=ier)

  else
    ! no ATTENUATION_VISCOACOUSTIC
    ! dummy arrays
    allocate(e1_acous(1,N_SLS), &
             e1_acous_temp(1,N_SLS), &
             dot_e1(1,N_SLS), &
             dot_e1_old(1,N_SLS), &
             A_newmark_e1(1,N_SLS), &
             B_newmark_e1(1,N_SLS), &
             A_newmark_e1_sf(1,1,1,1), &
             B_newmark_e1_sf(1,1,1,1),stat=ier)
  endif
  if (ier /= 0) call stop_the_code('Error allocating attenuation arrays')

  e1_acous(:,:) = 0._CUSTOM_REAL
  e1_acous_sf(:,:,:,:) = 0._CUSTOM_REAL

  dot_e1_old = 0._CUSTOM_REAL
  dot_e1     = 0._CUSTOM_REAL
  sum_forces_old = 0._CUSTOM_REAL

  ! backward/kernel simulations
  if (SIMULATION_TYPE == 3) then
    ! acoustic domains
    if (any_acoustic) then
      allocate(b_e1_acous_sf(N_SLS,NGLLX,NGLLZ,nspec_ATT_ac), &
               b_sum_forces_old(NGLLX,NGLLZ,nspec_ATT_ac),stat=ier)
      if (ier /= 0) call stop_the_code('Error allocating acoustic attenuation arrays')
      b_e1_acous_sf(:,:,:,:) = 0._CUSTOM_REAL
      b_sum_forces_old(:,:,:) = 0._CUSTOM_REAL
    endif
    ! elastic domains
    if (any_elastic) then
      allocate(b_e1(N_SLS,NGLLX,NGLLZ,nspec_ATT_el), &
               b_e11(N_SLS,NGLLX,NGLLZ,nspec_ATT_el), &
               b_e13(N_SLS,NGLLX,NGLLZ,nspec_ATT_el), &
               b_dux_dxl_old(NGLLX,NGLLZ,nspec_ATT_el), &
               b_duz_dzl_old(NGLLX,NGLLZ,nspec_ATT_el), &
               b_dux_dzl_plus_duz_dxl_old(NGLLX,NGLLZ,nspec_ATT_el),stat=ier)
      if (ier /= 0) call stop_the_code('Error allocating attenuation arrays')
      b_e1(:,:,:,:) = 0._CUSTOM_REAL
      b_e11(:,:,:,:) = 0._CUSTOM_REAL
      b_e13(:,:,:,:) = 0._CUSTOM_REAL
      b_dux_dxl_old(:,:,:) = 0._CUSTOM_REAL
      b_duz_dzl_old(:,:,:) = 0._CUSTOM_REAL
      b_dux_dzl_plus_duz_dxl_old(:,:,:) = 0._CUSTOM_REAL
    endif
  endif

  ! LDDRK time scheme
  if (time_stepping_scheme == 2) then
    ! LDDRK
    ! elastic domains
    if (ATTENUATION_VISCOELASTIC) then
      allocate(e1_LDDRK(NGLLX,NGLLZ,nspec_ATT_el,N_SLS))
      allocate(e11_LDDRK(NGLLX,NGLLZ,nspec_ATT_el,N_SLS))
      allocate(e13_LDDRK(NGLLX,NGLLZ,nspec_ATT_el,N_SLS))
    else
      allocate(e1_LDDRK(1,1,1,1))
      allocate(e11_LDDRK(1,1,1,1))
      allocate(e13_LDDRK(1,1,1,1))
    endif
    ! acoustic domains
    if (ATTENUATION_VISCOACOUSTIC) then
        allocate(e1_LDDRK_acous(nglob_att,N_SLS))
    else
        allocate(e1_LDDRK_acous(1,1))
    endif
  else
    ! dummy arrays
    allocate(e1_LDDRK(1,1,1,1))
    allocate(e11_LDDRK(1,1,1,1))
    allocate(e13_LDDRK(1,1,1,1))
    allocate(e1_LDDRK_acous(1,1))
  endif
  e1_LDDRK(:,:,:,:) = 0._CUSTOM_REAL
  e11_LDDRK(:,:,:,:) = 0._CUSTOM_REAL
  e13_LDDRK(:,:,:,:) = 0._CUSTOM_REAL

  e1_LDDRK_acous(:,:) = 0._CUSTOM_REAL

  ! RK time scheme
  if (time_stepping_scheme == 3) then
    ! RK scheme
    ! elastic domains
    allocate(e1_initial_rk(NGLLX,NGLLZ,nspec_ATT_el,N_SLS))
    allocate(e11_initial_rk(NGLLX,NGLLZ,nspec_ATT_el,N_SLS))
    allocate(e13_initial_rk(NGLLX,NGLLZ,nspec_ATT_el,N_SLS))
    allocate(e1_force_rk(NGLLX,NGLLZ,nspec_ATT_el,N_SLS,NSTAGE_TIME_SCHEME))
    allocate(e11_force_rk(NGLLX,NGLLZ,nspec_ATT_el,N_SLS,NSTAGE_TIME_SCHEME))
    allocate(e13_force_rk(NGLLX,NGLLZ,nspec_ATT_el,N_SLS,NSTAGE_TIME_SCHEME))
    ! acoustic domains
    if (ATTENUATION_VISCOACOUSTIC) then
      allocate(e1_initial_rk_acous(nglob_att,N_SLS))
      allocate(e1_force_rk_acous(nglob_att,N_SLS,NSTAGE_TIME_SCHEME))
    else
      allocate(e1_initial_rk_acous(1,1))
      allocate(e1_force_rk_acous(1,1,1))
    endif
  else
    ! dummy arrays
    allocate(e1_initial_rk(1,1,1,1))
    allocate(e11_initial_rk(1,1,1,1))
    allocate(e13_initial_rk(1,1,1,1))
    allocate(e1_force_rk(1,1,1,1,1))
    allocate(e11_force_rk(1,1,1,1,1))
    allocate(e13_force_rk(1,1,1,1,1))

    allocate(e1_initial_rk_acous(1,1))
    allocate(e1_force_rk_acous(1,1,1))
  endif
  e1_initial_rk(:,:,:,:) = 0._CUSTOM_REAL
  e11_initial_rk(:,:,:,:) = 0._CUSTOM_REAL
  e13_initial_rk(:,:,:,:) = 0._CUSTOM_REAL
  e1_force_rk(:,:,:,:,:) = 0._CUSTOM_REAL
  e11_force_rk(:,:,:,:,:) = 0._CUSTOM_REAL
  e13_force_rk(:,:,:,:,:) = 0._CUSTOM_REAL

  e1_initial_rk_acous(:,:) = 0._CUSTOM_REAL
  e1_force_rk_acous(:,:,:) = 0._CUSTOM_REAL

  ! attenuation arrays
  allocate(inv_tau_sigma_nu1(NGLLX,NGLLZ,max(nspec_ATT_el,nspec_ATT_ac),N_SLS), &
           inv_tau_sigma_nu2(NGLLX,NGLLZ,max(nspec_ATT_el,nspec_ATT_ac),N_SLS), &
           phi_nu1(NGLLX,NGLLZ,max(nspec_ATT_el,nspec_ATT_ac),N_SLS), &
           phi_nu2(NGLLX,NGLLZ,max(nspec_ATT_el,nspec_ATT_ac),N_SLS), &
           Mu_nu1(NGLLX,NGLLZ,max(nspec_ATT_el,nspec_ATT_ac)), &
           Mu_nu2(NGLLX,NGLLZ,max(nspec_ATT_el,nspec_ATT_ac)), &
! ZX ZX needed for further optimization with nspec_ATT_el replaced with nspec_PML
           tau_epsilon_nu1(NGLLX,NGLLZ,max(nspec_ATT_el,nspec_ATT_ac),N_SLS), &
           tau_epsilon_nu2(NGLLX,NGLLZ,max(nspec_ATT_el,nspec_ATT_ac),N_SLS),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating attenuation arrays')

  ! initialize to dummy values
  ! convention to indicate that Q = 9999 in that element i.e. that there is no viscoelasticity in that element
  inv_tau_sigma_nu1(:,:,:,:) = -1._CUSTOM_REAL
  inv_tau_sigma_nu2(:,:,:,:) = -1._CUSTOM_REAL

  tau_epsilon_nu1(:,:,:,:) = -1._CUSTOM_REAL
  tau_epsilon_nu2(:,:,:,:) = -1._CUSTOM_REAL

  phi_nu1(:,:,:,:) = 0._CUSTOM_REAL
  phi_nu2(:,:,:,:) = 0._CUSTOM_REAL

  ! do not change this, in the case of a viscoacoustic medium the mass matrix is multiplied by this,
  ! and thus the factor needs to be equal to +1 when QKappa = 9999 i.e. when viscoacousticity is turned off in parts of the medium
  Mu_nu1(:,:,:) = +1._CUSTOM_REAL
  Mu_nu2(:,:,:) = +1._CUSTOM_REAL

  ! allocate memory variables for viscous attenuation (poroelastic media)
  if (ATTENUATION_PORO_FLUID_PART) then
    allocate(rx_viscous(NGLLX,NGLLZ,nspec))
    allocate(rz_viscous(NGLLX,NGLLZ,nspec))
    allocate(viscox(NGLLX,NGLLZ,nspec))
    allocate(viscoz(NGLLX,NGLLZ,nspec))
    ! initialize memory variables for attenuation
    rx_viscous(:,:,:) = 0.d0
    rz_viscous(:,:,:) = 0.d0
    viscox(:,:,:) = 0.d0
    viscoz(:,:,:) = 0.d0

    if (time_stepping_scheme == 2) then
      allocate(rx_viscous_LDDRK(NGLLX,NGLLZ,nspec))
      allocate(rz_viscous_LDDRK(NGLLX,NGLLZ,nspec))
      rx_viscous_LDDRK(:,:,:) = 0.d0
      rz_viscous_LDDRK(:,:,:) = 0.d0
    endif

    if (time_stepping_scheme == 3) then
      allocate(rx_viscous_initial_rk(NGLLX,NGLLZ,nspec))
      allocate(rz_viscous_initial_rk(NGLLX,NGLLZ,nspec))
      allocate(rx_viscous_force_RK(NGLLX,NGLLZ,nspec,NSTAGE_TIME_SCHEME))
      allocate(rz_viscous_force_RK(NGLLX,NGLLZ,nspec,NSTAGE_TIME_SCHEME))
      rx_viscous_initial_rk(:,:,:) = 0.d0
      rz_viscous_initial_rk(:,:,:) = 0.d0
      rx_viscous_force_RK(:,:,:,:) = 0.d0
      rz_viscous_force_RK(:,:,:,:) = 0.d0
    endif
  endif

  ! allocate memory variables for permittivity attenuation (electromagnetic media)
  if (ATTENUATION_PERMITTIVITY) then
    allocate(rx_permattenuation(NGLLX,NGLLZ,nspec), &
             rz_permattenuation(NGLLX,NGLLZ,nspec), &
             permx(NGLLX,NGLLZ,nspec), &
             permz(NGLLX,NGLLZ,nspec), &
             tau_e(NGLLX,NGLLZ,nspec,2), &
             tau_d(NGLLX,NGLLZ,nspec,2), &
             alphaval_em(NGLLX,NGLLZ,nspec,2), &
             betaval_em(NGLLX,NGLLZ,nspec,2), &
             gammaval_em(NGLLX,NGLLZ,nspec,2),stat=ier)
  else
    ! dummy allocation
    ! (rx_permattenuation,rz_permattenuation allocation needed for subroutine call)
    allocate(rx_permattenuation(1,1,1), &
             rz_permattenuation(1,1,1), &
             permx(1,1,1), &
             permz(1,1,1), &
             tau_e(1,1,1,1), &
             tau_d(1,1,1,1), &
             alphaval_em(1,1,1,1), &
             betaval_em(1,1,1,1), &
             gammaval_em(1,1,1,1),stat=ier)
  endif
  if (ier /= 0) stop 'Error allocating rx_permattenuation,.. arrays'
  ! initialize memory variables for attenuation
  rx_permattenuation(:,:,:) = 0.d0
  rz_permattenuation(:,:,:) = 0.d0
  permx(:,:,:) = 0.d0
  permz(:,:,:) = 0.d0

  end subroutine prepare_attenuation_global_arrays

!
!-------------------------------------------------------------------------------------
!


  subroutine prepare_attenuation_with_PML()

  use constants, only: NGLLX,NGLLZ,IMAIN,myrank, &
    CPML_X_ONLY,CPML_XZ,CPML_Z_ONLY

  use specfem_par, only: N_SLS,ATTENUATION_VISCOELASTIC, &
    ispec_is_elastic,nspec, &
    qmu_attenuation_store,qkappa_attenuation_store, &
    inv_tau_sigma_nu1,inv_tau_sigma_nu2

  ! PML
  use specfem_par, only: spec_to_PML,ispec_is_PML,region_CPML,min_distance_between_CPML_parameter, &
    K_x_store,d_x_store,alpha_x_store,K_z_store,d_z_store,alpha_z_store


  implicit none

  ! local parameters
  double precision, dimension(2*N_SLS) :: tauinvnu
  double precision, dimension(2*N_SLS + 1) :: tauinvplusone
  double precision :: qkappal,qmul
  double precision :: d_x, d_z, K_x, K_z, alpha_x, alpha_z, beta_x, beta_z
  double precision :: const_for_separation_two

  integer :: i,j,ispec,i_sls,ispec_PML

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  preparing attenuation within PML region'
    call flush_IMAIN()
  endif

  const_for_separation_two = min_distance_between_CPML_parameter * 2.d0

  do ispec = 1,nspec

    ispec_PML = spec_to_PML(ispec)

    if (ispec_is_PML(ispec)) then
      do j = 1,NGLLZ
        do i = 1,NGLLX
          ! kappa & mu
          qkappal = qkappa_attenuation_store(i,j,ispec)
          qmul = qmu_attenuation_store(i,j,ispec)

          ! checks if anything to do; values of Q == 9999. mean no attenuation
          if (qkappal > 9998.999d0 .and. qmul > 9998.999d0) cycle

          K_x = K_x_store(i,j,ispec_PML)
          d_x = d_x_store(i,j,ispec_PML)
          alpha_x = alpha_x_store(i,j,ispec_PML)

          K_z = K_z_store(i,j,ispec_PML)
          d_z = d_z_store(i,j,ispec_PML)
          alpha_z = alpha_z_store(i,j,ispec_PML)

          ! visco-elastic
          if (ATTENUATION_VISCOELASTIC) then
            if (region_CPML(ispec) == CPML_XZ) then
              if (ispec_is_elastic(ispec)) then
                do i_sls = 1,N_SLS
                  tauinvnu(i_sls) = inv_tau_sigma_nu1(i,j,ispec,i_sls)
                  tauinvnu(i_sls+N_SLS) = inv_tau_sigma_nu2(i,j,ispec,i_sls)
                enddo

                call tauinvnu_arange_from_lt_to_gt(2*N_SLS,tauinvnu)

                do i_sls = 1,2*N_SLS
                  if (abs(alpha_x - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                    alpha_x = tauinvnu(i_sls) + const_for_separation_two
                  endif
                  if (abs(alpha_x - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                    call stop_the_code('error in separation of alpha_x, tauinvnu')
                  endif
                enddo
                do i_sls = 1,2*N_SLS
                  tauinvplusone(i_sls) = tauinvnu(i_sls)
                enddo
                tauinvplusone(2*N_SLS + 1) = alpha_x

                call tauinvnu_arange_from_lt_to_gt(2 * N_SLS + 1,tauinvplusone)

                do i_sls = 1,2*N_SLS + 1
                  if (abs(alpha_z - tauinvplusone(i_sls)) < min_distance_between_CPML_parameter) then
                    alpha_z = tauinvplusone(i_sls) + const_for_separation_two
                  endif
                  if (abs(alpha_z - tauinvplusone(i_sls)) < min_distance_between_CPML_parameter) then
                    call stop_the_code('error in separation of alpha_z, alpha_x,tauinvnu')
                  endif
                enddo

                beta_z = alpha_z + d_z / K_z

                do i_sls = 1,2*N_SLS + 1
                  if (abs(beta_z - tauinvplusone(i_sls)) < min_distance_between_CPML_parameter) then
                    beta_z = tauinvplusone(i_sls) + const_for_separation_two
                  endif
                  if (abs(beta_z - tauinvplusone(i_sls)) < min_distance_between_CPML_parameter) then
                    call stop_the_code('error in separation of beta_z, alpha_x,tauinvnu')
                  endif
                enddo

                do i_sls = 1,2*N_SLS
                  tauinvplusone(i_sls) = tauinvnu(i_sls)
                enddo
                tauinvplusone(2*N_SLS + 1) = alpha_z

                call tauinvnu_arange_from_lt_to_gt(2 * N_SLS + 1,tauinvplusone)
                beta_x = alpha_x + d_x / K_x
                do i_sls = 1,2*N_SLS + 1
                  if (abs(beta_x - tauinvplusone(i_sls)) < min_distance_between_CPML_parameter) then
                    beta_x = tauinvplusone(i_sls) + const_for_separation_two
                  endif
                  if (abs(beta_x - tauinvplusone(i_sls)) < min_distance_between_CPML_parameter) then
                    call stop_the_code('error in separation of beta_x, alpha_z,tauinvnu')
                  endif
                enddo

                d_x = (beta_x - alpha_x) * K_x
                d_z = (beta_z - alpha_z) * K_z

                d_x_store(i,j,ispec_PML) = d_x
                alpha_x_store(i,j,ispec_PML) = alpha_x
                d_z_store(i,j,ispec_PML) = d_z
                alpha_z_store(i,j,ispec_PML) = alpha_z
              endif
            endif

            if (region_CPML(ispec) == CPML_X_ONLY) then
              do i_sls = 1,2*N_SLS
                if (abs(alpha_x - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                  alpha_x = tauinvnu(i_sls) + const_for_separation_two
                endif
                if (abs(alpha_x - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                  call stop_the_code('error in separation of alpha_x, tauinvnu')
                endif
              enddo
              beta_x = alpha_x + d_x / K_x
              do i_sls = 1,2*N_SLS
                if (abs(beta_x - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                  beta_x = tauinvnu(i_sls) + const_for_separation_two
                endif
                if (abs(beta_x - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                  call stop_the_code('error in separation of beta_x, tauinvnu')
                endif
              enddo
              d_x = (beta_x - alpha_x) * K_x
              d_x_store(i,j,ispec_PML) = d_x
              alpha_x_store(i,j,ispec_PML) = alpha_x
            endif

            if (region_CPML(ispec) == CPML_Z_ONLY) then
              do i_sls = 1,2*N_SLS
                if (abs(alpha_z - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                  alpha_z = tauinvnu(i_sls) + const_for_separation_two
                endif
                if (abs(alpha_z - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                  call stop_the_code('error in separation of alpha_z, tauinvnu')
                endif
              enddo
              beta_z = alpha_z + d_z / K_z
              do i_sls = 1,2*N_SLS
                if (abs(beta_z - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                  beta_z = tauinvnu(i_sls) + const_for_separation_two
                endif
                if (abs(beta_z - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                  call stop_the_code('error in separation of beta_z, tauinvnu')
                endif
              enddo
              d_z = (beta_z - alpha_z) * K_z
              d_z_store(i,j,ispec_PML) = d_z
              alpha_z_store(i,j,ispec_PML) = alpha_z
            endif

          endif ! ATTENUATION_VISCOELASTIC

        enddo ! NGLLX
      enddo ! NGLLZ

    endif ! ispec_is_PML

  enddo ! ispec

  ! user output
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  done PML attenuation setup'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  contains

    subroutine tauinvnu_arange_from_lt_to_gt(siz,tauinvnu)

    implicit none
    integer, intent(in) :: siz
    double precision, dimension(siz), intent(inout) :: tauinvnu

    !local parameters
    integer :: i,j
    double precision :: temp

    do i = 1,siz
      do j = i+1,siz
        if (tauinvnu(i) > tauinvnu(j)) then
          temp = tauinvnu(i)
          tauinvnu(i) = tauinvnu(j)
          tauinvnu(j) =  temp
        endif
      enddo
    enddo
    end subroutine tauinvnu_arange_from_lt_to_gt

  end subroutine prepare_attenuation_with_PML




