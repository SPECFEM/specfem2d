
subroutine prepare_timerun_attenuation(e1,e11,e13,inv_tau_sigma_nu1, &
     inv_tau_sigma_nu2,Mu_nu1,Mu_nu2, phi_nu1, phi_nu2,nspec)

  ! allocate memory variables for attenuation
  allocate(e1(NGLLX,NGLLZ,nspec_allocate,N_SLS), &
       e11(NGLLX,NGLLZ,nspec_allocate,N_SLS), &
       e13(NGLLX,NGLLZ,nspec_allocate,N_SLS),stat=ier)
  if (ier /= 0) stop 'Error allocating attenuation arrays'

  e1(:,:,:,:) = 0._CUSTOM_REAL
  e11(:,:,:,:) = 0._CUSTOM_REAL
  e13(:,:,:,:) = 0._CUSTOM_REAL


  ! initialize to dummy values
  ! convention to indicate that Q = 9999 in that element i.e. that there is no viscoelasticity in that element
  inv_tau_sigma_nu1(:,:,:,:) = -1._CUSTOM_REAL
  phi_nu1(:,:,:,:) = -1._CUSTOM_REAL
  inv_tau_sigma_nu2(:,:,:,:) = -1._CUSTOM_REAL
  phi_nu2(:,:,:,:) = -1._CUSTOM_REAL
  Mu_nu1(:,:,:) = -1._CUSTOM_REAL
  Mu_nu2(:,:,:) = -1._CUSTOM_REAL

  do ispec = 1,nspec



     ! check that attenuation values entered by the user make sense
     if ((QKappa_attenuation(kmato(ispec)) <= 9998.999d0 .and. Qmu_attenuation(kmato(ispec)) > 9998.999d0) .or. &
          (QKappa_attenuation(kmato(ispec)) > 9998.999d0 .and. Qmu_attenuation(kmato(ispec)) <= 9998.999d0)) stop &
          'need to have Qkappa and Qmu both above or both below 9999 for a given material; trick: use 9998 if you want to turn off one'

     ! if no attenuation in that elastic element
     if (QKappa_attenuation(kmato(ispec)) > 9998.999d0) cycle

     call attenuation_model(QKappa_attenuation(kmato(ispec)),Qmu_attenuation(kmato(ispec)))
     ! si acosutique il faut appeler une autre fonction
     !   call attenuation_viscoacoustic_model ! todo

     do j = 1,NGLLZ
        do i = 1,NGLLX
           inv_tau_sigma_nu1(i,j,ispec,:) = inv_tau_sigma_nu1_sent(:)
           phi_nu1(i,j,ispec,:) = phi_nu1_sent(:)
           inv_tau_sigma_nu2(i,j,ispec,:) = inv_tau_sigma_nu2_sent(:) !0 si acoustique
           phi_nu2(i,j,ispec,:) = phi_nu2_sent(:)!0 si acoustique
           Mu_nu1(i,j,ispec) = Mu_nu1_sent
           Mu_nu2(i,j,ispec) = Mu_nu2_sent !0 si acoustique
        enddo
     enddo

     ! todo
     !   call shift_velocities_from_f0(vp,vs,rhol,mul,lambdal)

  enddo
end subroutine prepare_timerun_attenuation



subroutine compute_forces_viscoelastic(accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old, &
    e1,e11,e13)

  ! compute Grad(displ_elastic) at time step n for attenuation
  call compute_gradient_attenuation(displ_elastic,dux_dxl_n,duz_dxl_n, &
       dux_dzl_n,duz_dzl_n,xix,xiz,gammax,gammaz,ibool,ispec_is_elastic,hprime_xx,hprime_zz,nspec,nglob)

  ! je vais avoir besoin de :
  !  displ_elastic_old(:,:) = displ_elastic(:,:) + deltatsquareover2/TWO * accel_elastic(:,:) !cf update_displacement scheme l=518

  ! compute Grad(disp_elastic_old) at time step n-1 for attenuation
  call compute_gradient_attenuation(displ_elastic_old,dux_dxl_nsub1,duz_dxl_nsub1, &
       dux_dzl_nsub1,duz_dzl_nsub1,xix,xiz,gammax,gammaz,ibool,ispec_is_elastic,hprime_xx,hprime_zz,nspec,nglob)

  !update memory variable in viscoelastic simulation
  ! loop over spectral elements and GLL points
  do ispec = 1,nspec
     do j = 1,NGLLZ
        do i = 1,NGLLX

           theta_n_u = dux_dxl_n(i,j,ispec) + duz_dzl_n(i,j,ispec)
           theta_nsub1_u = dux_dxl_nsub1(i,j,ispec) + duz_dzl_nsub1(i,j,ispec)

           do i_sls = 1,N_SLS
              phinu1 = phi_nu1(i,j,ispec,i_sls)
              tauinvnu1 = inv_tau_sigma_nu1(i,j,ispec,i_sls)
              phinu2 = phi_nu2(i,j,ispec,i_sls)
              tauinvnu2 = inv_tau_sigma_nu2(i,j,ispec,i_sls)

              ! update e1, e11, e13 in convolution formation with modified recursive convolution scheme on basis of
              ! second-order accurate convolution term calculation from equation (21) of
              ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
              ! Anisotropic-Medium PML for Vector FETD With Modified Basis functions,
              ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

              if (CONVOLUTION_MEMORY_VARIABLES) then
                 call compute_coef_convolution(tauinvnu1,deltat,coef0,coef1,coef2)

                 e1(i,j,ispec,i_sls) = coef0 * e1(i,j,ispec,i_sls) + &
                      phinu1 * (coef1 * theta_n_u + coef2 * theta_nsub1_u)

                 call compute_coef_convolution(tauinvnu2,deltat,coef0,coef1,coef2)

                 e11(i,j,ispec,i_sls) = coef0 * e11(i,j,ispec,i_sls) + &
                      phinu2 * (coef1 * (dux_dxl_n(i,j,ispec)-theta_n_u/TWO) + &
                      coef2 * (dux_dxl_nsub1(i,j,ispec)-theta_nsub1_u/TWO))

                 e13(i,j,ispec,i_sls) = coef0 * e13(i,j,ispec,i_sls) + &
                      phinu2 * (coef1 * (dux_dzl_n(i,j,ispec) + duz_dxl_n(i,j,ispec)) + &
                      coef2 * (dux_dzl_nsub1(i,j,ispec) + duz_dxl_nsub1(i,j,ispec)))
              else
                 e1(i,j,ispec,i_sls) = e1(i,j,ispec,i_sls) + deltat * &
                      (- e1(i,j,ispec,i_sls)*tauinvnu1 + phinu1 * theta_n_u)

                 e11(i,j,ispec,i_sls) = e11(i,j,ispec,i_sls) + deltat * &
                      (- e11(i,j,ispec,i_sls)*tauinvnu2 + phinu2 * (dux_dxl_n(i,j,ispec)-theta_n_u/TWO))

                 e13(i,j,ispec,i_sls) = e13(i,j,ispec,i_sls) + deltat * &
                      (- e13(i,j,ispec,i_sls)*tauinvnu2 + phinu2 * (dux_dzl_n(i,j,ispec) + duz_dxl_n(i,j,ispec)))
              endif

            endif
        enddo
     enddo
  enddo  ! end of update memory variable in viscoelastic simulation

  do ispec = 1,nspec
     do j = 1,NGLLZ
        do i = 1,NGLLX

           e1_sum = 0._CUSTOM_REAL; e11_sum = 0._CUSTOM_REAL;  e13_sum = 0._CUSTOM_REAL
           do i_sls = 1,N_SLS
              e1_sum = e1_sum + e1(i,j,ispec,i_sls)
              e11_sum = e11_sum + e11(i,j,ispec,i_sls)
              e13_sum = e13_sum + e13(i,j,ispec,i_sls)
           enddo

           sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
                    sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
                    sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dx

           ! use the right formula with 1/N included
           ! i.e. use the unrelaxed moduli here (see Carcione's book, third edition, equation (3.189))
           sigma_xx = sigma_xx + lambdalplusmul_unrelaxed_elastic * e1_sum + TWO * mul_unrelaxed_elastic * e11_sum
           sigma_xz = sigma_xz + mul_unrelaxed_elastic * e13_sum
           sigma_zz = sigma_zz + lambdalplusmul_unrelaxed_elastic * e1_sum - TWO * mul_unrelaxed_elastic * e11_sum
           sigma_zx = sigma_xz

        enddo
     enddo
  enddo


end subroutine compute_forces_viscoelastic








subroutine attenuation_model(QKappa_att,QMu_att)

  ! define the attenuation constants

  implicit none

  double precision :: QKappa_att,QMu_att

  double precision, dimension(N_SLS) :: tau_epsilon_nu1_d,tau_sigma_nu1,tau_epsilon_nu2_d,tau_sigma_nu2

  double precision :: f_min_attenuation, f_max_attenuation

  ! attenuation constants for standard linear solids
  ! nu1 is the dilatation/incompressibility mode (QKappa)
  ! nu2 is the shear mode (Qmu)
  ! array index (1) is the first standard linear solid, (2) is the second etc.

  ! use a wide bandwidth (always OK when using three or more Standard Linear Solids, can be a bit inaccurate if using only two)
  f_min_attenuation = f0_attenuation / 10.d0
  f_max_attenuation = f0_attenuation * 10.d0

  ! si mileu acoustique cette fonction est apellee q'une seule fois avec QKappa uniquement et nu=nu1

  call compute_attenuation_coeffs(N_SLS,QKappa_att,f0_attenuation,f_min_attenuation,f_max_attenuation, &
       tau_epsilon_nu1_d,tau_sigma_nu1)

  call compute_attenuation_coeffs(N_SLS,QMu_att,f0_attenuation,f_min_attenuation,f_max_attenuation, &
       tau_epsilon_nu2_d,tau_sigma_nu2)



  if (any(tau_sigma_nu1 < 0.d0) .or. any(tau_sigma_nu2 < 0.d0) .or. &
       any(tau_epsilon_nu1_d < 0.d0) .or. any(tau_epsilon_nu2_d < 0.d0)) &
       stop 'error: negative relaxation time found for a viscoelastic material'

  ! in the old formulation of Carcione 1993, which is based on Liu et al. 1976, the 1/N factor is missing
  ! and thus this does not apply; it only applies to the right formula with 1/N included


  ! When implementing viscoelasticity according to the Carcione 1993 paper, attenuation is
  ! non-causal rather than causal i.e. wave speed up instead of slowing down
  ! when attenuation is turned on. We fixed that issue (which is not incorrect but non traditional)
  ! by taking the unrelaxed state (infinite frequency) as a reference instead of the relaxed state (zero frequency)
  ! and also using equations in Carcione's 2007 book.
  ! See file doc/old_problem_attenuation_reference_Specfem2D_fixed_by_Xie_Zhinan.pdf
  ! and doc/how_we_modified_Carcione_1993_to_make_it_causal_and_include_the_missing_1_over_L_factor.pdf

  ! See also J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
  ! vol. 58(1), p. 110-120 (1993) for two memory-variable mechanisms (page 112).

  ! and J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation
  ! in a linear viscoelastic medium, Geophysical Journal International,
  ! vol. 95, p. 597-611 (1988) for two memory-variable mechanisms (page 604).

  !
  !--- other constants computed from the parameters above, do not modify
  !

  tau_epsilon_nu1(:) = real(tau_epsilon_nu1_d(:),kind=CUSTOM_REAL)
  tau_epsilon_nu2(:) = real(tau_epsilon_nu2_d(:),kind=CUSTOM_REAL) !0 en acoustique

  inv_tau_sigma_nu1_sent(:) = real(dble(ONE) / tau_sigma_nu1(:),kind=CUSTOM_REAL)
  inv_tau_sigma_nu2_sent(:) = real(dble(ONE) / tau_sigma_nu2(:),kind=CUSTOM_REAL) ! 0 en acoustique

  ! use the right formula with 1/N included
  phi_nu1_sent(:) = real((dble(ONE) - tau_epsilon_nu1_d(:)/tau_sigma_nu1(:)) / tau_sigma_nu1(:) &
       / sum(tau_epsilon_nu1_d/tau_sigma_nu1),kind=CUSTOM_REAL)
  phi_nu2_sent(:) = real((dble(ONE) - tau_epsilon_nu2_d(:)/tau_sigma_nu2(:)) / tau_sigma_nu2(:) &
       / sum(tau_epsilon_nu2_d/tau_sigma_nu2),kind=CUSTOM_REAL) !zero en acoustique

  Mu_nu1_sent = real(sum(tau_epsilon_nu1_d/tau_sigma_nu1) / dble(N_SLS),kind=CUSTOM_REAL)
  Mu_nu2_sent = real(sum(tau_epsilon_nu2_d/tau_sigma_nu2) / dble(N_SLS),kind=CUSTOM_REAL) !zero en acoustique

  if (Mu_nu1_sent < 1. .or. Mu_nu2_sent < 1.) &
       stop 'error in Zener viscoelasticity: must have Mu_nu1 and Mu_nu2 both greater than one'

end subroutine attenuation_model
