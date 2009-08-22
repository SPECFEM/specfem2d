
!========================================================================
!
!                   S P E C F E M 2 D  Version 6.0
!                   ------------------------------
!
! Copyright Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France,
! and Princeton University, USA.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
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

  subroutine compute_pressure_whole_medium(potential_dot_dot_acoustic,displ_elastic,&
         displs_poroelastic,displw_poroelastic,elastic,poroelastic,vector_field_display, &
         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,assign_external_model, &
         numat,kmato,density,porosity,tortuosity,elastcoef,vpext,vsext,rhoext,e1,e11, &
         TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON,Mu_nu1,Mu_nu2,N_SLS)

! compute pressure in acoustic elements and in elastic elements

  implicit none

  include "constants.h"

  integer :: nspec,npoin,numat

  integer, dimension(nspec) :: kmato
  integer, dimension(NGLLX,NGLLX,nspec) :: ibool

  double precision, dimension(2,numat) :: density
  double precision, dimension(numat) :: porosity,tortuosity
  double precision, dimension(4,3,numat) :: elastcoef
  double precision, dimension(NGLLX,NGLLX,nspec) :: vpext,vsext,rhoext

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz

  logical, dimension(nspec) :: elastic,poroelastic
  real(kind=CUSTOM_REAL), dimension(npoin) :: potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(NDIM,npoin) :: displ_elastic,displs_poroelastic,displw_poroelastic
  double precision, dimension(NDIM,npoin) :: vector_field_display

! array with derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz

  logical :: assign_external_model,TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON

  integer :: N_SLS
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec,N_SLS) :: e1,e11
  double precision, dimension(NGLLX,NGLLZ,nspec) :: Mu_nu1,Mu_nu2

! local variables
  integer :: i,j,ispec,iglob

! pressure in this element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: pressure_element

! loop over spectral elements
  do ispec = 1,nspec

! compute pressure in this element
    call compute_pressure_one_element(pressure_element,potential_dot_dot_acoustic,displ_elastic,&
         displs_poroelastic,displw_poroelastic,elastic,poroelastic,&
         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,assign_external_model, &
         numat,kmato,density,porosity,tortuosity,elastcoef,vpext,vsext,rhoext,ispec,e1,e11, &
         TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON,Mu_nu1,Mu_nu2,N_SLS)

! use vector_field_display as temporary storage, store pressure in its second component
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        vector_field_display(2,iglob) = pressure_element(i,j)
      enddo
    enddo

  enddo

  end subroutine compute_pressure_whole_medium

!
!=====================================================================
!

  subroutine compute_pressure_one_element(pressure_element,potential_dot_dot_acoustic,displ_elastic,&
         displs_poroelastic,displw_poroelastic,elastic,poroelastic,&
         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,assign_external_model, &
         numat,kmato,density,porosity,tortuosity,elastcoef,vpext,vsext,rhoext,ispec,e1,e11, &
         TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON,Mu_nu1,Mu_nu2,N_SLS)

! compute pressure in acoustic elements and in elastic elements

  implicit none

  include "constants.h"

  integer nspec,npoin,numat,ispec

  integer, dimension(nspec) :: kmato
  integer, dimension(NGLLX,NGLLX,nspec) :: ibool

  double precision, dimension(2,numat) :: density
  double precision, dimension(numat) :: porosity,tortuosity
  double precision, dimension(4,3,numat) :: elastcoef
  double precision, dimension(NGLLX,NGLLX,nspec) :: vpext,vsext,rhoext

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz

! pressure in this element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: pressure_element

  logical, dimension(nspec) :: elastic,poroelastic
  real(kind=CUSTOM_REAL), dimension(npoin) :: potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(NDIM,npoin) :: displ_elastic,displs_poroelastic,displw_poroelastic

! array with derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz

  logical :: assign_external_model,TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON

  integer :: N_SLS
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec,N_SLS) :: e1,e11
  real(kind=CUSTOM_REAL) :: e1_sum,e11_sum
  double precision, dimension(NGLLX,NGLLZ,nspec) :: Mu_nu1,Mu_nu2
  integer :: i_sls

! local variables
  integer :: i,j,k,iglob

! jacobian
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl

! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_zz,sigmap
  real(kind=CUSTOM_REAL) :: dwx_dxi,dwx_dgamma,dwz_dxi,dwz_dgamma
  real(kind=CUSTOM_REAL) :: dwx_dxl,dwz_dzl

! material properties of the elastic medium
  real(kind=CUSTOM_REAL) :: mul_relaxed,lambdal_relaxed,lambdalplus2mul_relaxed,denst
  real(kind=CUSTOM_REAL) :: mul_unrelaxed,lambdal_unrelaxed,lambdalplus2mul_unrelaxed,cpl,csl

  real(kind=CUSTOM_REAL) :: mul_s,kappal_s,rhol_s
  real(kind=CUSTOM_REAL) :: kappal_f,rhol_f
  real(kind=CUSTOM_REAL) :: mul_fr,kappal_fr,phil,tortl
  real(kind=CUSTOM_REAL) :: D_biot,H_biot,C_biot,M_biot,rhol_bar
  real(kind=CUSTOM_REAL) :: mul_G,lambdal_G,lambdalplus2mul_G

! if elastic element
!
! from L. S. Bennethum, Compressibility Moduli for Porous Materials Incorporating Volume Fraction,
! J. Engrg. Mech., vol. 132(11), p. 1205-1214 (2006), below equation (5):
! for a 3D isotropic solid, pressure is defined in terms of the trace of the stress tensor as
! p = -1/3 (t11 + t22 + t33) where t is the Cauchy stress tensor.

! to compute pressure in 3D in an elastic solid, one uses pressure = - trace(sigma) / 3
! sigma_ij = lambda delta_ij trace(epsilon) + 2 mu epsilon_ij
!          = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_ij
! sigma_xx = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_xx
! sigma_yy = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_yy
! sigma_zz = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_zz
! pressure = - trace(sigma) / 3 = - (lambda + 2/3 mu) trace(epsilon) = - kappa * trace(epsilon)
!
! to compute pressure in 2D in an elastic solid, one uses pressure = - trace(sigma) / 2
! sigma_ij = lambda delta_ij trace(epsilon) + 2 mu epsilon_ij
!          = lambda (epsilon_xx + epsilon_yy) + 2 mu epsilon_ij
! sigma_xx = lambda (epsilon_xx + epsilon_yy) + 2 mu epsilon_xx
! sigma_yy = lambda (epsilon_xx + epsilon_yy) + 2 mu epsilon_yy
! pressure = - trace(sigma) / 2 = - (lambda + mu) trace(epsilon)
!
  if(elastic(ispec)) then

! get relaxed elastic parameters of current spectral element
    lambdal_relaxed = elastcoef(1,1,kmato(ispec))
    mul_relaxed = elastcoef(2,1,kmato(ispec))
    lambdalplus2mul_relaxed = elastcoef(3,1,kmato(ispec))

    do j = 1,NGLLZ
      do i = 1,NGLLX

!--- if external medium, get elastic parameters of current grid point
        if(assign_external_model) then
          cpl = vpext(i,j,ispec)
          csl = vsext(i,j,ispec)
          denst = rhoext(i,j,ispec)
          mul_relaxed = denst*csl*csl
          lambdal_relaxed = denst*cpl*cpl - TWO*mul_relaxed
        endif

! derivative along x and along z
        dux_dxi = ZERO
        duz_dxi = ZERO

        dux_dgamma = ZERO
        duz_dgamma = ZERO

! first double loop over GLL points to compute and store gradients
! we can merge the two loops because NGLLX == NGLLZ
        do k = 1,NGLLX
          dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
          duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
          dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
          duz_dgamma = duz_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

! derivatives of displacement
        dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
        duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

! compute diagonal components of the stress tensor (include attenuation or anisotropy if needed)

  if(TURN_ATTENUATION_ON) then

! attenuation is implemented following the memory variable formulation of
! J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993). More details can be found in
! J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation in a linear
! viscoelastic medium, Geophysical Journal International, vol. 95, p. 597-611 (1988).

! compute unrelaxed elastic coefficients from formulas in Carcione 1993 page 111
    lambdal_unrelaxed = (lambdal_relaxed + mul_relaxed) * Mu_nu1(i,j,ispec) - mul_relaxed * Mu_nu2(i,j,ispec)
    mul_unrelaxed = mul_relaxed * Mu_nu2(i,j,ispec)
    lambdalplus2mul_unrelaxed = lambdal_unrelaxed + TWO*mul_unrelaxed

! compute the stress using the unrelaxed Lame parameters (Carcione 1993, page 111)
    sigma_xx = lambdalplus2mul_unrelaxed*dux_dxl + lambdal_unrelaxed*duz_dzl
    sigma_zz = lambdalplus2mul_unrelaxed*duz_dzl + lambdal_unrelaxed*dux_dxl

! add the memory variables using the relaxed parameters (Carcione 1993, page 111)
! beware: there is a bug in Carcione's equation (2c) for sigma_zz, we fixed it in the code below
    e1_sum = 0._CUSTOM_REAL
    e11_sum = 0._CUSTOM_REAL

    do i_sls = 1,N_SLS
      e1_sum = e1_sum + e1(i,j,ispec,i_sls)
      e11_sum = e11_sum + e11(i,j,ispec,i_sls)
    enddo

    sigma_xx = sigma_xx + (lambdal_relaxed + mul_relaxed) * e1_sum + TWO * mul_relaxed * e11_sum
    sigma_zz = sigma_zz + (lambdal_relaxed + mul_relaxed) * e1_sum - TWO * mul_relaxed * e11_sum

  else

! no attenuation
    sigma_xx = lambdalplus2mul_relaxed*dux_dxl + lambdal_relaxed*duz_dzl
    sigma_zz = lambdalplus2mul_relaxed*duz_dzl + lambdal_relaxed*dux_dxl

  endif

! full anisotropy
  if(TURN_ANISOTROPY_ON) then

! implement anisotropy in 2D
     sigma_xx = c11val*dux_dxl + c15val*(duz_dxl + dux_dzl) + c13val*duz_dzl
     sigma_zz = c13val*dux_dxl + c35val*(duz_dxl + dux_dzl) + c33val*duz_dzl

  endif

! store pressure
        pressure_element(i,j) = - (sigma_xx + sigma_zz) / 2.d0

      enddo
    enddo

  elseif(poroelastic(ispec)) then

! get poroelastic parameters of current spectral element
    phil = porosity(kmato(ispec))
    tortl = tortuosity(kmato(ispec))
!solid properties
    mul_s = elastcoef(2,1,kmato(ispec))
    kappal_s = elastcoef(3,1,kmato(ispec)) - FOUR_THIRDS*mul_s
    rhol_s = density(1,kmato(ispec))
!fluid properties
    kappal_f = elastcoef(1,2,kmato(ispec))
    rhol_f = density(2,kmato(ispec))
!frame properties
    mul_fr = elastcoef(2,3,kmato(ispec))
    kappal_fr = elastcoef(3,3,kmato(ispec)) - FOUR_THIRDS*mul_fr
    rhol_bar =  (1.d0 - phil)*rhol_s + phil*rhol_f
!Biot coefficients for the input phi
      D_biot = kappal_s*(1.d0 + phil*(kappal_s/kappal_f - 1.d0))
      H_biot = (kappal_s - kappal_fr)*(kappal_s - kappal_fr)/(D_biot - kappal_fr) + kappal_fr + FOUR_THIRDS*mul_fr
      C_biot = kappal_s*(kappal_s - kappal_fr)/(D_biot - kappal_fr)
      M_biot = kappal_s*kappal_s/(D_biot - kappal_fr)
!where T = G:grad u_s + C div w I
!and T_f = C div u_s I + M div w I
!we are expressing lambdaplus2mu, lambda, and mu for G, C, and M
      mul_G = mul_fr
      lambdal_G = H_biot - TWO*mul_fr
      lambdalplus2mul_G = lambdal_G + TWO*mul_G

    do j = 1,NGLLZ
      do i = 1,NGLLX

! derivative along x and along z
        dux_dxi = ZERO
        duz_dxi = ZERO

        dux_dgamma = ZERO
        duz_dgamma = ZERO

        dwx_dxi = ZERO
        dwz_dxi = ZERO

        dwx_dgamma = ZERO
        dwz_dgamma = ZERO

! first double loop over GLL points to compute and store gradients
! we can merge the two loops because NGLLX == NGLLZ
        do k = 1,NGLLX
          dux_dxi = dux_dxi + displs_poroelastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
          duz_dxi = duz_dxi + displs_poroelastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
          dux_dgamma = dux_dgamma + displs_poroelastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
          duz_dgamma = duz_dgamma + displs_poroelastic(2,ibool(i,k,ispec))*hprime_zz(j,k)

          dwx_dxi = dwx_dxi + displw_poroelastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
          dwz_dxi = dwz_dxi + displw_poroelastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
          dwx_dgamma = dwx_dgamma + displw_poroelastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
          dwz_dgamma = dwz_dgamma + displw_poroelastic(2,ibool(i,k,ispec))*hprime_zz(j,k)

        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

! derivatives of displacement
        dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
        duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

        dwx_dxl = dwx_dxi*xixl + dwx_dgamma*gammaxl
        dwz_dzl = dwz_dxi*xizl + dwz_dgamma*gammazl

! compute diagonal components of the stress tensor (include attenuation if needed)

  if(TURN_ATTENUATION_ON) then
!-------------------- ATTENTION TO BE DEFINED ------------------------------!

! attenuation is implemented following the memory variable formulation of
! J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993). More details can be found in
! J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation in a linear
! viscoelastic medium, Geophysical Journal International, vol. 95, p. 597-611 (1988).

! compute unrelaxed elastic coefficients from formulas in Carcione 1993 page 111
    lambdal_unrelaxed = (lambdal_relaxed + mul_relaxed) * Mu_nu1(i,j,ispec) - mul_relaxed * Mu_nu2(i,j,ispec)
    mul_unrelaxed = mul_relaxed * Mu_nu2(i,j,ispec)
    lambdalplus2mul_unrelaxed = lambdal_unrelaxed + TWO*mul_unrelaxed

! compute the stress using the unrelaxed Lame parameters (Carcione 1993, page 111)
    sigma_xx = lambdalplus2mul_unrelaxed*dux_dxl + lambdal_unrelaxed*duz_dzl
    sigma_zz = lambdalplus2mul_unrelaxed*duz_dzl + lambdal_unrelaxed*dux_dxl

! add the memory variables using the relaxed parameters (Carcione 1993, page 111)
! beware: there is a bug in Carcione's equation (2c) for sigma_zz, we fixed it in the code below
    e1_sum = 0._CUSTOM_REAL
    e11_sum = 0._CUSTOM_REAL

    do i_sls = 1,N_SLS
      e1_sum = e1_sum + e1(i,j,ispec,i_sls)
      e11_sum = e11_sum + e11(i,j,ispec,i_sls)
    enddo

    sigma_xx = sigma_xx + (lambdal_relaxed + mul_relaxed) * e1_sum + TWO * mul_relaxed * e11_sum
    sigma_zz = sigma_zz + (lambdal_relaxed + mul_relaxed) * e1_sum - TWO * mul_relaxed * e11_sum

  else

! no attenuation
    sigma_xx = lambdalplus2mul_G*dux_dxl + lambdal_G*duz_dzl + C_biot*(dwx_dxl + dwz_dzl)
    sigma_zz = lambdalplus2mul_G*duz_dzl + lambdal_G*dux_dxl + C_biot*(dwx_dxl + dwz_dzl)

    sigmap = C_biot*(dux_dxl + duz_dzl) + M_biot*(dwx_dxl + dwz_dzl)

  endif

! store pressure
        pressure_element(i,j) = - (sigma_xx + sigma_zz) / 2.d0
!        pressure_element2(i,j) = - sigmap
      enddo
    enddo

! pressure = - Chi_dot_dot if acoustic element
  else

    do j = 1,NGLLZ
      do i = 1,NGLLX

        iglob = ibool(i,j,ispec)

! store pressure
        pressure_element(i,j) = - potential_dot_dot_acoustic(iglob)

      enddo
    enddo

  endif ! end of test if acoustic or elastic element

  end subroutine compute_pressure_one_element

