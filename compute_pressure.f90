
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
!                         Dimitri Komatitsch
!                     University of Pau, France
!
!                          (c) April 2007
!
!========================================================================

  subroutine compute_pressure_whole_medium(potential_dot_dot_acoustic,displ_elastic,elastic,vector_field_display, &
         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,assign_external_model, &
         numat,kmato,density,elastcoef,vpext,vsext,rhoext,e1_mech1,e11_mech1, &
         e1_mech2,e11_mech2,TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON)

! compute pressure in acoustic elements and in elastic elements

  implicit none

  include "constants.h"

  integer :: nspec,npoin,numat

  integer, dimension(nspec) :: kmato
  integer, dimension(NGLLX,NGLLX,nspec) :: ibool

  double precision, dimension(numat) :: density
  double precision, dimension(4,numat) :: elastcoef
  double precision, dimension(NGLLX,NGLLX,nspec) :: vpext,vsext,rhoext

  double precision, dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz

  logical, dimension(nspec) :: elastic
  double precision, dimension(npoin) :: potential_dot_dot_acoustic
  double precision, dimension(NDIM,npoin) :: displ_elastic,vector_field_display

! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLLX,NGLLX) :: hprime_xx
  double precision, dimension(NGLLZ,NGLLZ) :: hprime_zz

  logical :: assign_external_model,TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON

  double precision, dimension(NGLLX,NGLLZ,nspec) :: e1_mech1,e11_mech1,e1_mech2,e11_mech2

! local variables
  integer :: i,j,ispec,iglob

! pressure in this element
  double precision, dimension(NGLLX,NGLLX) :: pressure_element

! loop over spectral elements
  do ispec = 1,nspec

! compute pressure in this element
    call compute_pressure_one_element(pressure_element,potential_dot_dot_acoustic,displ_elastic,elastic, &
         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,assign_external_model, &
         numat,kmato,density,elastcoef,vpext,vsext,rhoext,ispec,e1_mech1,e11_mech1, &
         e1_mech2,e11_mech2,TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON)

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

  subroutine compute_pressure_one_element(pressure_element,potential_dot_dot_acoustic,displ_elastic,elastic, &
         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,assign_external_model, &
         numat,kmato,density,elastcoef,vpext,vsext,rhoext,ispec,e1_mech1,e11_mech1, &
         e1_mech2,e11_mech2,TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON)

! compute pressure in acoustic elements and in elastic elements

  implicit none

  include "constants.h"

  integer nspec,npoin,numat,ispec

  integer, dimension(nspec) :: kmato
  integer, dimension(NGLLX,NGLLX,nspec) :: ibool

  double precision, dimension(numat) :: density
  double precision, dimension(4,numat) :: elastcoef
  double precision, dimension(NGLLX,NGLLX,nspec) :: vpext,vsext,rhoext

  double precision, dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz

! pressure in this element
  double precision, dimension(NGLLX,NGLLX) :: pressure_element

  logical, dimension(nspec) :: elastic
  double precision, dimension(npoin) :: potential_dot_dot_acoustic
  double precision, dimension(NDIM,npoin) :: displ_elastic

! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLLX,NGLLX) :: hprime_xx
  double precision, dimension(NGLLZ,NGLLZ) :: hprime_zz

  logical :: assign_external_model,TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON

  double precision, dimension(NGLLX,NGLLZ,nspec) :: e1_mech1,e11_mech1,e1_mech2,e11_mech2

! local variables
  integer :: i,j,k,iglob

! jacobian
  double precision :: xixl,xizl,gammaxl,gammazl

! spatial derivatives
  double precision :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  double precision :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  double precision :: sigma_xx,sigma_zz

! material properties of the elastic medium
  integer :: material
  double precision :: mul_relaxed,lambdal_relaxed,lambdalplus2mul_relaxed,denst
  double precision :: mul_unrelaxed,lambdal_unrelaxed,lambdalplus2mul_unrelaxed,cpl,csl

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
    lambdal_relaxed = elastcoef(1,kmato(ispec))
    mul_relaxed = elastcoef(2,kmato(ispec))
    lambdalplus2mul_relaxed = elastcoef(3,kmato(ispec))

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
    lambdal_unrelaxed = (lambdal_relaxed + mul_relaxed) * Mu_nu1 - mul_relaxed * Mu_nu2
    mul_unrelaxed = mul_relaxed * Mu_nu2
    lambdalplus2mul_unrelaxed = lambdal_unrelaxed + TWO*mul_unrelaxed

! compute the stress using the unrelaxed Lame parameters (Carcione 1993, page 111)
    sigma_xx = lambdalplus2mul_unrelaxed*dux_dxl + lambdal_unrelaxed*duz_dzl
    sigma_zz = lambdalplus2mul_unrelaxed*duz_dzl + lambdal_unrelaxed*dux_dxl

! add the memory variables using the relaxed parameters (Carcione 1993, page 111)
! beware: there is a bug in Carcione's equation (2c) for sigma_zz, we fixed it in the code below
    sigma_xx = sigma_xx + (lambdal_relaxed + mul_relaxed)* &
      (e1_mech1(i,j,ispec) + e1_mech2(i,j,ispec)) + TWO * mul_relaxed * (e11_mech1(i,j,ispec) + e11_mech2(i,j,ispec))
    sigma_zz = sigma_zz + (lambdal_relaxed + mul_relaxed)* &
      (e1_mech1(i,j,ispec) + e1_mech2(i,j,ispec)) - TWO * mul_relaxed * (e11_mech1(i,j,ispec) + e11_mech2(i,j,ispec))

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

! pressure = - rho * Chi_dot_dot if acoustic element
  else

    do j = 1,NGLLZ
      do i = 1,NGLLX

        iglob = ibool(i,j,ispec)

        material = kmato(ispec)
        denst = density(material)
        if(assign_external_model) denst = rhoext(i,j,ispec)

! store pressure
        pressure_element(i,j) = - denst * potential_dot_dot_acoustic(iglob)

      enddo
    enddo

  endif ! end of test if acoustic or elastic element

  end subroutine compute_pressure_one_element

