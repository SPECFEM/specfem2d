
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

  subroutine compute_forces_elastic(npoin,nspec,nelemabs,numat,iglob_source, &
               ispec_selected_source,source_type,it,NSTEP,anyabs,assign_external_model, &
               initialfield,TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON,angleforce,deltatcube, &
               deltatfourth,twelvedeltat,fourdeltatsquare,ibool,kmato,numabs,elastic,codeabs, &
               accel_elastic,veloc_elastic,displ_elastic,density,elastcoef,xix,xiz,gammax,gammaz, &
               jacobian,vpext,vsext,rhoext,source_time_function,sourcearray,e1_mech1,e11_mech1, &
               e13_mech1,e1_mech2,e11_mech2,e13_mech2,dux_dxl_n,duz_dzl_n,duz_dxl_n,dux_dzl_n, &
               dux_dxl_np1,duz_dzl_np1,duz_dxl_np1,dux_dzl_np1,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll)

! compute forces for the elastic elements

  implicit none

  include "constants.h"

  integer :: npoin,nspec,nelemabs,numat,iglob_source,ispec_selected_source,source_type,it,NSTEP

  logical :: anyabs,assign_external_model,initialfield,TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON

  double precision :: angleforce,deltatcube,deltatfourth,twelvedeltat,fourdeltatsquare

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: kmato
  integer, dimension(nelemabs) :: numabs

  logical, dimension(nspec) :: elastic
  logical, dimension(4,nelemabs)  :: codeabs

  double precision, dimension(NDIM,npoin) :: accel_elastic,veloc_elastic,displ_elastic
  double precision, dimension(numat) :: density
  double precision, dimension(4,numat) :: elastcoef
  double precision, dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz,jacobian
  double precision, dimension(NGLLX,NGLLZ,nspec) :: vpext,vsext,rhoext
  double precision, dimension(NSTEP) :: source_time_function
  double precision, dimension(NDIM,NGLLX,NGLLZ) :: sourcearray

  double precision, dimension(NGLLX,NGLLZ,nspec) :: &
    e1_mech1,e11_mech1,e13_mech1,e1_mech2,e11_mech2,e13_mech2, &
    dux_dxl_n,duz_dzl_n,duz_dxl_n,dux_dzl_n,dux_dxl_np1,duz_dzl_np1,duz_dxl_np1,dux_dzl_np1

! derivatives of Lagrange polynomials
  double precision, dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  double precision, dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz

! Gauss-Lobatto-Legendre weights
  double precision, dimension(NGLLX) :: wxgll
  double precision, dimension(NGLLZ) :: wzgll

!---
!--- local variables
!---

  integer :: ispec,i,j,k,iglob,ispecabs,ibegin,iend

! spatial derivatives
  double precision :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  double precision :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  double precision :: sigma_xx,sigma_xz,sigma_zz
  double precision :: nx,nz,vx,vz,vn,rho_vp,rho_vs,tx,tz,weight,xxi,zxi,xgamma,zgamma,jacobian1D

  double precision, dimension(NGLLX,NGLLZ) :: tempx1,tempx2,tempz1,tempz2

! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl,jacobianl

! material properties of the elastic medium
  double precision :: mul_relaxed,lambdal_relaxed,lambdalplus2mul_relaxed,kappal,cpl,csl,rhol, &
                      lambdal_unrelaxed,mul_unrelaxed,lambdalplus2mul_unrelaxed

! for attenuation
  double precision :: Un,Unp1,tauinv,Sn,Snp1,theta_n,theta_np1,tauinvsquare,tauinvcube,tauinvUn

! compute Grad(displ_elastic) at time step n for attenuation
  if(TURN_ATTENUATION_ON) call compute_gradient_attenuation(displ_elastic,dux_dxl_n,duz_dxl_n, &
      dux_dzl_n,duz_dzl_n,xix,xiz,gammax,gammaz,ibool,elastic,hprime_xx,hprime_zz,nspec,npoin)

!   loop over spectral elements
    do ispec = 1,nspec

!---
!--- elastic spectral element
!---
    if(elastic(ispec)) then

! get relaxed elastic parameters of current spectral element
      lambdal_relaxed = elastcoef(1,kmato(ispec))
      mul_relaxed = elastcoef(2,kmato(ispec))
      lambdalplus2mul_relaxed = elastcoef(3,kmato(ispec))

! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX

!--- if external medium, get elastic parameters of current grid point
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
            csl = vsext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
            mul_relaxed = rhol*csl*csl
            lambdal_relaxed = rhol*cpl*cpl - TWO*mul_relaxed
            lambdalplus2mul_relaxed = lambdal_relaxed + TWO*mul_relaxed
          endif

! derivative along x and along z
          dux_dxi = ZERO
          duz_dxi = ZERO

          dux_dgamma = ZERO
          duz_dgamma = ZERO

! first double loop over GLL points to compute and store gradients
! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(k,i)
            duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(k,i)
            dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(k,j)
            duz_dgamma = duz_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(k,j)
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

! derivatives of displacement
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

! compute stress tensor (include attenuation or anisotropy if needed)

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
    sigma_xz = mul_unrelaxed*(duz_dxl + dux_dzl)
    sigma_zz = lambdalplus2mul_unrelaxed*duz_dzl + lambdal_unrelaxed*dux_dxl

! add the memory variables using the relaxed parameters (Carcione 1993, page 111)
! beware: there is a bug in Carcione's equation (2c) for sigma_zz, we fixed it in the code below
    sigma_xx = sigma_xx + (lambdal_relaxed + mul_relaxed)* &
      (e1_mech1(i,j,ispec) + e1_mech2(i,j,ispec)) + TWO * mul_relaxed * (e11_mech1(i,j,ispec) + e11_mech2(i,j,ispec))
    sigma_xz = sigma_xz + mul_relaxed * (e13_mech1(i,j,ispec) + e13_mech2(i,j,ispec))
    sigma_zz = sigma_zz + (lambdal_relaxed + mul_relaxed)* &
      (e1_mech1(i,j,ispec) + e1_mech2(i,j,ispec)) - TWO * mul_relaxed * (e11_mech1(i,j,ispec) + e11_mech2(i,j,ispec))

  else

! no attenuation
    sigma_xx = lambdalplus2mul_relaxed*dux_dxl + lambdal_relaxed*duz_dzl
    sigma_xz = mul_relaxed*(duz_dxl + dux_dzl)
    sigma_zz = lambdalplus2mul_relaxed*duz_dzl + lambdal_relaxed*dux_dxl

  endif

! full anisotropy
  if(TURN_ANISOTROPY_ON) then

! implement anisotropy in 2D
     sigma_xx = c11val*dux_dxl + c15val*(duz_dxl + dux_dzl) + c13val*duz_dzl
     sigma_zz = c13val*dux_dxl + c35val*(duz_dxl + dux_dzl) + c33val*duz_dzl
     sigma_xz = c15val*dux_dxl + c55val*(duz_dxl + dux_dzl) + c35val*duz_dzl

  endif

          jacobianl = jacobian(i,j,ispec)

! weak formulation term based on stress tensor (non-symmetric form)
! also add GLL integration weights
          tempx1(i,j) = wzgll(j)*jacobianl*(sigma_xx*xixl+sigma_xz*xizl)
          tempz1(i,j) = wzgll(j)*jacobianl*(sigma_xz*xixl+sigma_zz*xizl)

          tempx2(i,j) = wxgll(i)*jacobianl*(sigma_xx*gammaxl+sigma_xz*gammazl)
          tempz2(i,j) = wxgll(i)*jacobianl*(sigma_xz*gammaxl+sigma_zz*gammazl)

        enddo
      enddo

!
! second double-loop over GLL to compute all the terms
!
      do j = 1,NGLLZ
        do i = 1,NGLLX

          iglob = ibool(i,j,ispec)

! along x direction and z direction
! and assemble the contributions
! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tempx1(k,j)*hprimewgll_xx(i,k) + tempx2(i,k)*hprimewgll_zz(j,k))
            accel_elastic(2,iglob) = accel_elastic(2,iglob) - (tempz1(k,j)*hprimewgll_xx(i,k) + tempz2(i,k)*hprimewgll_zz(j,k))
          enddo

        enddo ! second loop over the GLL points
      enddo

    endif ! end of test if elastic element

    enddo ! end of loop over all spectral elements

!
!--- absorbing boundaries
!
  if(anyabs) then

    do ispecabs = 1,nelemabs

      ispec = numabs(ispecabs)

! get elastic parameters of current spectral element
      lambdal_relaxed = elastcoef(1,kmato(ispec))
      mul_relaxed = elastcoef(2,kmato(ispec))
      rhol  = density(kmato(ispec))
      kappal  = lambdal_relaxed + TWO*mul_relaxed/3.d0
      cpl = sqrt((kappal + 4.d0*mul_relaxed/3.d0)/rhol)
      csl = sqrt(mul_relaxed/rhol)

!--- left absorbing boundary
      if(codeabs(ILEFT,ispecabs)) then

        i = 1

        do j = 1,NGLLZ

          iglob = ibool(i,j,ispec)

! external velocity model
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
            csl = vsext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)
          nx = + zgamma / jacobian1D
          nz = - xgamma / jacobian1D

          weight = jacobian1D * wzgll(j)

! Clayton-Engquist condition if elastic
          if(elastic(ispec)) then
            vx = veloc_elastic(1,iglob)
            vz = veloc_elastic(2,iglob)

            vn = nx*vx+nz*vz

            tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
            tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

            accel_elastic(1,iglob) = accel_elastic(1,iglob) - tx*weight
            accel_elastic(2,iglob) = accel_elastic(2,iglob) - tz*weight
          endif

        enddo

      endif  !  end of left absorbing boundary

!--- right absorbing boundary
      if(codeabs(IRIGHT,ispecabs)) then

        i = NGLLX

        do j = 1,NGLLZ

          iglob = ibool(i,j,ispec)

! external velocity model
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
            csl = vsext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)
          nx = + zgamma / jacobian1D
          nz = - xgamma / jacobian1D

          weight = jacobian1D * wzgll(j)

! Clayton-Engquist condition if elastic
          if(elastic(ispec)) then
            vx = veloc_elastic(1,iglob)
            vz = veloc_elastic(2,iglob)

            vn = nx*vx+nz*vz

            tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
            tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

            accel_elastic(1,iglob) = accel_elastic(1,iglob) - tx*weight
            accel_elastic(2,iglob) = accel_elastic(2,iglob) - tz*weight
          endif

        enddo

      endif  !  end of right absorbing boundary

!--- bottom absorbing boundary
      if(codeabs(IBOTTOM,ispecabs)) then

        j = 1

! exclude corners to make sure there is no contradiction on the normal
        ibegin = 1
        iend = NGLLX
        if(codeabs(ILEFT,ispecabs)) ibegin = 2
        if(codeabs(IRIGHT,ispecabs)) iend = NGLLX-1

        do i = ibegin,iend

          iglob = ibool(i,j,ispec)

! external velocity model
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
            csl = vsext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)
          nx = + zxi / jacobian1D
          nz = - xxi / jacobian1D

          weight = jacobian1D * wxgll(i)

! Clayton-Engquist condition if elastic
          if(elastic(ispec)) then
            vx = veloc_elastic(1,iglob)
            vz = veloc_elastic(2,iglob)

            vn = nx*vx+nz*vz

            tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
            tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

            accel_elastic(1,iglob) = accel_elastic(1,iglob) - tx*weight
            accel_elastic(2,iglob) = accel_elastic(2,iglob) - tz*weight
          endif

        enddo

      endif  !  end of bottom absorbing boundary

!--- top absorbing boundary
      if(codeabs(ITOP,ispecabs)) then

        j = NGLLZ

! exclude corners to make sure there is no contradiction on the normal
        ibegin = 1
        iend = NGLLX
        if(codeabs(ILEFT,ispecabs)) ibegin = 2
        if(codeabs(IRIGHT,ispecabs)) iend = NGLLX-1

        do i = ibegin,iend

          iglob = ibool(i,j,ispec)

! external velocity model
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
            csl = vsext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)
          nx = + zxi / jacobian1D
          nz = - xxi / jacobian1D

          weight = jacobian1D * wxgll(i)

! Clayton-Engquist condition if elastic
          if(elastic(ispec)) then
            vx = veloc_elastic(1,iglob)
            vz = veloc_elastic(2,iglob)

            vn = nx*vx+nz*vz

            tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
            tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

            accel_elastic(1,iglob) = accel_elastic(1,iglob) - tx*weight
            accel_elastic(2,iglob) = accel_elastic(2,iglob) - tz*weight
          endif

        enddo

      endif  !  end of top absorbing boundary

    enddo

  endif  ! end of absorbing boundaries


! --- add the source
  if(.not. initialfield) then

! collocated force
! beware, for acoustic medium, source is a potential, therefore source time function
! gives shape of velocity, not displacement
  if(source_type == 1) then
    if(elastic(ispec_selected_source)) then
      accel_elastic(1,iglob_source) = accel_elastic(1,iglob_source) - sin(angleforce)*source_time_function(it)
      accel_elastic(2,iglob_source) = accel_elastic(2,iglob_source) + cos(angleforce)*source_time_function(it)
    endif

! moment tensor
  else if(source_type == 2) then

    if(elastic(ispec_selected_source)) then
! add source array
      do j=1,NGLLZ
        do i=1,NGLLX
          iglob = ibool(i,j,ispec_selected_source)
          accel_elastic(:,iglob) = accel_elastic(:,iglob) + sourcearray(:,i,j)*source_time_function(it)
        enddo
      enddo
    endif

  endif

  else
    stop 'wrong source type'
  endif

! implement attenuation
  if(TURN_ATTENUATION_ON) then

! compute Grad(displ_elastic) at time step n+1 for attenuation
    call compute_gradient_attenuation(displ_elastic,dux_dxl_np1,duz_dxl_np1, &
      dux_dzl_np1,duz_dzl_np1,xix,xiz,gammax,gammaz,ibool,elastic,hprime_xx,hprime_zz,nspec,npoin)

! update memory variables with fourth-order Runge-Kutta time scheme for attenuation
! loop over spectral elements
  do ispec = 1,nspec

  do j=1,NGLLZ
  do i=1,NGLLX

  theta_n   = dux_dxl_n(i,j,ispec) + duz_dzl_n(i,j,ispec)
  theta_np1 = dux_dxl_np1(i,j,ispec) + duz_dzl_np1(i,j,ispec)

! evolution e1_mech1
  Un = e1_mech1(i,j,ispec)
  tauinv = - inv_tau_sigma_nu1_mech1
  tauinvsquare = tauinv * tauinv
  tauinvcube = tauinvsquare * tauinv
  tauinvUn = tauinv * Un
  Sn   = theta_n * phi_nu1_mech1
  Snp1 = theta_np1 * phi_nu1_mech1
  Unp1 = Un + (deltatfourth*tauinvcube*(Sn + tauinvUn) + &
      twelvedeltat*(Sn + Snp1 + 2*tauinvUn) + &
      fourdeltatsquare*tauinv*(2*Sn + Snp1 + 3*tauinvUn) + &
      deltatcube*tauinvsquare*(3*Sn + Snp1 + 4*tauinvUn))* ONE_OVER_24
  e1_mech1(i,j,ispec) = Unp1

! evolution e1_mech2
  Un = e1_mech2(i,j,ispec)
  tauinv = - inv_tau_sigma_nu1_mech2
  tauinvsquare = tauinv * tauinv
  tauinvcube = tauinvsquare * tauinv
  tauinvUn = tauinv * Un
  Sn   = theta_n * phi_nu1_mech2
  Snp1 = theta_np1 * phi_nu1_mech2
  Unp1 = Un + (deltatfourth*tauinvcube*(Sn + tauinvUn) + &
      twelvedeltat*(Sn + Snp1 + 2*tauinvUn) + &
      fourdeltatsquare*tauinv*(2*Sn + Snp1 + 3*tauinvUn) + &
      deltatcube*tauinvsquare*(3*Sn + Snp1 + 4*tauinvUn))* ONE_OVER_24
  e1_mech2(i,j,ispec) = Unp1

! evolution e11_mech1
  Un = e11_mech1(i,j,ispec)
  tauinv = - inv_tau_sigma_nu2_mech1
  tauinvsquare = tauinv * tauinv
  tauinvcube = tauinvsquare * tauinv
  tauinvUn = tauinv * Un
  Sn   = (dux_dxl_n(i,j,ispec) - theta_n/TWO) * phi_nu2_mech1
  Snp1 = (dux_dxl_np1(i,j,ispec) - theta_np1/TWO) * phi_nu2_mech1
  Unp1 = Un + (deltatfourth*tauinvcube*(Sn + tauinvUn) + &
      twelvedeltat*(Sn + Snp1 + 2*tauinvUn) + &
      fourdeltatsquare*tauinv*(2*Sn + Snp1 + 3*tauinvUn) + &
      deltatcube*tauinvsquare*(3*Sn + Snp1 + 4*tauinvUn))* ONE_OVER_24
  e11_mech1(i,j,ispec) = Unp1

! evolution e11_mech2
  Un = e11_mech2(i,j,ispec)
  tauinv = - inv_tau_sigma_nu2_mech2
  tauinvsquare = tauinv * tauinv
  tauinvcube = tauinvsquare * tauinv
  tauinvUn = tauinv * Un
  Sn   = (dux_dxl_n(i,j,ispec) - theta_n/TWO) * phi_nu2_mech2
  Snp1 = (dux_dxl_np1(i,j,ispec) - theta_np1/TWO) * phi_nu2_mech2
  Unp1 = Un + (deltatfourth*tauinvcube*(Sn + tauinvUn) + &
      twelvedeltat*(Sn + Snp1 + 2*tauinvUn) + &
      fourdeltatsquare*tauinv*(2*Sn + Snp1 + 3*tauinvUn) + &
      deltatcube*tauinvsquare*(3*Sn + Snp1 + 4*tauinvUn))* ONE_OVER_24
  e11_mech2(i,j,ispec) = Unp1

! evolution e13_mech1
  Un = e13_mech1(i,j,ispec)
  tauinv = - inv_tau_sigma_nu2_mech1
  tauinvsquare = tauinv * tauinv
  tauinvcube = tauinvsquare * tauinv
  tauinvUn = tauinv * Un
  Sn   = (dux_dzl_n(i,j,ispec) + duz_dxl_n(i,j,ispec)) * phi_nu2_mech1
  Snp1 = (dux_dzl_np1(i,j,ispec) + duz_dxl_np1(i,j,ispec)) * phi_nu2_mech1
  Unp1 = Un + (deltatfourth*tauinvcube*(Sn + tauinvUn) + &
      twelvedeltat*(Sn + Snp1 + 2*tauinvUn) + &
      fourdeltatsquare*tauinv*(2*Sn + Snp1 + 3*tauinvUn) + &
      deltatcube*tauinvsquare*(3*Sn + Snp1 + 4*tauinvUn))* ONE_OVER_24
  e13_mech1(i,j,ispec) = Unp1

! evolution e13_mech2
  Un = e13_mech2(i,j,ispec)
  tauinv = - inv_tau_sigma_nu2_mech2
  tauinvsquare = tauinv * tauinv
  tauinvcube = tauinvsquare * tauinv
  tauinvUn = tauinv * Un
  Sn   = (dux_dzl_n(i,j,ispec) + duz_dxl_n(i,j,ispec)) * phi_nu2_mech2
  Snp1 = (dux_dzl_np1(i,j,ispec) + duz_dxl_np1(i,j,ispec)) * phi_nu2_mech2
  Unp1 = Un + (deltatfourth*tauinvcube*(Sn + tauinvUn) + &
      twelvedeltat*(Sn + Snp1 + 2*tauinvUn) + &
      fourdeltatsquare*tauinv*(2*Sn + Snp1 + 3*tauinvUn) + &
      deltatcube*tauinvsquare*(3*Sn + Snp1 + 4*tauinvUn))* ONE_OVER_24
  e13_mech2(i,j,ispec) = Unp1

  enddo
  enddo
  enddo

  endif ! end of test on attenuation

  end subroutine compute_forces_elastic

