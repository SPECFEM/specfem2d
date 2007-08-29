
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
               ispec_selected_source,is_proc_source,source_type,it,NSTEP,anyabs,assign_external_model, &
               initialfield,TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON,angleforce,deltatcube, &
               deltatfourth,twelvedeltat,fourdeltatsquare,ibool,kmato,numabs,elastic,codeabs, &
               accel_elastic,veloc_elastic,displ_elastic,density,elastcoef,xix,xiz,gammax,gammaz, &
               jacobian,vpext,vsext,rhoext,source_time_function,sourcearray,e1,e11, &
               e13,dux_dxl_n,duz_dzl_n,duz_dxl_n,dux_dzl_n, &
               dux_dxl_np1,duz_dzl_np1,duz_dxl_np1,dux_dzl_np1,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll,inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2,Mu_nu1,Mu_nu2)

! compute forces for the elastic elements

  implicit none

  include "constants.h"

  integer :: npoin,nspec,nelemabs,numat,iglob_source,ispec_selected_source,is_proc_source,source_type,it,NSTEP

  logical :: anyabs,assign_external_model,initialfield,TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON

  double precision :: angleforce,deltatcube,deltatfourth,twelvedeltat,fourdeltatsquare

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: kmato
  integer, dimension(nelemabs) :: numabs

  logical, dimension(nspec) :: elastic
  logical, dimension(4,nelemabs)  :: codeabs

  real(kind=CUSTOM_REAL), dimension(NDIM,npoin) :: accel_elastic,veloc_elastic,displ_elastic
  double precision, dimension(numat) :: density
  double precision, dimension(4,numat) :: elastcoef
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz,jacobian
  double precision, dimension(NGLLX,NGLLZ,nspec) :: vpext,vsext,rhoext
  real(kind=CUSTOM_REAL), dimension(NSTEP) :: source_time_function
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: sourcearray

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec,N_SLS) :: e1,e11,e13
  double precision, dimension(N_SLS) :: inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2
  double precision :: Mu_nu1,Mu_nu2
  real(kind=CUSTOM_REAL) :: e1_sum,e11_sum,e13_sum
  integer :: i_sls

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: &
    dux_dxl_n,duz_dzl_n,duz_dxl_n,dux_dzl_n,dux_dxl_np1,duz_dzl_np1,duz_dxl_np1,dux_dzl_np1

! derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz

! Gauss-Lobatto-Legendre weights
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: wxgll
  real(kind=CUSTOM_REAL), dimension(NGLLZ) :: wzgll

!---
!--- local variables
!---

  integer :: ispec,i,j,k,iglob,ispecabs,ibegin,iend

! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_xz,sigma_zz
  real(kind=CUSTOM_REAL) :: nx,nz,vx,vz,vn,rho_vp,rho_vs,tx,tz,weight,xxi,zxi,xgamma,zgamma,jacobian1D

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2,tempz1,tempz2

! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

! material properties of the elastic medium
  real(kind=CUSTOM_REAL) :: mul_relaxed,lambdal_relaxed,lambdalplus2mul_relaxed,kappal,cpl,csl,rhol, &
                      lambdal_unrelaxed,mul_unrelaxed,lambdalplus2mul_unrelaxed

! for attenuation
  real(kind=CUSTOM_REAL) :: Un,Unp1,tauinv,Sn,Snp1,theta_n,theta_np1,tauinvsquare,tauinvcube,tauinvUn

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
    e1_sum = 0._CUSTOM_REAL
    e11_sum = 0._CUSTOM_REAL
    e13_sum = 0._CUSTOM_REAL

    do i_sls = 1,N_SLS
      e1_sum = e1_sum + e1(i,j,ispec,i_sls)
      e11_sum = e11_sum + e11(i,j,ispec,i_sls)
      e13_sum = e13_sum + e13(i,j,ispec,i_sls)
    enddo

    sigma_xx = sigma_xx + (lambdal_relaxed + mul_relaxed) * e1_sum + TWO * mul_relaxed * e11_sum
    sigma_xz = sigma_xz + mul_relaxed * e13_sum
    sigma_zz = sigma_zz + (lambdal_relaxed + mul_relaxed) * e1_sum - TWO * mul_relaxed * e11_sum

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
            accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tempx1(k,j)*hprimewgll_xx(k,i) + tempx2(i,k)*hprimewgll_zz(k,j))
            accel_elastic(2,iglob) = accel_elastic(2,iglob) - (tempz1(k,j)*hprimewgll_xx(k,i) + tempz2(i,k)*hprimewgll_zz(k,j))
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
      kappal  = lambdal_relaxed + TWO*mul_relaxed/3._CUSTOM_REAL
      cpl = sqrt((kappal + 4._CUSTOM_REAL*mul_relaxed/3._CUSTOM_REAL)/rhol)
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

     if ( is_proc_source == 1 ) then
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
     endif

  else
     call exit_MPI('wrong source type')
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

! loop on all the standard linear solids
  do i_sls = 1,N_SLS

! evolution e1
  Un = e1(i,j,ispec,i_sls)
  tauinv = - inv_tau_sigma_nu1(i_sls)
  tauinvsquare = tauinv * tauinv
  tauinvcube = tauinvsquare * tauinv
  tauinvUn = tauinv * Un
  Sn   = theta_n * phi_nu1(i_sls)
  Snp1 = theta_np1 * phi_nu1(i_sls)
  Unp1 = Un + (deltatfourth*tauinvcube*(Sn + tauinvUn) + &
      twelvedeltat*(Sn + Snp1 + 2*tauinvUn) + &
      fourdeltatsquare*tauinv*(2*Sn + Snp1 + 3*tauinvUn) + &
      deltatcube*tauinvsquare*(3*Sn + Snp1 + 4*tauinvUn))* ONE_OVER_24
  e1(i,j,ispec,i_sls) = Unp1

! evolution e11
  Un = e11(i,j,ispec,i_sls)
  tauinv = - inv_tau_sigma_nu2(i_sls)
  tauinvsquare = tauinv * tauinv
  tauinvcube = tauinvsquare * tauinv
  tauinvUn = tauinv * Un
  Sn   = (dux_dxl_n(i,j,ispec) - theta_n/TWO) * phi_nu2(i_sls)
  Snp1 = (dux_dxl_np1(i,j,ispec) - theta_np1/TWO) * phi_nu2(i_sls)
  Unp1 = Un + (deltatfourth*tauinvcube*(Sn + tauinvUn) + &
      twelvedeltat*(Sn + Snp1 + 2*tauinvUn) + &
      fourdeltatsquare*tauinv*(2*Sn + Snp1 + 3*tauinvUn) + &
      deltatcube*tauinvsquare*(3*Sn + Snp1 + 4*tauinvUn))* ONE_OVER_24
  e11(i,j,ispec,i_sls) = Unp1

! evolution e13
  Un = e13(i,j,ispec,i_sls)
  tauinv = - inv_tau_sigma_nu2(i_sls)
  tauinvsquare = tauinv * tauinv
  tauinvcube = tauinvsquare * tauinv
  tauinvUn = tauinv * Un
  Sn   = (dux_dzl_n(i,j,ispec) + duz_dxl_n(i,j,ispec)) * phi_nu2(i_sls)
  Snp1 = (dux_dzl_np1(i,j,ispec) + duz_dxl_np1(i,j,ispec)) * phi_nu2(i_sls)
  Unp1 = Un + (deltatfourth*tauinvcube*(Sn + tauinvUn) + &
      twelvedeltat*(Sn + Snp1 + 2*tauinvUn) + &
      fourdeltatsquare*tauinv*(2*Sn + Snp1 + 3*tauinvUn) + &
      deltatcube*tauinvsquare*(3*Sn + Snp1 + 4*tauinvUn))* ONE_OVER_24
  e13(i,j,ispec,i_sls) = Unp1

  enddo

  enddo
  enddo
  enddo

  endif ! end of test on attenuation

  end subroutine compute_forces_elastic

