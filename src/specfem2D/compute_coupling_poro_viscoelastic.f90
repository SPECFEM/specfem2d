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

! for poroelastic solver

  subroutine compute_coupling_poro_viscoelastic(displ_elastic,displs_poroelastic,displw_poroelastic,accels_poroelastic)

  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLZ,ZERO,TWO,IRIGHT,ILEFT,IBOTTOM,ITOP

  use specfem_par, only: num_solid_poro_edges, &
                         ibool,wxgll,wzgll,xix,xiz,gammax,gammaz,jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse, &
                         hprime_xx,hprime_zz, &
                         solid_poro_elastic_ispec,solid_poro_elastic_iedge, &
                         solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge, &
                         kmato,poroelastcoef, &
                         assign_external_model,c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,anisotropy, &
                         nglob_elastic,nglob_poroelastic

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_elastic), intent(in) :: displ_elastic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic), intent(in) :: displs_poroelastic,displw_poroelastic

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic), intent(inout) :: accels_poroelastic

  !local variables
  integer :: inum,ispec_elastic,iedge_elastic,ispec_poroelastic,iedge_poroelastic, &
             i,j,k,ipoin1D,iglob

  double precision :: kappa_s,kappa_f,kappa_fr,mu_s,mu_fr,rho_s,rho_f,eta_f,phi,tort,rho_bar
  double precision :: D_biot,H_biot,C_biot,M_biot

  double precision :: c11,c13,c15,c33,c35,c55,c12,c23,c25
  double precision :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic,lambdaplus2mu_unrelaxed_elastic

  real(kind=CUSTOM_REAL) :: mu_G,lambdal_G,lambdalplus2mul_G
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: dwx_dxi,dwx_dgamma,dwz_dxi,dwz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  real(kind=CUSTOM_REAL) :: dwx_dxl,dwz_dxl,dwx_dzl,dwz_dzl
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_xz,sigma_zz
  real(kind=CUSTOM_REAL) :: xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl
  real(kind=CUSTOM_REAL) :: sigmap

  ! loop on all the coupling edges
  do inum = 1,num_solid_poro_edges

    ! get the edge of the elastic element
    ispec_elastic = solid_poro_elastic_ispec(inum)
    iedge_elastic = solid_poro_elastic_iedge(inum)

    ! get the corresponding edge of the poroelastic element
    ispec_poroelastic = solid_poro_poroelastic_ispec(inum)
    iedge_poroelastic = solid_poro_poroelastic_iedge(inum)

    ! implement 1D coupling along the edge
    do ipoin1D = 1,NGLLX

      ! get point values for the elastic side, which matches our side in the inverse direction
      i = ivalue_inverse(ipoin1D,iedge_elastic)
      j = jvalue_inverse(ipoin1D,iedge_elastic)
      iglob = ibool(i,j,ispec_elastic)

      ! get elastic properties
      lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec_elastic))
      mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec_elastic))
      lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec_elastic))

      ! derivative along x and along z for u_s and w
      dux_dxi = ZERO
      duz_dxi = ZERO

      dux_dgamma = ZERO
      duz_dgamma = ZERO

      ! first double loop over GLL points to compute and store gradients
      ! we can merge the two loops because NGLLX == NGLLZ
      do k = 1,NGLLX
        dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec_elastic))*hprime_xx(i,k)
        duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec_elastic))*hprime_xx(i,k)
        dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec_elastic))*hprime_zz(j,k)
        duz_dgamma = duz_dgamma + displ_elastic(2,ibool(i,k,ispec_elastic))*hprime_zz(j,k)
      enddo

      xixl = xix(i,j,ispec_elastic)
      xizl = xiz(i,j,ispec_elastic)
      gammaxl = gammax(i,j,ispec_elastic)
      gammazl = gammaz(i,j,ispec_elastic)

      ! derivatives of displacement
      dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
      dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

      duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
      duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

      ! compute stress tensor
      ! full anisotropy
      if (kmato(ispec_elastic) == 2) then
        ! implement anisotropy in 2D
        if (assign_external_model) then
          c11 = c11ext(i,j,ispec_elastic)
          c13 = c13ext(i,j,ispec_elastic)
          c15 = c15ext(i,j,ispec_elastic)
          c33 = c33ext(i,j,ispec_elastic)
          c35 = c35ext(i,j,ispec_elastic)
          c55 = c55ext(i,j,ispec_elastic)
          c12 = c12ext(i,j,ispec_elastic)
          c23 = c23ext(i,j,ispec_elastic)
          c25 = c25ext(i,j,ispec_elastic)
        else
          c11 = anisotropy(1,kmato(ispec_elastic))
          c13 = anisotropy(2,kmato(ispec_elastic))
          c15 = anisotropy(3,kmato(ispec_elastic))
          c33 = anisotropy(4,kmato(ispec_elastic))
          c35 = anisotropy(5,kmato(ispec_elastic))
          c55 = anisotropy(6,kmato(ispec_elastic))
          c12 = anisotropy(7,kmato(ispec_elastic))
          c23 = anisotropy(8,kmato(ispec_elastic))
          c25 = anisotropy(9,kmato(ispec_elastic))
        endif
        sigma_xx = c11*dux_dxl + c15*(duz_dxl + dux_dzl) + c13*duz_dzl
        sigma_zz = c13*dux_dxl + c35*(duz_dxl + dux_dzl) + c33*duz_dzl
        sigma_xz = c15*dux_dxl + c55*(duz_dxl + dux_dzl) + c35*duz_dzl
      else
        ! no attenuation
        sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
        sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
        sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl
      endif

      ! get point values for the poroelastic side
      i = ivalue(ipoin1D,iedge_poroelastic)
      j = jvalue(ipoin1D,iedge_poroelastic)
      iglob = ibool(i,j,ispec_poroelastic)

      ! gets poroelastic material
      call get_poroelastic_material(ispec_poroelastic,phi,tort,mu_s,kappa_s,rho_s, &
                                    kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar)

      ! Biot coefficients for the input phi
      call get_poroelastic_Biot_coeff(phi,kappa_s,kappa_f,kappa_fr,mu_fr,D_biot,H_biot,C_biot,M_biot)

      mu_G = mu_fr
      lambdal_G = H_biot - 2._CUSTOM_REAL*mu_fr
      lambdalplus2mul_G = lambdal_G + TWO*mu_G

      ! derivative along x and along z for u_s and w
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
        dux_dxi = dux_dxi + displs_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
        duz_dxi = duz_dxi + displs_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
        dux_dgamma = dux_dgamma + displs_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
        duz_dgamma = duz_dgamma + displs_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)

        dwx_dxi = dwx_dxi + displw_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
        dwz_dxi = dwz_dxi + displw_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
        dwx_dgamma = dwx_dgamma + displw_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
        dwz_dgamma = dwz_dgamma + displw_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
      enddo

      xixl = xix(i,j,ispec_poroelastic)
      xizl = xiz(i,j,ispec_poroelastic)
      gammaxl = gammax(i,j,ispec_poroelastic)
      gammazl = gammaz(i,j,ispec_poroelastic)

      ! derivatives of displacement
      dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
      dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

      duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
      duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

      dwx_dxl = dwx_dxi*xixl + dwx_dgamma*gammaxl
      dwx_dzl = dwx_dxi*xizl + dwx_dgamma*gammazl

      dwz_dxl = dwz_dxi*xixl + dwz_dgamma*gammaxl
      dwz_dzl = dwz_dxi*xizl + dwz_dgamma*gammazl

      ! compute stress tensor

      ! no attenuation
      sigma_xx = sigma_xx + lambdalplus2mul_G*dux_dxl + lambdal_G*duz_dzl + C_biot*(dwx_dxl + dwz_dzl)
      sigma_xz = sigma_xz + mu_G*(duz_dxl + dux_dzl)
      sigma_zz = sigma_zz + lambdalplus2mul_G*duz_dzl + lambdal_G*dux_dxl + C_biot*(dwx_dxl + dwz_dzl)

      sigmap = C_biot*(dux_dxl + duz_dzl) + M_biot*(dwx_dxl + dwz_dzl)

      ! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
      ! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
      ! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
      ! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
      ! Blackwell Science, page 110, equation (4.60).
      if (iedge_poroelastic == ITOP) then
        xxi = + gammaz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        zxi = - gammax(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D
        weight = jacobian1D * wxgll(i)
      else if (iedge_poroelastic == IBOTTOM) then
        xxi = + gammaz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        zxi = - gammax(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D
        weight = jacobian1D * wxgll(i)
      else if (iedge_poroelastic == ILEFT) then
        xgamma = - xiz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        zgamma = + xix(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D
        weight = jacobian1D * wzgll(j)
      else if (iedge_poroelastic == IRIGHT) then
        xgamma = - xiz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        zgamma = + xix(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D
        weight = jacobian1D * wzgll(j)
      endif

      ! contribution to the solid phase
      accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) + &
        weight*((sigma_xx*nx + sigma_xz*nz)/2.d0 -phi/tort*sigmap*nx)

      accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) + &
        weight*((sigma_xz*nx + sigma_zz*nz)/2.d0 -phi/tort*sigmap*nz)

      ! contribution to the fluid phase
      ! w = 0
    enddo

  enddo

  end subroutine compute_coupling_poro_viscoelastic

!
!========================================================================
!

  subroutine compute_coupling_poro_viscoelastic_for_stabilization(veloc_elastic,accel_elastic, &
                                                                  velocs_poroelastic,accels_poroelastic, &
                                                                  velocw_poroelastic,accelw_poroelastic,deltatover2)

! Explanation of the code below, from Christina Morency and Yang Luo, January 2012:
!
! Coupled elastic-poroelastic simulations imply continuity of traction and
! displacement at the interface.
! For the traction we pass on both sides n*(T + Te)/2 , that is, the average
! between the total stress (from the poroelastic part) and the elastic stress.
! For the displacement, we enforce its continuity in the assembling stage,
! realizing that continuity of displacement correspond to the continuity of
! the acceleration we have:
!
! accel_elastic = rmass_inverse_elastic * force_elastic
! accels_poroelastic = rmass_s_inverse_poroelastic * force_poroelastic
!
! Therefore, continuity of acceleration gives
!
! accel = (force_elastic + force_poroelastic)/
!     (1/rmass_inverse_elastic + 1/rmass_inverse_poroelastic)
!
! Then
!
! accel_elastic = accel
! accels_poroelastic = accel
! accelw_poroelastic = 0
!
! From there, the velocity and displacement are updated.
! Note that force_elastic and force_poroelastic are the right hand sides of
! the equations we solve, that is, the acceleration terms before the
! division by the inverse of the mass matrices. This is why in the code below
! we first need to recover the accelerations (which are then
! the right hand sides forces) and the velocities before the update.
!
! This implementation highly helped stability especially with unstructured meshes.

  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLZ,ZERO

  use specfem_par, only: num_solid_poro_edges, &
                         ibool,ivalue,jvalue, &
                         solid_poro_elastic_ispec,solid_poro_elastic_iedge, &
                         solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge, &
                         rmass_inverse_elastic, &
                         rmass_s_inverse_poroelastic, &
                         time_stepping_scheme,nglob, &
                         nglob_elastic,nglob_poroelastic


  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_elastic), intent(inout) :: veloc_elastic,accel_elastic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic), intent(inout) :: velocs_poroelastic,accels_poroelastic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic), intent(inout) :: velocw_poroelastic,accelw_poroelastic

  double precision :: deltatover2

  !local variables
  integer :: inum,ispec_elastic,iedge_elastic,ispec_poroelastic,iedge_poroelastic, &
             i,j,ipoin1D,iglob
  logical,dimension(nglob) :: mask_ibool

  ! initializes
  mask_ibool(:) = .false.

  ! loop on all the coupling edges
  do inum = 1,num_solid_poro_edges
    ! get the edge of the elastic element
    ispec_elastic = solid_poro_elastic_ispec(inum)
    iedge_elastic = solid_poro_elastic_iedge(inum)

    ! get the corresponding edge of the poroelastic element
    ispec_poroelastic = solid_poro_poroelastic_ispec(inum)
    iedge_poroelastic = solid_poro_poroelastic_iedge(inum)

    do ipoin1D = 1,NGLLX
      ! recovering original velocities and accelerations on boundaries (elastic side)
      i = ivalue(ipoin1D,iedge_poroelastic)
      j = jvalue(ipoin1D,iedge_poroelastic)

      ! gets global boundary node
      iglob = ibool(i,j,ispec_poroelastic)

      ! stabilization imposes continuity
      if (.not. mask_ibool(iglob)) then
        ! only do this once on a global node
        mask_ibool(iglob) = .true.

        if (time_stepping_scheme == 1) then
          veloc_elastic(1,iglob) = veloc_elastic(1,iglob) - deltatover2*accel_elastic(1,iglob)
          veloc_elastic(2,iglob) = veloc_elastic(2,iglob) - deltatover2*accel_elastic(2,iglob)

          accel_elastic(1,iglob) = accel_elastic(1,iglob) / rmass_inverse_elastic(1,iglob)
          accel_elastic(2,iglob) = accel_elastic(2,iglob) / rmass_inverse_elastic(2,iglob)

          ! recovering original velocities and accelerations on boundaries (poro side)
          velocs_poroelastic(1,iglob) = velocs_poroelastic(1,iglob) - deltatover2*accels_poroelastic(1,iglob)
          velocs_poroelastic(2,iglob) = velocs_poroelastic(2,iglob) - deltatover2*accels_poroelastic(2,iglob)

          accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) / rmass_s_inverse_poroelastic(iglob)
          accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) / rmass_s_inverse_poroelastic(iglob)

          ! assembling accelerations
          accel_elastic(1,iglob) = ( accel_elastic(1,iglob) + accels_poroelastic(1,iglob) ) / &
                                  ( 1.0/rmass_inverse_elastic(1,iglob) + 1.0/rmass_s_inverse_poroelastic(iglob) )
          accel_elastic(2,iglob) = ( accel_elastic(2,iglob) + accels_poroelastic(2,iglob) ) / &
                                  ( 1.0/rmass_inverse_elastic(2,iglob) + 1.0/rmass_s_inverse_poroelastic(iglob) )

          ! imposes continuity
          accels_poroelastic(1,iglob) = accel_elastic(1,iglob)
          accels_poroelastic(2,iglob) = accel_elastic(2,iglob)

          ! updating velocities
          velocs_poroelastic(1,iglob) = velocs_poroelastic(1,iglob) + deltatover2*accels_poroelastic(1,iglob)
          velocs_poroelastic(2,iglob) = velocs_poroelastic(2,iglob) + deltatover2*accels_poroelastic(2,iglob)

          veloc_elastic(1,iglob) = veloc_elastic(1,iglob) + deltatover2*accel_elastic(1,iglob)
          veloc_elastic(2,iglob) = veloc_elastic(2,iglob) + deltatover2*accel_elastic(2,iglob)

          ! zeros w
          accelw_poroelastic(1,iglob) = 0._CUSTOM_REAL
          accelw_poroelastic(2,iglob) = 0._CUSTOM_REAL
          velocw_poroelastic(1,iglob) = 0._CUSTOM_REAL
          velocw_poroelastic(2,iglob) = 0._CUSTOM_REAL
        endif

!         if (time_stepping_scheme == 2) then
!        recovering original velocities and accelerations on boundaries (elastic side)
!      veloc_elastic = veloc_elastic - BETA_LDDRK(i_stage) * veloc_elastic_LDDRK
!      displ_elastic = displ_elastic - BETA_LDDRK(i_stage) * displ_elastic_LDDRK
!      veloc_elastic_LDDRK = (veloc_elastic_LDDRK - deltat * accel_elastic) / ALPHA_LDDRK(i_stage)
!      displ_elastic_LDDRK = (displ_elastic_LDDRK - deltat * veloc_elastic) / ALPHA_LDDRK(i_stage)
!            accel_elastic(1,iglob) = accel_elastic(1,iglob) / rmass_inverse_elastic(1,iglob)
!            accel_elastic(2,iglob) = accel_elastic(2,iglob) / rmass_inverse_elastic(2,iglob)

            ! recovering original velocities and accelerations on boundaries (poro side)
!      velocs_poroelastic = velocs_poroelastic - BETA_LDDRK(i_stage) * velocs_poroelastic_LDDRK
!      displs_poroelastic = displs_poroelastic - BETA_LDDRK(i_stage) * displs_poroelastic_LDDRK
!      velocs_poroelastic_LDDRK = (velocs_poroelastic_LDDRK - deltat * accels_poroelastic) / ALPHA_LDDRK(i_stage)
!      displs_poroelastic_LDDRK = (velocs_poroelastic_LDDRK - deltat * velocs_poroelastic) / ALPHA_LDDRK(i_stage)
!            accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) / rmass_s_inverse_poroelastic(iglob)
!            accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) / rmass_s_inverse_poroelastic(iglob)

            ! assembling accelerations
!            accel_elastic(1,iglob) = ( accel_elastic(1,iglob) + accels_poroelastic(1,iglob) ) / &
!                                   ( 1.0/rmass_inverse_elastic(1,iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
!            accel_elastic(2,iglob) = ( accel_elastic(2,iglob) + accels_poroelastic(2,iglob) ) / &
!                                   ( 1.0/rmass_inverse_elastic(2,iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
!            accels_poroelastic(1,iglob) = accel_elastic(1,iglob)
!            accels_poroelastic(2,iglob) = accel_elastic(2,iglob)

      ! updating velocities
            ! updating velocities(elastic side)
!      veloc_elastic_LDDRK = ALPHA_LDDRK(i_stage) * veloc_elastic_LDDRK + deltat * accel_elastic
!      displ_elastic_LDDRK = ALPHA_LDDRK(i_stage) * displ_elastic_LDDRK + deltat * veloc_elastic
!      veloc_elastic = veloc_elastic + BETA_LDDRK(i_stage) * veloc_elastic_LDDRK
!      displ_elastic = displ_elastic + BETA_LDDRK(i_stage) * displ_elastic_LDDRK
            ! updating velocities(poro side)
!      velocs_poroelastic_LDDRK = ALPHA_LDDRK(i_stage) * velocs_poroelastic_LDDRK + deltat * accels_poroelastic
!      displs_poroelastic_LDDRK = ALPHA_LDDRK(i_stage) * displs_poroelastic_LDDRK + deltat * velocs_poroelastic
!      velocs_poroelastic = velocs_poroelastic + BETA_LDDRK(i_stage) * velocs_poroelastic_LDDRK
!      displs_poroelastic = displs_poroelastic + BETA_LDDRK(i_stage) * displs_poroelastic_LDDRK

            ! zeros w
!            accelw_poroelastic(1,iglob) = ZERO
!            accelw_poroelastic(2,iglob) = ZERO
!            velocw_poroelastic(1,iglob) = ZERO
!            velocw_poroelastic(2,iglob) = ZERO
!            endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      if (time_stepping_scheme == 3) then

        ! recovering original velocities and accelerations on boundaries (elastic side)
!        if (i_stage==1 .or. i_stage==2 .or. i_stage==3) then

!        if (i_stage == 1)weight_rk = 0.5d0
!        if (i_stage == 2)weight_rk = 0.5d0
!        if (i_stage == 3)weight_rk = 1.0d0

!        veloc_elastic(1,iglob) = veloc_elastic_initial_rk(1,iglob) - weight_rk * accel_elastic_rk(1,iglob,i_stage)
!        veloc_elastic(2,iglob) = veloc_elastic_initial_rk(2,iglob) - weight_rk * accel_elastic_rk(2,iglob,i_stage)
!        displ_elastic(1,iglob) = displ_elastic_initial_rk(1,iglob) - weight_rk * veloc_elastic_rk(1,iglob,i_stage)
!        displ_elastic(2,iglob) = displ_elastic_initial_rk(2,iglob) - weight_rk * veloc_elastic_rk(2,iglob,i_stage)


!        else if (i_stage==4) then

!        veloc_elastic(1,iglob) = veloc_elastic_initial_rk(1,iglob) - 1.0d0 / 6.0d0 * &
!        (accel_elastic_rk(1,iglob,1) + 2.0d0 * accel_elastic_rk(1,iglob,2) + &
!         2.0d0 * accel_elastic_rk(1,iglob,3) + accel_elastic_rk(1,iglob,4))

!        veloc_elastic(2,iglob) = veloc_elastic_initial_rk(2,iglob) - 1.0d0 / 6.0d0 * &
!        (accel_elastic_rk(2,iglob,1) + 2.0d0 * accel_elastic_rk(2,iglob,2) + &
!         2.0d0 * accel_elastic_rk(2,iglob,3) + accel_elastic_rk(2,iglob,4))

!        displ_elastic(1,iglob) = displ_elastic_initial_rk(1,iglob) - 1.0d0 / 6.0d0 * &
!        (veloc_elastic_rk(1,iglob,1) + 2.0d0 * veloc_elastic_rk(1,iglob,2) + &
!         2.0d0 * veloc_elastic_rk(1,iglob,3) + veloc_elastic_rk(1,iglob,4))

!        displ_elastic(2,iglob) = displ_elastic_initial_rk(2,iglob) - 1.0d0 / 6.0d0 * &
!        (veloc_elastic_rk(2,iglob,1) + 2.0d0 * veloc_elastic_rk(2,iglob,2) + &
!         2.0d0 * veloc_elastic_rk(2,iglob,3) + veloc_elastic_rk(2,iglob,4))

!        endif

!        accel_elastic(1,iglob) = accel_elastic(1,iglob) / rmass_inverse_elastic(1,iglob)
!        accel_elastic(2,iglob) = accel_elastic(2,iglob) / rmass_inverse_elastic(2,iglob)

!        accel_elastic_rk(1,iglob,i_stage) = accel_elastic(1,iglob) / deltat
!        accel_elastic_rk(2,iglob,i_stage) = accel_elastic(2,iglob) / deltat
!        veloc_elastic_rk(1,iglob,i_stage) = veloc_elastic(1,iglob) / deltat
!        veloc_elastic_rk(2,iglob,i_stage) = veloc_elastic(2,iglob) / deltat


        ! recovering original velocities and accelerations on boundaries (poro side)
!        if (i_stage==1 .or. i_stage==2 .or. i_stage==3) then

!        if (i_stage == 1)weight_rk = 0.5d0
!        if (i_stage == 2)weight_rk = 0.5d0
!        if (i_stage == 3)weight_rk = 1.0d0

!        velocs_poroelastic(1,iglob) = velocs_poroelastic_initial_rk(1,iglob) - weight_rk * accels_poroelastic_rk(1,iglob,i_stage)
!  velocs_poroelastic(2,iglob) = velocs_poroelastic_initial_rk(2,iglob) - weight_rk * accels_poroelastic_rk(2,iglob,i_stage)
!        displs_poroelastic(1,iglob) = displs_poroelastic_initial_rk(1,iglob) - weight_rk * velocs_poroelastic_rk(1,iglob,i_stage)
!  displs_poroelastic(2,iglob) = displs_poroelastic_initial_rk(2,iglob) - weight_rk * velocs_poroelastic_rk(2,iglob,i_stage)


!        else if (i_stage==4) then

!        velocs_poroelastic(1,iglob) = velocs_poroelastic_initial_rk(1,iglob) - 1.0d0 / 6.0d0 * &
!        (accels_poroelastic_rk(1,iglob,1) + 2.0d0 * accels_poroelastic_rk(1,iglob,2) + &
!         2.0d0 * accels_poroelastic_rk(1,iglob,3) + accels_poroelastic_rk(1,iglob,4))

!        velocs_poroelastic(2,iglob) = velocs_poroelastic_initial_rk(2,iglob) - 1.0d0 / 6.0d0 * &
!        (accels_poroelastic_rk(2,iglob,1) + 2.0d0 * accels_poroelastic_rk(2,iglob,2) + &
!         2.0d0 * accels_poroelastic_rk(2,iglob,3) + accels_poroelastic_rk(2,iglob,4))

!        displs_poroelastic(1,iglob) = displs_poroelastic_initial_rk(1,iglob) - 1.0d0 / 6.0d0 * &
!        (velocs_poroelastic_rk(1,iglob,1) + 2.0d0 * velocs_poroelastic_rk(1,iglob,2) + &
!         2.0d0 * velocs_poroelastic_rk(1,iglob,3) + velocs_poroelastic_rk(1,iglob,4))

!        displs_poroelastic(2,iglob) = displs_poroelastic_initial_rk(2,iglob) - 1.0d0 / 6.0d0 * &
!        (velocs_poroelastic_rk(2,iglob,1) + 2.0d0 * velocs_poroelastic_rk(2,iglob,2) + &
!         2.0d0 * velocs_poroelastic_rk(2,iglob,3) + velocs_poroelastic_rk(2,iglob,4))

!        endif

!        accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) / rmass_s_inverse_poroelastic(iglob)
!        accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) / rmass_s_inverse_poroelastic(iglob)

!        accels_poroelastic_rk(1,iglob,i_stage) = accels_poroelastic(1,iglob) / deltat
!        accels_poroelastic_rk(2,iglob,i_stage) = accels_poroelastic(2,iglob) / deltat
!        velocs_poroelastic_rk(1,iglob,i_stage) = velocs_poroelastic(1,iglob) / deltat
!        velocs_poroelastic_rk(2,iglob,i_stage) = velocs_poroelastic(2,iglob) / deltat


        ! assembling accelerations
!            accel_elastic(1,iglob) = ( accel_elastic(1,iglob) + accels_poroelastic(1,iglob) ) / &
!                                   ( 1.0/rmass_inverse_elastic(1,iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
!            accel_elastic(2,iglob) = ( accel_elastic(2,iglob) + accels_poroelastic(2,iglob) ) / &
!                                   ( 1.0/rmass_inverse_elastic(2,iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
!            accels_poroelastic(1,iglob) = accel_elastic(1,iglob)
!            accels_poroelastic(2,iglob) = accel_elastic(2,iglob)

   ! updating velocities
        ! updating velocities(elastic side)

 !       accel_elastic_rk(1,iglob,i_stage) = accel_elastic(1,iglob) * deltat
 !       accel_elastic_rk(2,iglob,i_stage) = accel_elastic(2,iglob) * deltat

 !       if (i_stage==1 .or. i_stage==2 .or. i_stage==3) then

 !       if (i_stage == 1)weight_rk = 0.5d0
 !       if (i_stage == 2)weight_rk = 0.5d0
 !       if (i_stage == 3)weight_rk = 1.0d0

 !       veloc_elastic(1,iglob) = veloc_elastic_initial_rk(1,iglob) + weight_rk * accel_elastic_rk(1,iglob,i_stage)
 !       veloc_elastic(2,iglob) = veloc_elastic_initial_rk(2,iglob) + weight_rk * accel_elastic_rk(2,iglob,i_stage)
 !       displ_elastic(1,iglob) = displ_elastic_initial_rk(1,iglob) + weight_rk * veloc_elastic_rk(1,iglob,i_stage)
 !       displ_elastic(2,iglob) = displ_elastic_initial_rk(2,iglob) + weight_rk * veloc_elastic_rk(2,iglob,i_stage)


 !       else if (i_stage==4) then

 !       veloc_elastic(1,iglob) = veloc_elastic_initial_rk(1,iglob) + 1.0d0 / 6.0d0 * &
 !       (accel_elastic_rk(1,iglob,1) + 2.0d0 * accel_elastic_rk(1,iglob,2) + &
 !        2.0d0 * accel_elastic_rk(1,iglob,3) + accel_elastic_rk(1,iglob,4))
!
 !       veloc_elastic(2,iglob) = veloc_elastic_initial_rk(2,iglob) + 1.0d0 / 6.0d0 * &
 !       (accel_elastic_rk(2,iglob,1) + 2.0d0 * accel_elastic_rk(2,iglob,2) + &
 !        2.0d0 * accel_elastic_rk(2,iglob,3) + accel_elastic_rk(2,iglob,4))

 !       displ_elastic(1,iglob) = displ_elastic_initial_rk(1,iglob) + 1.0d0 / 6.0d0 * &
 !       (veloc_elastic_rk(1,iglob,1) + 2.0d0 * veloc_elastic_rk(1,iglob,2) + &
 !        2.0d0 * veloc_elastic_rk(1,iglob,3) + veloc_elastic_rk(1,iglob,4))

 !       displ_elastic(2,iglob) = displ_elastic_initial_rk(2,iglob) + 1.0d0 / 6.0d0 * &
 !       (veloc_elastic_rk(2,iglob,1) + 2.0d0 * veloc_elastic_rk(2,iglob,2) + &
 !        2.0d0 * veloc_elastic_rk(2,iglob,3) + veloc_elastic_rk(2,iglob,4))

 !       endif
        ! updating velocities(poro side)

 !       accels_poroelastic_rk(1,iglob,i_stage) = deltat * accels_poroelastic(1,iglob)
 !       accels_poroelastic_rk(2,iglob,i_stage) = deltat * accels_poroelastic(2,iglob)
 !       velocs_poroelastic_rk(1,iglob,i_stage) = deltat * velocs_poroelastic(1,iglob)
 !       velocs_poroelastic_rk(2,iglob,i_stage) = deltat * velocs_poroelastic(2,iglob)


 !       if (i_stage==1 .or. i_stage==2 .or. i_stage==3) then

 !       if (i_stage == 1)weight_rk = 0.5d0
 !       if (i_stage == 2)weight_rk = 0.5d0
 !       if (i_stage == 3)weight_rk = 1.0d0

 !       velocs_poroelastic(1,iglob) = velocs_poroelastic_initial_rk(1,iglob) + weight_rk * accels_poroelastic_rk(1,iglob,i_stage)
 ! velocs_poroelastic(2,iglob) = velocs_poroelastic_initial_rk(2,iglob) + weight_rk * accels_poroelastic_rk(2,iglob,i_stage)
 !       displs_poroelastic(1,iglob) = displs_poroelastic_initial_rk(1,iglob) + weight_rk * velocs_poroelastic_rk(1,iglob,i_stage)
 ! displs_poroelastic(2,iglob) = displs_poroelastic_initial_rk(2,iglob) + weight_rk * velocs_poroelastic_rk(2,iglob,i_stage)


 !       else if (i_stage==4) then

 !       velocs_poroelastic(1,iglob) = velocs_poroelastic_initial_rk(1,iglob) + 1.0d0 / 6.0d0 * &
 !       (accels_poroelastic_rk(1,iglob,1) + 2.0d0 * accels_poroelastic_rk(1,iglob,2) + &
 !        2.0d0 * accels_poroelastic_rk(1,iglob,3) + accels_poroelastic_rk(1,iglob,4))

 !       velocs_poroelastic(2,iglob) = velocs_poroelastic_initial_rk(2,iglob) + 1.0d0 / 6.0d0 * &
 !       (accels_poroelastic_rk(2,iglob,1) + 2.0d0 * accels_poroelastic_rk(2,iglob,2) + &
 !        2.0d0 * accels_poroelastic_rk(2,iglob,3) + accels_poroelastic_rk(2,iglob,4))
!
 !       displs_poroelastic(1,iglob) = displs_poroelastic_initial_rk(1,iglob) + 1.0d0 / 6.0d0 * &
 !       (velocs_poroelastic_rk(1,iglob,1) + 2.0d0 * velocs_poroelastic_rk(1,iglob,2) + &
 !        2.0d0 * velocs_poroelastic_rk(1,iglob,3) + velocs_poroelastic_rk(1,iglob,4))
!
 !       displs_poroelastic(2,iglob) = displs_poroelastic_initial_rk(2,iglob) + 1.0d0 / 6.0d0 * &
 !       (velocs_poroelastic_rk(2,iglob,1) + 2.0d0 * velocs_poroelastic_rk(2,iglob,2) + &
 !        2.0d0 * velocs_poroelastic_rk(2,iglob,3) + velocs_poroelastic_rk(2,iglob,4))

 !       endif

 !     endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      endif ! mask
    enddo
  enddo

  end subroutine compute_coupling_poro_viscoelastic_for_stabilization


